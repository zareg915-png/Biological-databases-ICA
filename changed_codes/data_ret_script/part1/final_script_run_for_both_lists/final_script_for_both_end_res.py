#!/usr/bin/env python3
# =============================================================================
# GenAI declaration: Tool: Claude (Anthropic)
# Use: review and corrected comments
# =============================================================================
"""Part 1 (p1_basic_anno): resolve input gene labels to current human NCBI Gene records.

For each dataset (development, challenge) every input label (symbol, alias or
human Ensembl gene ID) is searched in NCBI Gene via E-utilities, checked to be a
current human gene, and annotated with its NCBI summary fields and Ensembl
(ENSG) cross-references. One row is written per input, in input order, so
unresolved, ambiguous and repeated inputs are kept and explained rather than
dropped.

Output: <project root>/OP1_data/part1_data/final_data/<dataset>/ncbi_gene_summary2.tsv
(plain, unquoted TSV; columns as in COLUMNS below).

This module also provides the shared HTTP, NCBI and file helpers and the
run() entry point that Part 2 imports, so keep its file name p1_basic_anno.py.

Exit codes: 0 = all outputs written cleanly; 2 = outputs written but some
lookups need attention (see Resolution_Status / Resolution_Note);
1 = a fatal error stopped the run.

Optional environment variables: NCBI_EMAIL and NCBI_API_KEY (sent to NCBI as
recommended in the E-utilities usage policy) and BDB_PROJECT_ROOT.
"""
import argparse
import csv
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import timezone
from email.utils import parsedate_to_datetime
from functools import cache
from http.client import HTTPException
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

# Input gene lists from the assignment brief. 
GENES = {
    "development": "BRCA1 BRCA2 TP53 EGFR MYC APOE CFTR ACE2 TNF VEGFA".split(),
    "challenge": "MYC c-Myc Blah ENSG00000136997 AADACL2 Plet1 HYDIN ENSG00000261846 XIST MALAT1".split(),
}
# Output column order. CHANGED in review: reordered to the agreed pipeline
# order and Chromosome removed. The MySQL loader rejects any other header
# order, and the database schema/dump follows this order.
COLUMNS = """Query_Symbol Official_Symbol Official_Name Gene_ID ENSG_IDs Map_Location
             Other_Aliases Summary Resolution_Status Resolution_Note Same_Gene_Note""".split()


DATASETS = ("development", "challenge")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
# NCBI Datasets v2 gene endpoint. 
NCBI_GENES = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/"

# Exceptions that indicate a bad or unexpected API response for one gene. They are caught per input, so one bad lookup is recorded in that row instead of stopping the whole dataset. ServiceError is a ValueError.
API_ERRORS = (ValueError, KeyError, TypeError, IndexError, AttributeError, ET.ParseError)


class ServiceError(ValueError):
    """A remote service failed or returned a response this script cannot trust.

    `status` holds the HTTP status code when there is one.
    """

    def __init__(self, message, status=None):
        super().__init__(message)
        self.status = status


# HTTP helpers
def retry_delay(value, attempt):
    """Return seconds to wait before a retry.

    Honours a Retry-After header given either as seconds or as an HTTP date,
    but never waits less than the exponential back-off (1, 2, 4 s ...).
    """
    try:
        return max(2 ** attempt, int(value))
    except (TypeError, ValueError):
        try:
            when = parsedate_to_datetime(value)
            return max(2 ** attempt, when.astimezone(timezone.utc).timestamp() - time.time())
        except (TypeError, ValueError, OverflowError):
            return 2 ** attempt


def fetch(url, params=None, form=None):
    """GET (or POST when `form` is given) a URL and return (body bytes, headers).

    E-utilities requests are tagged with the tool name and, if set, the user's
    email and API key. Each request is preceded by a 0.4 s pause to stay under
    NCBI's limit of 3 requests per second without a key. Transient failures
    (timeouts, 429 and 5xx responses, network errors) are retried up to three
    attempts in total; POST requests are never retried. Any failure is raised
    as ServiceError with a short, readable message.
    """
    params = dict(params or {})
    headers = {"User-Agent": "BDB-ICA", "Accept": "application/json"}
    if url.startswith(EUTILS):
        params.update(tool="bdb_ica")
        for key, variable in (("email", "NCBI_EMAIL"), ("api_key", "NCBI_API_KEY")):
            if os.getenv(variable):
                params[key] = os.environ[variable]
        if params.get("retmode") == "xml":
            headers["Accept"] = "application/xml"
    if params:
        url += ("&" if "?" in url else "?") + urlencode(params)
    request = Request(url, data=urlencode(form).encode() if form else None, headers=headers)
    host = urlsplit(url).hostname
    for attempt in range(3):
        time.sleep(0.4)  # Rate limit: applies to every attempt, including retries.
        try:
            with urlopen(request, timeout=30) as response:
                return response.read(), response.headers
        except HTTPError as error:
            problem = ServiceError(f"{host} returned HTTP {error.code}.", error.code)
            delay = retry_delay(error.headers.get("Retry-After"), attempt)
            error.close()
            # Fail at once on non-transient errors, or if the server asks us
            # to wait longer than a minute.
            if error.code not in (408, 429, 500, 502, 503, 504) or delay > 60:
                raise problem from None
        except (URLError, OSError, HTTPException) as error:
            problem = ServiceError(f"Could not reach {host} ({type(error).__name__}).")
            delay = 2 ** attempt
        if attempt == 2 or form:
            raise problem from None
        time.sleep(delay)


def decode(body):
    """Parse a JSON body and raise ServiceError if it is a top-level error object."""
    data = json.loads(body)
    if isinstance(data, dict) and (data.get("error") or data.get("errors") or data.get("ERROR")
                                  or str(data.get("code", 200)) != "200"):
        raise ServiceError("The service returned an error response.")
    return data


def get_json(url, params=None, form=None):
    """Fetch a URL and return its decoded JSON body."""
    return decode(fetch(url, params, form)[0])


# NCBI lookups
@cache  # One efetch per Gene ID per run, even if the gene is checked twice.
def gene_ensembl_ids(gene_id, prefix):
    """Return the sorted, unversioned Ensembl gene IDs linked to an NCBI Gene ID.

    ESummary does not include Ensembl cross-references, so the full Entrezgene
    XML record is fetched and its Gene-ref database tags are read. `prefix`
    restricts the result to one species (e.g. "ENSG" for human). The record is
    checked to be the gene that was requested before anything is returned.
    """
    body, _ = fetch(EUTILS + "efetch.fcgi", {"db": "gene", "id": gene_id, "retmode": "xml"})
    root = ET.fromstring(body)
    # Strip any XML namespace so the plain NCBI tag names below can be used.
    for node in root.iter():
        node.tag = node.tag.rsplit("}", 1)[-1]
    record = root if root.tag == "Entrezgene" else root.find("Entrezgene")
    if record is None or record.findtext("Entrezgene_track-info/Gene-track/Gene-track_geneid") != gene_id:
        raise ServiceError("NCBI XML did not identify the requested gene.")
    gene = record.find("Entrezgene_gene/Gene-ref")
    if gene is None:
        raise ServiceError("NCBI XML is missing the gene record.")
    ids = set()
    for tag in gene.findall("Gene-ref_db/Dbtag"):
        value = tag.findtext("Dbtag_tag/Object-id/Object-id_str", "").upper()
        # Keep only Ensembl tags whose value is a well-formed ID for this
        # species; drop any version suffix (ENSG...12.3 -> ENSG...12).
        if tag.findtext("Dbtag_db", "").casefold().startswith("ensembl") and re.fullmatch(prefix + r"\d+(\.\d+)?", value):
            ids.add(value.split(".")[0])
    return sorted(ids)


# File helpers

def part_dir(data_root, part, dataset):
    """Return the output folder for one part and dataset (Part 1 uses final_data/)."""
    folder = data_root / f"part{part}_data"
    return (folder / "final_data" if part == 1 else folder) / dataset


def read_table(path, plain=False):
    """Read a TSV into a list of dicts, failing on an empty or ragged table.

    Not used in Part 1; Part 2 uses it to read the Part 1 TSV. Use plain=True
    for files written with write_tsv(plain=True), so quotation marks in the
    text are read literally rather than as CSV quoting.
    """
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", quoting=csv.QUOTE_NONE if plain else csv.QUOTE_MINIMAL)
        rows = list(reader)
    # A None key or value means a row had more or fewer fields than the header.
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError(f"Empty or malformed input table: {path}")
    return rows


def write_tsv(path, columns, rows, plain=False):
    """Write rows to a TSV atomically.

    The file is written to a temporary folder beside the target and then moved
    into place, so an interrupted run never leaves a half-written table.
    With plain=True (Part 1), values are written unquoted and any tab or line
    break inside a value is replaced by a space, so every row stays on one
    line. Otherwise (Part 2) standard CSV quoting is used. Missing values are
    written as empty cells.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        with staged.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for row in rows:
                values = [row.get(column) for column in columns]
                if plain:
                    cells = [re.sub(r"[\r\n\t]", " ", str(v)) if v is not None else "" for v in values]
                    stream.write("\t".join(cells) + "\n")
                else:
                    writer.writerow(values)
        staged.replace(path)


# Command-line entry point
def run(process):
    """Parse arguments, run `process` for each selected dataset, and set the exit code.

    `process(data_root, dataset)` returns True when the dataset was saved but
    has lookup problems; the run then continues and exits with code 2.
    A fatal error (file, CSV or unexpected data error) stops the run with code 1.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", "--dataset", choices=("both", *DATASETS), default="both")
    parser.add_argument("--project-root", type=Path, default=Path(os.getenv("BDB_PROJECT_ROOT", "~/BDB_ICA")))
    args = parser.parse_args()
    failed = False
    try:
        for dataset in DATASETS if args.run == "both" else (args.run,):
            if process(args.project_root.expanduser() / "OP1_data", dataset):
                print(f"[WARN] {dataset} saved with lookup problems. See its status and note columns.", file=sys.stderr)
                failed = True
    except (KeyError, ValueError, OSError, csv.Error) as error:
        sys.exit(f"[ERROR] {error}")  # Message to stderr, exit code 1.
    sys.exit(2 if failed else 0)


# Gene resolution
def gene_summaries(ids, known):
    """Return NCBI ESummary records for `ids`, in the same order.

    `known` is a per-dataset cache (Gene ID -> record): only IDs not already in
    it are requested, in a single batched ESummary call. Each returned record
    must carry the requested UID and no error field.
    """
    missing = [gid for gid in ids if gid not in known]
    if missing:
        params = {"db": "gene", "id": ",".join(missing), "retmode": "json"}
        result = get_json(EUTILS + "esummary.fcgi", params)["result"]
        for gid in missing:
            record = result[gid]
            if str(record.get("uid")) != gid or record.get("error"):
                raise ServiceError(f"No usable NCBI summary for Gene ID {gid}.")
            known[gid] = record
    return [known[gid] for gid in ids]


def official(record):
    """Return the official (HGNC) symbol, falling back to the NCBI gene name."""
    return record.get("nomenclaturesymbol") or record.get("name") or ""


def current_human(record):
    """True if the record is a live human gene (taxon 9606, not replaced or withdrawn).

    In ESummary, status "0" means live and a non-zero currentid means the
    record has been merged into another Gene ID.
    """
    taxon = (record.get("organism") or {}).get("taxid") or record.get("taxid")
    if not taxon:
        raise ServiceError("NCBI did not supply a taxonomy ID.")
    return (str(taxon) == "9606" and str(record.get("status") or "0") == "0"
            and str(record.get("currentid") or "0") == "0")


def resolve(label, known):
    """Resolve one input label to a single current human NCBI Gene record.

    Returns (status, note, record), where status is "resolved", "ambiguous"
    or "not found" and record is the ESummary record (None unless resolved).

    Strategy:
      * Ensembl IDs (ENSG..., optionally versioned) are searched in all fields,
        and a candidate is accepted only if its own Ensembl cross-references
        contain that exact ID.
      * Other labels are searched as gene names, then matched in order of
        confidence: official symbol, then alias, then (if the search returned
        exactly one current human gene) that single hit.
    A label is never guessed when more than one candidate fits equally well.
    """
    # Reject labels that cannot be a gene symbol or ID before querying NCBI
    if not re.fullmatch(r"[\w.-]+", label):
        return "not found", "Expected a gene symbol, alias or human Ensembl gene ID.", None
    ensembl = re.fullmatch(r"(ENSG\d+)(?:\.\d+)?", label.upper())
    query = f'"{ensembl[1]}"[All Fields]' if ensembl else f'{label}[Gene Name]'
    params = {"db": "gene", "term": query + ' AND "Homo sapiens"[Organism] AND alive[prop]',
              "retmode": "json", "retmax": 100}
    search = get_json(EUTILS + "esearch.fcgi", params)["esearchresult"]
    ids, count = search["idlist"], int(search["count"])

    errors = search.get("errorlist") or {}
    if search.get("error") or any(v for k, v in errors.items() if k != "phrasesnotfound") or (count and errors.get("phrasesnotfound")):
        raise ServiceError("NCBI could not interpret the full search.")
    if not isinstance(ids, list) or any(not re.fullmatch(r"[1-9]\d*", gid) for gid in ids) or count < len(ids):
        raise ServiceError("Invalid NCBI search result.")
    # More hits than retmax: the true match may not be in the returned page.
    if count > len(ids):
        return "ambiguous", "Too many candidates to select one safely.", None
    candidates = [r for r in gene_summaries(ids, known) if current_human(r)]
    if ensembl:
        matches = [r for r in candidates if ensembl[1] in gene_ensembl_ids(str(r["uid"]), "ENSG")]
        reason = "Verified the exact Ensembl gene cross-reference."
    else:
        matches = [r for r in candidates if official(r).casefold() == label.casefold()]
        reason = "Matched the official symbol."
        if not matches:
            # otheraliases is a comma-separated string, e.g. "MRTL, MYCC, c-Myc".
            matches = [r for r in candidates if label.casefold() in
                       {alias.strip().casefold() for alias in (r.get("otheraliases") or "").split(",")}]
            reason = "Matched an alias."
        if not matches and len(candidates) == 1:
            matches, reason = candidates, "One current human gene matched the gene-name search."
    if len(matches) == 1:
        return "resolved", reason, matches[0]
    if len(matches) > 1 or (not ensembl and len(candidates) > 1):
        return "ambiguous", "Several candidates without a unique match.", None
    return "not found", "No matching current human gene.", None


def annotate(label, known):
    """Build the output row for one input label.

    Any API problem is recorded in the row as status "failed" with the error
    message, so the dataset still produces one row per input. If only the
    Ensembl lookup fails, the gene stays "resolved" and the failure is
    appended to Resolution_Note.
    """
    row = {"Query_Symbol": label}
    try:
        status, note, record = resolve(label, known)
        row.update(Resolution_Status=status, Resolution_Note=note)
        if record is None:
            return row
        if not official(record):
            raise ServiceError("NCBI returned a gene without a symbol.")
        gid = str(record["uid"])
        row.update(Gene_ID=gid, Official_Symbol=official(record),
                   Official_Name=record.get("nomenclaturename") or record.get("description"),
                   Map_Location=record.get("maplocation"), Other_Aliases=record.get("otheraliases"),
                   Summary=record.get("summary"))
        try:
            # Cached, so Ensembl-ID inputs do not trigger a second request here.
            row["ENSG_IDs"] = ";".join(gene_ensembl_ids(gid, "ENSG"))
        except API_ERRORS as error:
            row["Resolution_Note"] += f" Ensembl cross-reference lookup failed: {error}"
    except API_ERRORS as error:
        row.update(Resolution_Status="failed", Resolution_Note=str(error))
    return row


def process(data_root, dataset):
    """Annotate every input in one dataset, write its TSV, and report problems.

    Repeated genes keep their own row; later repeats get a Same_Gene_Note
    pointing to the first input row (1-based) with the same Gene ID.

    Returns True if the saved output needs attention: for the development
    list, any input that is not resolved or lacks an ENSG ID (every gene there
    is expected to resolve); for the challenge list, only lookups that failed,
    since unresolved and ambiguous inputs are expected outcomes there.
    """
    known, first_input, rows = {}, {}, []
    for label in GENES[dataset]:
        row = annotate(label, known)
        gid = row.get("Gene_ID")
        if gid in first_input:
            row["Same_Gene_Note"] = f"Same Gene ID as input {first_input[gid]}."
        elif gid:
            first_input[gid] = len(rows) + 1
        rows.append(row)
    write_tsv(part_dir(data_root, 1, dataset) / "ncbi_gene_summary2.tsv", COLUMNS, rows, plain=True)
    resolved = sum(row["Resolution_Status"] == "resolved" for row in rows)
    print(f"Part 1 {dataset}: {resolved}/{len(rows)} inputs resolved to {len(first_input)} distinct genes.")
    if dataset == "development" and (resolved != len(rows) or any(not row.get("ENSG_IDs") for row in rows)):
        return True
    return any(row["Resolution_Status"] == "failed" or "lookup failed:" in row["Resolution_Note"] for row in rows)


if __name__ == "__main__":
    run(process)
