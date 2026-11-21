#number1 
# #!/usr/bin/env python3
"""BDB ICA Part 2: mouse orthologs from the final Part 1 annotation results.

BASE AND ATTRIBUTION
Base: Simon Tomlinson's supplied OrthologMap3_SRT.py (November 2026).
The supplied header states: "This script was generated using ELM
(elm.edina.ed.ac.uk) and edited. Simon Tomlinson Nov 2026".
GenAI disclosure: OpenAI ChatGPT assisted with adapting and combining the
supplied ortholog scripts and the final Part 1 script, including validation,
request handling, output tracking, comments and checks.

WHY THIS BASE AND WHAT WAS COMBINED
- SRT base: keep its callable connect_ncbi_ortholog(), returned rows,
  configurable output prefix and NCBI ortholog/Ensembl workflow.
- OrthologMap3.py: restore an editable configuration and a guarded main()
  so the same file also runs directly from PyCharm, VS Code or a terminal.
- Final Part 1 script: read its exact JSON keys or TSV columns, reuse its
  resolved Gene IDs, preserve resolution notes, pace every request, and use
  finite timeouts, bounded retries and environment settings.
- Earlier edited ortholog scripts: add pagination, per-input outcomes,
  one-to-many rows and the current nested ESummary organism.taxid field.
- The lecturer's symbol-resolution step is replaced by reading Part 1.
  Searching symbols again could undo Part 1's careful alias/ENSG decisions.

RUNNING
Dependency: requests (the same external dependency used by Part 1).
    python -m pip install requests
    python p2_orthologs_tailored.py
The default run reads development and challenge JSON files at INPUT_ROOT,
then writes ncbi_ortholog.tsv and ncbi_ortholog.json in each OUTPUT_ROOT
subfolder. Part 1 files are read only. Edit the configuration below if needed.

For one selected file, including a Part 1 TSV:
    python p2_orthologs_tailored.py --input path/to/ncbi_gene_summary2.tsv \
        --output-prefix path/to/part2/ncbi_ortholog
Use --run development or --run challenge to run just that configured list.
For import use, pass a Part 1 filename OR its returned Python result records:
    rows = connect_ncbi_ortholog(part1_results, "output/ncbi_ortholog")
The first argument now contains Part 1 records, not a raw list of symbols.

OUTPUT
The first eight column names match the lecturer's ortholog table. Additional
columns explain the original input, its Part 1 outcome and its Part 2 outcome.
An unsuccessful input still has a row; multiple orthologs have separate rows.
No mouse ortholog returned means NCBI returned none in this lookup; it does
not establish that no biological ortholog exists.

SOURCES
- Supplied OrthologMap3_SRT.py and OrthologMap3.py.
- Supplied BD_ICA_ExtraTips_2per.pdf, slides 13-14.
- User's final_script_for_both_end_res (Part 1 input schema and file paths).
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/openapi3/openapi3.docs.yaml
"""

import argparse
import csv
import json
import logging
import math
import os
import re
import sys
import tempfile
import time
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import requests


# ---------------------------------------------------------------------------
# Configuration: paths match main() in the supplied final Part 1 script.
# These lists check input identity/order only; they never supply gene mappings.
# A custom --input run can contain any list and does not use these list checks.
# ---------------------------------------------------------------------------

INPUT_ROOT = Path("~/BDB_ICA/OP1_data/part1_data/final_data").expanduser()
OUTPUT_ROOT = Path("~/BDB_ICA/OP2_data/part2_data").expanduser()
PART1_FILENAME = "ncbi_gene_summary2.json"  # Change to .tsv to use that export.
OUTPUT_NAME = "ncbi_ortholog"             # Prefix, without .tsv or .json.
EXPECTED_INPUTS = {
    "development": [
        "BRCA1", "BRCA2", "TP53", "EGFR", "MYC",
        "APOE", "CFTR", "ACE2", "TNF", "VEGFA",
    ],
    "challenge": [
        "MYC", "c-Myc", "Blah", "ENSG00000136997", "AADACL2",
        "Plet1", "HYDIN", "ENSG00000261846", "XIST", "MALAT1",
    ],
}

HUMAN_TAX_ID = "9606"
MOUSE_TAX_ID = "10090"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_DATASETS_ORTHOLOGS_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/orthologs"
)
SLEEPBETWEEN = 2.0
MAX_ATTEMPTS = 3
TIMEOUT = 30
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}
LOG_FILE = "ncbi_ortholog.log"
LOGGER = logging.getLogger("ncbi_ortholog")
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False

# Retain the lecturer's names so subsequent readers can still find the IDs.
# Resolution_* comes from Part 1; Mapping_* describes this ortholog lookup.
OUTPUT_COLUMNS = [
    "Human_Input_ID", "Human_Gene_ID", "Human_Symbol", "Human_Ensembl_Gene_ID",
    "Mouse_Gene_ID", "Mouse_Symbol", "Mouse_Tax_ID", "Mouse_Ensembl_Gene_ID",
    "Mouse_Official_Name", "Input_Row", "Same_Gene_Note",
    "Resolution_Status", "Resolution_Note", "Mouse_Ortholog_Count",
    "Mapping_Status", "Mapping_Note", "Human_Ensembl_Status", "Mouse_Ensembl_Status",
]
SUCCESS_STATUSES = {
    "ortholog_found", "multiple_mouse_orthologs", "no_mouse_ortholog_returned",
}


class OrthologError(ValueError):
    """A validation or lookup error whose message is safe to display."""


def report(message: str, level: int = logging.WARNING) -> None:
    """Log fixed explanations, never raw request errors or credential URLs."""
    LOGGER.log(level, message)
    print(message, file=sys.stderr)


def validate_gene_id(value: Any) -> str:
    """Gene IDs must be positive integers, rather than symbols or queries."""
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise OrthologError("An NCBI Gene ID must be a positive whole number.")
    value = str(value).strip()
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise OrthologError("An NCBI Gene ID must be a positive whole number.")
    return value


def normalise_ensembl_ids(values: Any, prefix: str) -> List[str]:
    """Keep unique gene-level IDs, stripping optional version suffixes.

    Reject a transcript/protein ID or an ID for the wrong species instead of
    treating it as a usable gene cross-reference. Missing IDs remain [].
    """
    if not isinstance(values, list) or any(not isinstance(v, str) for v in values):
        raise OrthologError("An Ensembl gene ID field must be a list of text IDs.")
    result = []
    for value in values:
        match = re.fullmatch(rf"({prefix}[0-9]+)(?:\.[0-9]+)?", value.strip(), re.I)
        if not match:
            raise OrthologError("An Ensembl gene cross-reference has an invalid ID or species.")
        identifier = match.group(1).upper()
        if identifier not in result:
            result.append(identifier)
    return result


# ---------------------------------------------------------------------------
# Generic HTTP helper: SRT request wrapper strengthened using Part 1's policy.
# Put pacing here, before every request, rather than after processing rows.
# Both numeric and HTTP-date Retry-After headers are understood.
# ---------------------------------------------------------------------------

def retry_after_seconds(response: requests.Response) -> float:
    value = response.headers.get("Retry-After", "0")
    try:
        seconds = float(value)
    except (TypeError, ValueError):
        try:
            when = parsedate_to_datetime(value)
            if when.tzinfo is None:
                when = when.replace(tzinfo=timezone.utc)
            seconds = (when - datetime.now(timezone.utc)).total_seconds()
        except (TypeError, ValueError, OverflowError):
            return 0.0
    return max(0.0, seconds) if math.isfinite(seconds) else 0.0


def http_get(
    url: str, params: Optional[Dict[str, Any]] = None, timeout: float = TIMEOUT,
) -> requests.Response:
    """Return a successful response or raise; a failure is never an empty result.

    Retry read-only GETs for temporary failures only. Certificate failures,
    permanent HTTP errors and invalid data are not retried. Callers close the
    returned response after reading it, using a with statement.
    """
    datasets_request = bool(re.fullmatch(
        r"https://api\.ncbi\.nlm\.nih\.gov/datasets/v2/gene/id/[1-9][0-9]*/orthologs", url
    ))
    if not datasets_request and url not in {
        f"{NCBI_EUTILS_BASE}/esummary.fcgi", f"{NCBI_EUTILS_BASE}/efetch.fcgi",
    }:
        raise OrthologError("Unsupported NCBI request endpoint.")
    if params is not None and not isinstance(params, dict):
        raise OrthologError("NCBI request parameters must be a dictionary.")
    if isinstance(MAX_ATTEMPTS, bool) or not isinstance(MAX_ATTEMPTS, int) or not 1 <= MAX_ATTEMPTS <= 5:
        raise OrthologError("Use between one and five request attempts.")
    for value, name in [(SLEEPBETWEEN, "request delay"), (timeout, "timeout")]:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise OrthologError(f"The {name} must be a finite number.")
    if SLEEPBETWEEN < 0 or timeout <= 0:
        raise OrthologError("Use a nonnegative request delay and a positive timeout.")

    params = dict(params or {})
    headers = {"User-Agent": "bdb_ica_orthologs"}
    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if datasets_request:
        headers["Accept"] = "application/json"
        if api_key:
            headers["api-key"] = api_key
    else:
        params.setdefault("tool", "bdb_ica_orthologs")
        email = os.environ.get("NCBI_EMAIL", "").strip()
        if email:
            if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email):
                raise OrthologError("NCBI_EMAIL must contain a valid email address.")
            params.setdefault("email", email)
        if api_key:
            params["api_key"] = api_key

    for attempt in range(1, MAX_ATTEMPTS + 1):
        time.sleep(SLEEPBETWEEN)
        wait = 0.0
        try:
            response = requests.get(
                url, params=params, headers=headers, timeout=timeout, allow_redirects=False,
            )
            if response.status_code == 200:
                return response
            status = response.status_code
            wait = retry_after_seconds(response)
            response.close()
            if status not in RETRYABLE_STATUS_CODES:
                raise OrthologError(f"NCBI returned HTTP {status}; this response was not retried.")
            report(f"[WARN] Temporary NCBI HTTP {status}, attempt {attempt}/{MAX_ATTEMPTS}.")
        except requests.exceptions.SSLError:
            raise OrthologError("HTTPS certificate verification failed; request not retried.") from None
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            report(f"[WARN] NCBI connection/timeout failure, attempt {attempt}/{MAX_ATTEMPTS}.")
        except requests.exceptions.RequestException:
            raise OrthologError("The NCBI request failed; request not retried.") from None
        if attempt < MAX_ATTEMPTS:
            wait = max(wait, float(2 ** (attempt - 1)))
            if wait > 60:
                raise OrthologError("NCBI requested a wait over 60 seconds; try again later.")
            time.sleep(wait)
    raise OrthologError("NCBI lookup failed after the configured request attempts.")


def response_json(response: requests.Response) -> Dict[str, Any]:
    """Validate JSON separately from networking so malformed data is not retried."""
    try:
        data = response.json()
    except ValueError:
        raise OrthologError("NCBI returned invalid JSON.") from None
    if not isinstance(data, dict) or any(data.get(k) for k in ("error", "errors", "error_message", "ERROR")):
        raise OrthologError("NCBI returned an API error or an unexpected JSON structure.")
    return data


# ---------------------------------------------------------------------------
# Step 0: read Part 1 instead of resolving its symbols/aliases/ENSG IDs again.
# The original source row is retained even when Part 1 could not resolve it.
# ---------------------------------------------------------------------------

def read_part1_results(filename: Union[str, Path]) -> List[Dict[str, Any]]:
    """Read the exact JSON schema or TSV columns written by final Part 1."""
    path = Path(filename).expanduser()
    if not path.is_file():
        raise OrthologError("Part 1 input file was not found. Check INPUT_ROOT/PART1_FILENAME or select a file with --input.")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        if path.suffix.lower() == ".json":
            records = json.load(handle)
        elif path.suffix.lower() == ".tsv":
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"Query_Symbol", "Gene_ID", "Official_Symbol", "ENSG_IDs", "Resolution_Status", "Resolution_Note", "Same_Gene_Note"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise OrthologError("Part 1 TSV is missing required final-script columns.")
            if len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise OrthologError("Part 1 TSV has duplicate column names.")
            records = []
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise OrthologError("Part 1 TSV has a row with an incorrect number of fields.")
                records.append({
                    "query_symbol": row["Query_Symbol"], "gene_id": row["Gene_ID"],
                    "summary": {"nomenclaturesymbol": row["Official_Symbol"]},
                    "ensembl_gene_ids": [v.strip() for v in row["ENSG_IDs"].split(";") if v.strip()],
                    "resolution_status": row["Resolution_Status"],
                    "resolution_note": row["Resolution_Note"],
                    "same_gene_note": row["Same_Gene_Note"],
                })
        else:
            raise OrthologError("Select the .json or .tsv file written by Part 1.")
    if not isinstance(records, list) or not records:
        raise OrthologError("Part 1 input must contain a nonempty list of annotation records.")
    return records


def prepare_part1_records(
    records: List[Dict[str, Any]], expected_inputs: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Check the input before requests; preserve order, aliases and repetitions."""
    if not isinstance(records, list) or not records:
        raise OrthologError("Supply a nonempty list of Part 1 annotation records.")
    prepared = []
    for position, record in enumerate(records, 1):
        if not isinstance(record, dict) or not {
            "query_symbol", "gene_id", "summary", "ensembl_gene_ids", "resolution_status",
        }.issubset(record):
            raise OrthologError("A Part 1 record is missing its required final-script fields.")
        if not isinstance(record["summary"], dict):
            raise OrthologError("A Part 1 summary must be a dictionary.")
        for key in ["query_symbol", "resolution_status", "resolution_note", "same_gene_note"]:
            if not isinstance(record.get(key, ""), str):
                raise OrthologError("Part 1 query labels, statuses and notes must contain text.")
        if not record["query_symbol"].strip() or not record["resolution_status"].strip():
            raise OrthologError("A Part 1 query label or resolution status is empty.")
        gene_id = record["gene_id"]
        if gene_id is None:
            gene_id = ""
        if isinstance(gene_id, bool) or not isinstance(gene_id, (str, int)):
            raise OrthologError("A Part 1 Gene ID must be text, an integer or empty.")
        summary = record["summary"]
        symbol = summary.get("nomenclaturesymbol") or summary.get("name") or ""
        if not isinstance(symbol, str):
            raise OrthologError("A Part 1 official symbol must contain text.")
        prepared.append({
            "input_id": record["query_symbol"], "gene_id": str(gene_id).strip(),
            "symbol": symbol, "input_row": position,
            "ensembl_gene_ids": normalise_ensembl_ids(record["ensembl_gene_ids"], "ENSG"),
            "resolution_status": record["resolution_status"],
            "resolution_note": record.get("resolution_note", ""),
            "same_gene_note": record.get("same_gene_note", ""),
        })
    if expected_inputs is not None and [r["input_id"] for r in prepared] != expected_inputs:
        raise OrthologError("Part 1 results do not contain the expected input rows in order.")
    return prepared


def fetch_gene_summary_by_id(gene_id: str) -> Dict[str, str]:
    """Retain SRT's species check, reading the current nested organism.taxid.

    We confirm the existing ID rather than select another ID. A replaced or
    discontinued record is reported and never silently redirected to a new one.
    """
    gene_id = validate_gene_id(gene_id)
    with http_get(f"{NCBI_EUTILS_BASE}/esummary.fcgi", {
        "db": "gene", "id": gene_id, "retmode": "json",
    }) as response:
        data = response_json(response)
    result = data.get("result")
    record = result.get(gene_id) if isinstance(result, dict) else None
    if not isinstance(record, dict) or not record or record.get("error"):
        raise OrthologError("NCBI did not return a usable summary for the supplied Gene ID.")
    if str(record.get("uid", "")) != gene_id:
        raise OrthologError("The returned NCBI summary has a different Gene ID.")
    organism = record.get("organism") or {}
    if not isinstance(organism, dict):
        raise OrthologError("NCBI returned an invalid organism record.")
    taxid = str(organism.get("taxid") or record.get("taxid") or "")
    if not taxid.isdigit():
        raise OrthologError("The NCBI summary is missing its taxonomy ID.")
    symbol = record.get("nomenclaturesymbol") or record.get("name") or ""
    if not isinstance(symbol, str):
        raise OrthologError("NCBI returned an invalid gene symbol.")
    return {
        "gene_id": gene_id, "symbol": symbol, "taxid": taxid,
        "status": str(record.get("status") or ""),
        "currentid": str(record.get("currentid") or ""),
    }


# ---------------------------------------------------------------------------
# Step 1: get mouse orthologs via the same NCBI Datasets endpoint as SRT.
# Ask for mouse on the server, check species locally, and follow every page.
# A successful empty response is different from a failed/malformed response.
# ---------------------------------------------------------------------------

def fetch_orthologs_for_gene(gene_id: str) -> Dict[str, Any]:
    """Return a validated, combined reports list for this human Gene ID."""
    gene_id = validate_gene_id(gene_id)
    url = NCBI_DATASETS_ORTHOLOGS_URL.format(gene_id=gene_id)
    params: Dict[str, Any] = {"taxon_filter": MOUSE_TAX_ID, "page_size": 100}
    reports_by_id: Dict[str, Dict[str, Any]] = {}
    expected_count = None
    seen_tokens = set()
    for _ in range(100):
        with http_get(url, params=params) as response:
            data = response_json(response)
        if "total_count" in data:
            value = data["total_count"]
            if isinstance(value, bool) or not isinstance(value, (str, int)) or not re.fullmatch(r"[0-9]+", str(value)):
                raise OrthologError("NCBI returned an invalid ortholog result count.")
            count = int(value)
            if expected_count is not None and count != expected_count:
                raise OrthologError("NCBI returned an inconsistent ortholog result count.")
            expected_count = count
        reports = data.get("reports")
        # Protobuf JSON can omit default-valued fields; {} is a valid zero
        # result. A nonempty unexplained response is not accepted as zero.
        if reports is None and (not data or expected_count == 0):
            reports = []
        if not isinstance(reports, list):
            raise OrthologError("NCBI response is missing its ortholog report list.")
        for item in reports:
            gene = item.get("gene") if isinstance(item, dict) else None
            if not isinstance(gene, dict) or item.get("errors"):
                raise OrthologError("NCBI returned an invalid or unsuccessful ortholog report.")
            identifier = validate_gene_id(gene.get("gene_id"))
            if str(gene.get("tax_id", "")) != MOUSE_TAX_ID:
                raise OrthologError("NCBI returned a species outside the requested mouse filter.")
            if identifier in reports_by_id and reports_by_id[identifier] != item:
                raise OrthologError("NCBI returned conflicting reports for the same mouse gene.")
            reports_by_id[identifier] = item
        token = data.get("next_page_token", "")
        if not isinstance(token, str):
            raise OrthologError("NCBI returned an invalid pagination token.")
        if not token:
            break
        if token in seen_tokens:
            raise OrthologError("NCBI repeated a pagination token.")
        seen_tokens.add(token)
        params["page_token"] = token
    else:
        raise OrthologError("NCBI ortholog pagination exceeded the configured page limit.")
    if expected_count is not None and len(reports_by_id) != expected_count:
        raise OrthologError("NCBI did not return the complete ortholog result set.")
    return {"reports": list(reports_by_id.values())}


def extract_mouse_orthologs_from_json(
    orthologs_json: Dict[str, Any], mouse_tax_id: str = MOUSE_TAX_ID,
) -> List[Dict[str, Any]]:
    """Retain all reported mouse genes, including their names and Ensembl IDs."""
    mice = []
    for item in orthologs_json["reports"]:
        gene = item["gene"]
        if str(gene["tax_id"]) != str(mouse_tax_id):
            raise OrthologError("An ortholog report does not belong to mouse.")
        symbol, name = gene.get("symbol", ""), gene.get("description", "")
        if not isinstance(symbol, str) or not isinstance(name, str):
            raise OrthologError("NCBI returned an invalid mouse symbol or name.")
        mice.append({
            "mouse_gene_id": validate_gene_id(gene["gene_id"]),
            "mouse_symbol": symbol, "mouse_name": name,
            "mouse_tax_id": MOUSE_TAX_ID,
            "ensembl_gene_ids": normalise_ensembl_ids(gene.get("ensembl_gene_ids", []), "ENSMUSG"),
        })
    return mice


def build_human_mouse_ortholog_rows(
    resolved_genes: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Preserve SRT's row builder, adding explicit outcomes instead of skipping.

    Its input includes unresolved Part 1 records as well as resolved genes.
    Cache each unique Gene ID within this call: MYC, c-Myc and a matching ENSG
    share requests, while each original input still receives its own rows.
    """
    rows = []
    summaries: Dict[str, Any] = {}
    orthologs: Dict[str, Any] = {}
    for record in resolved_genes:
        row = {key: "" for key in OUTPUT_COLUMNS}
        ensembl = ";".join(record["ensembl_gene_ids"])
        row.update({
            "Human_Input_ID": record["input_id"], "Human_Gene_ID": record["gene_id"],
            "Human_Symbol": record["symbol"], "Human_Ensembl_Gene_ID": ensembl,
            "Input_Row": record["input_row"], "Same_Gene_Note": record["same_gene_note"],
            "Resolution_Status": record["resolution_status"], "Resolution_Note": record["resolution_note"],
            "Human_Ensembl_Status": "from_part1" if ensembl else "not_requested",
            "Mouse_Ensembl_Status": "not_requested",
        })
        print(f"[INFO] Part 2 input row {record['input_row']}: {record['input_id']}")
        if record["resolution_status"] != "resolved":
            row.update(Mapping_Status="input_not_resolved", Mapping_Note="Not queried: Part 1 did not mark this input as resolved; see Resolution_Status and Resolution_Note.")
            rows.append(row)
            continue
        try:
            gene_id = validate_gene_id(record["gene_id"])
        except OrthologError:
            row.update(Mapping_Status="invalid_human_gene_id", Mapping_Note="Part 1 marked this input resolved but did not provide a positive NCBI Gene ID.")
            rows.append(row)
            continue
        # Store either a successful response or its safe exception, so an alias
        # does not immediately repeat a failed request three more times.
        if gene_id not in summaries:
            try:
                summaries[gene_id] = fetch_gene_summary_by_id(gene_id)
            except OrthologError as error:
                summaries[gene_id] = error
        summary = summaries[gene_id]
        if isinstance(summary, OrthologError):
            row.update(Mapping_Status="human_lookup_failed", Mapping_Note=str(summary))
        elif summary["taxid"] != HUMAN_TAX_ID:
            row.update(Mapping_Status="input_gene_is_not_human", Mapping_Note="The supplied Gene ID does not belong to Homo sapiens; no ortholog query was made.")
        elif summary["status"] not in {"", "0"} or summary["currentid"] not in {"", "0"}:
            row.update(Mapping_Status="human_record_not_current", Mapping_Note="The supplied human record is discontinued or replaced; no alternative ID was selected.")
        else:
            if not row["Human_Symbol"]:
                row["Human_Symbol"] = summary["symbol"]
            if gene_id not in orthologs:
                try:
                    orthologs[gene_id] = extract_mouse_orthologs_from_json(fetch_orthologs_for_gene(gene_id))
                except OrthologError as error:
                    orthologs[gene_id] = error
            mice = orthologs[gene_id]
            if isinstance(mice, OrthologError):
                row.update(Mapping_Status="ortholog_lookup_failed", Mapping_Note=str(mice))
            elif not mice:
                row.update(Mouse_Ortholog_Count=0, Mapping_Status="no_mouse_ortholog_returned", Mapping_Note="NCBI successfully returned no mouse ortholog for this human Gene ID.")
            else:
                for mouse in mice:
                    item = dict(row)
                    mouse_ensembl = ";".join(mouse["ensembl_gene_ids"])
                    item.update({
                        "Mouse_Gene_ID": mouse["mouse_gene_id"], "Mouse_Symbol": mouse["mouse_symbol"],
                        "Mouse_Tax_ID": MOUSE_TAX_ID, "Mouse_Official_Name": mouse["mouse_name"],
                        "Mouse_Ensembl_Gene_ID": mouse_ensembl, "Mouse_Ortholog_Count": len(mice),
                        "Mouse_Ensembl_Status": "from_ncbi_datasets" if mouse_ensembl else "not_requested",
                        "Mapping_Status": "ortholog_found" if len(mice) == 1 else "multiple_mouse_orthologs",
                        "Mapping_Note": "Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.",
                    })
                    rows.append(item)
                continue
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Step 2: retain the lecturer's EFetch/XML route as a fallback for missing IDs.
# Reuse Part 1 human ENSG IDs and mouse IDs already provided by NCBI Datasets.
# Only remaining gaps need an additional, batched NCBI request.
# ---------------------------------------------------------------------------

def strip_namespaces(root: ET.Element) -> None:
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]


def efetch_gene_xml(gene_ids: List[str]) -> Dict[str, str]:
    """Fetch distinct Gene IDs in batches and verify the returned record IDs."""
    unique = list(dict.fromkeys(validate_gene_id(gid) for gid in gene_ids))
    xml_by_id = {}
    for start in range(0, len(unique), 100):
        batch = unique[start:start + 100]
        with http_get(f"{NCBI_EUTILS_BASE}/efetch.fcgi", {
            "db": "gene", "id": ",".join(batch), "retmode": "xml",
        }) as response:
            text = response.text
        try:
            root = ET.fromstring(text)
        except ET.ParseError:
            raise OrthologError("NCBI returned invalid Gene XML.") from None
        strip_namespaces(root)
        entries = [root] if root.tag == "Entrezgene" else root.findall(".//Entrezgene")
        batch_records = {}
        for entry in entries:
            gid = (entry.findtext(".//Gene-track_geneid") or "").strip()
            if gid not in batch or gid in batch_records:
                raise OrthologError("NCBI XML contained an unexpected or repeated Gene ID.")
            batch_records[gid] = ET.tostring(entry, encoding="unicode")
        if set(batch_records) != set(batch):
            raise OrthologError("NCBI did not return every requested Gene XML record.")
        xml_by_id.update(batch_records)
    return xml_by_id


def extract_ensembl_ids_from_gene_xml(xml_text: str) -> List[str]:
    """Read gene-level Ensembl cross-references, excluding linked-gene noise.

    Searching every Dbtag in the full record can pick up IDs from other linked
    genes. The own-gene Gene-ref_db path follows Part 1's stronger ENSG check.
    Namespace removal supports either namespaced or plain Gene XML.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        raise OrthologError("NCBI returned invalid Gene XML.") from None
    strip_namespaces(root)
    if root.find(".//Entrezgene_gene/Gene-ref") is None:
        raise OrthologError("NCBI XML did not contain the expected own-gene record.")
    ids = []
    for dbtag in root.findall(".//Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag"):
        database = (dbtag.findtext("Dbtag_db") or "").strip()
        if not database.casefold().startswith("ensembl"):
            continue
        identifier = dbtag.findtext("Dbtag_tag/Object-id/Object-id_str", default="").strip()
        match = re.fullmatch(r"(ENSG[0-9]+|ENSMUSG[0-9]+)(?:\.[0-9]+)?", identifier, re.I)
        if match and match.group(1).upper() not in ids:
            ids.append(match.group(1).upper())
    return ids


def add_ensembl_ids_to_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Fill missing identifiers; keep relationship and metadata failures separate."""
    updated = [dict(row) for row in rows]
    needed: Dict[str, str] = {}
    for row in updated:
        if row["Mapping_Status"] not in SUCCESS_STATUSES:
            continue
        for species, prefix in [("Human", "ENSG"), ("Mouse", "ENSMUSG")]:
            gid = row[f"{species}_Gene_ID"]
            if gid and not row[f"{species}_Ensembl_Gene_ID"]:
                needed[gid] = prefix
    if not needed:
        return updated
    # One lookup per distinct ID, shared by all rows for the same gene.
    lookups: Dict[str, Any] = {}
    try:
        xml_by_id = efetch_gene_xml(list(needed))
    except OrthologError as error:
        report(f"[WARN] Ensembl fallback: {error}")
        lookups = {gid: error for gid in needed}
    else:
        for gid, prefix in needed.items():
            try:
                ids = extract_ensembl_ids_from_gene_xml(xml_by_id[gid])
                lookups[gid] = [identifier for identifier in ids if identifier.startswith(prefix)]
            except OrthologError as error:
                lookups[gid] = error
    for row in updated:
        if row["Mapping_Status"] not in SUCCESS_STATUSES:
            continue
        for species in ["Human", "Mouse"]:
            gid = row[f"{species}_Gene_ID"]
            if gid not in lookups or row[f"{species}_Ensembl_Gene_ID"]:
                continue
            value = lookups[gid]
            if isinstance(value, OrthologError):
                row[f"{species}_Ensembl_Status"] = "lookup_failed"
                row["Mapping_Note"] += f" {species} Ensembl lookup failed: {value}"
            else:
                row[f"{species}_Ensembl_Gene_ID"] = ";".join(value)
                row[f"{species}_Ensembl_Status"] = "from_ncbi_xml" if value else "no_id_returned"
    return updated


# ---------------------------------------------------------------------------
# Output and reusable entry point: SRT wrapper plus standalone-script strengths.
# Prepare both exports before replacing either; each replacement is atomic.
# The pair of replacements is not a database transaction.
# ---------------------------------------------------------------------------

def write_results(rows: List[Dict[str, Any]], prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ortholog_", dir=prefix.parent) as folder:
        temporary = Path(folder)
        with (temporary / "results.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with (temporary / "results.json").open("w", encoding="utf-8") as handle:
            json.dump(rows, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        for extension in [".tsv", ".json"]:
            destination = Path(str(prefix) + extension)
            os.replace(temporary / ("results" + extension), destination)
            print(f"[INFO] Saved {destination.resolve()}")


def connect_ncbi_ortholog(
    part1_results: Union[str, Path, List[Dict[str, Any]]],
    ncbi_ortholog_file: Union[str, Path],
    expected_inputs: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Read/accept final Part 1 results, map them, save TSV/JSON and return rows.

    Unlike the original wrapper, the first argument is Part 1 records or a
    Part 1 filename. Raw symbol lists are rejected: ID resolution belongs to
    Part 1. All-unresolved input is valid and still produces explanatory rows.
    """
    prefix = Path(ncbi_ortholog_file).expanduser()
    if isinstance(part1_results, (str, Path)):
        source = Path(part1_results).expanduser().resolve()
        protected = {source.with_suffix(".json"), source.with_suffix(".tsv")}
        if any(Path(str(prefix) + ext).resolve() in protected for ext in [".json", ".tsv"]):
            raise OrthologError("Choose a Part 2 output prefix that cannot overwrite either Part 1 export.")
        records = read_part1_results(source)
    else:
        records = part1_results
    prepared = prepare_part1_records(records, expected_inputs)
    rows = add_ensembl_ids_to_rows(build_human_mouse_ortholog_rows(prepared))
    # Check source-row coverage, allowing more than one output per input.
    covered = list(dict.fromkeys(row["Input_Row"] for row in rows))
    if covered != list(range(1, len(prepared) + 1)):
        raise OrthologError("Ortholog results lost or reordered an input row; no output was written.")
    write_results(rows, prefix)
    pairs = {(row["Human_Gene_ID"], row["Mouse_Gene_ID"]) for row in rows if row["Mouse_Gene_ID"]}
    print(f"[INFO] Preserved {len(prepared)} input rows; returned {len(pairs)} distinct human-mouse pairs in {len(rows)} output rows.")
    return rows


def has_lookup_failures(rows: List[Dict[str, Any]]) -> bool:
    """A missing biological match is allowed; a technical failure is incomplete."""
    return any(
        row["Mapping_Status"] not in SUCCESS_STATUSES | {"input_not_resolved"}
        or row["Resolution_Status"] == "failed"
        or row["Human_Ensembl_Status"] == "lookup_failed"
        or row["Mouse_Ensembl_Status"] == "lookup_failed"
        for row in rows
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="BDB ICA Part 2: use final Part 1 records to retrieve NCBI mouse orthologs.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--run", choices=["both", "development", "challenge"], help="Configured lists; default: both")
    group.add_argument("--input", help="A custom Part 1 JSON or TSV file")
    parser.add_argument("--output-prefix", help="Custom Part 2 prefix without extension; required with --input")
    parser.add_argument("--input-root", default=str(INPUT_ROOT), help="Parent of the Part 1 development/challenge folders")
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT), help="Parent of the Part 2 development/challenge folders")
    args = parser.parse_args()
    if args.input:
        if not args.output_prefix:
            parser.error("--input requires --output-prefix")
        rows = connect_ncbi_ortholog(args.input, args.output_prefix)
        if has_lookup_failures(rows):
            raise OrthologError("Results were saved with lookup failures; inspect the status/note columns before retrying.")
        return
    if args.output_prefix:
        parser.error("--output-prefix is used with --input; otherwise select --output-root")

    runs = ["development", "challenge"] if args.run in {None, "both"} else [args.run]
    sources = {run: Path(args.input_root).expanduser() / run / PART1_FILENAME for run in runs}
    # Check both configured input files before making the first network call.
    for run in runs:
        prepare_part1_records(read_part1_results(sources[run]), EXPECTED_INPUTS[run])
    for run in runs:
        print(f"[INFO] Starting Part 2 {run} run.")
        prefix = Path(args.output_root).expanduser() / run / OUTPUT_NAME
        rows = connect_ncbi_ortholog(sources[run], prefix, EXPECTED_INPUTS[run])
        incomplete = has_lookup_failures(rows)
        if run == "development":
            incomplete = incomplete or any(row["Resolution_Status"] != "resolved" for row in rows)
        if incomplete:
            # Save explanations first. Do not let an incomplete development
            # run trigger an automatic challenge run, matching Part 1's gate.
            raise OrthologError(f"Part 2 {run} results were saved but contain incomplete lookups; further runs were not started.")
        print(f"[INFO] Part 2 {run} checks passed. Successful empty ortholog results are allowed.")


if __name__ == "__main__":
    try:
        handler = logging.FileHandler(LOG_FILE, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        LOGGER.setLevel(logging.INFO)
        LOGGER.addHandler(handler)
        main()
    except KeyboardInterrupt:
        report("[ERROR] Run cancelled.", logging.ERROR)
        sys.exit(130)
    except Exception as error:
        LOGGER.error("Run stopped: %s", type(error).__name__)
        for frame in traceback.extract_tb(error.__traceback__):
            LOGGER.error("%s:%s in %s", os.path.basename(frame.filename), frame.lineno, frame.name)
        message = str(error) if isinstance(error, OrthologError) else "Run stopped. Check the input, output paths and ncbi_ortholog.log."
        report(f"[ERROR] {message}", logging.ERROR)
        sys.exit(1)