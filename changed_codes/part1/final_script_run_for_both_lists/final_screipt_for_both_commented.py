#!/usr/bin/env python3
"""
NCBI Gene connector: basic human gene annotation for the BDB ICA.

Resolves each input (gene symbol, alias or Ensembl gene ID) to a current NCBI
Gene record, retrieves its ESummary annotation and the Ensembl gene IDs listed
in the Gene XML cross-references, and saves the results as TSV and JSON.

Provenance: based on the course example p1_basic_anno.py (S. Tomlinson),
combined with p1_basic_anno_ref.py.
GenAI declaration: Claude (Anthropic) was used to help debug, refine and
condense this connector, as permitted for API connectors in this ICA.

Running the script annotates the development list first, checks the
results, and only then annotates the challenging list (see main()).

Usage:
    export NCBI_EMAIL="you@example.ac.uk"     # requested by NCBI usage policy
    export NCBI_API_KEY="..."                 # optional, raises the rate limit
    python ConnecttoNCBI.py
"""
# ---------------------------------------------------------------------------
# Revision markup (remove these markers once the changes are accepted):
#   ####       code that is not needed; it can be deleted without changing
#              behaviour.
#   ########   code that can be written more concisely; the replacement
#              follows directly underneath.
#   FIX:       a bug fix made in this revision.
# ---------------------------------------------------------------------------

import sys
import time
import json
import os
import re
import math
import logging
import traceback
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from typing import List, Dict, Any, Optional
import xml.etree.ElementTree as ET
import requests


# NCBI E-utilities settings. NCBI asks for at most 3 requests per second
# without an API key; a 2-second pause before each request stays well within it.
SLEEPBETWEEN = 2.0
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_TOOL = "bdb_ica_gene_annotation"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
NCBI_API_KEY: Optional[str] = os.environ.get("NCBI_API_KEY", "").strip() or None
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}
LOG_FILE = "ncbi_annotation.log"

######## # Define the NCBI service, contact details and a two-second request delay.
######## # Set NCBI_EMAIL in your environment before running; NCBI_API_KEY is optional.
######## SLEEPBETWEEN = 2.0
######## NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
######## NCBI_TOOL = "bdb_ica_gene_annotation"
######## NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
######## NCBI_API_KEY: Optional[str] = os.environ.get("NCBI_API_KEY", "").strip() or None
######## RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}
######## LOG_FILE = "ncbi_annotation.log"
# Diagnostic logger: silent until a file handler is attached in __main__.
LOGGER = logging.getLogger("ncbi_annotation")
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False

# Cache of Gene XML records already downloaded in this run (Gene ID -> XML).
_GENE_XML_CACHE: Dict[str, str] = {}


class AnnotationError(ValueError):
    """Invalid input or an unsuccessful run; the message is safe to display."""


# ---------------------------------------------------------------------------
# Input validation and reporting
# ---------------------------------------------------------------------------

######## # New: validate inputs before making requests. Do not accept query expressions.
def validate_gene_symbol(symbol: str) -> str:
    """Return a trimmed gene symbol, rejecting anything that could alter the
    Entrez query (spaces, brackets, quotes or Boolean operators)."""
    if not isinstance(symbol, str):
        raise AnnotationError("Each gene symbol must be text.")
    symbol = symbol.strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", symbol):
        raise AnnotationError("Gene symbols must be nonempty and use letters, digits, dots, underscores or hyphens.")
    return symbol


def validate_gene_id(gene_id: str) -> str:
    """Return the NCBI Gene ID as a string if it is a positive integer."""
    gene_id = str(gene_id)
    if not re.fullmatch(r"[1-9][0-9]*", gene_id):
        raise AnnotationError("An NCBI Gene ID must be a positive whole number.")
    return gene_id


######## # New: send the same safe message to the console and diagnostic log.
######## # Never pass raw exception text, response bodies, URLs or credentials here.
def report(message: str, level: int = logging.WARNING) -> None:
    """Write a message to stderr and the diagnostic log.

    Only pass fixed, safe text: never raw exception messages, response bodies,
    URLs or credentials (a request URL can contain the API key)."""
    LOGGER.log(level, message)
    print(message, file=sys.stderr)


######## def retry_after_seconds(response: requests.Response) -> float:
########     value = response.headers.get("Retry-After")
########     if value is None:
########         return 0.0
########     try:
########         seconds = float(value)
########         return max(0.0, seconds) if math.isfinite(seconds) else 0.0
########     except ValueError:
########         try:
########             when = parsedate_to_datetime(value)
########             if when.tzinfo is None:
########                 when = when.replace(tzinfo=timezone.utc)
########             return max(0.0, (when - datetime.now(timezone.utc)).total_seconds())
########         except (TypeError, ValueError, OverflowError):
########             return 0.0
def retry_after_seconds(response: requests.Response) -> float:
    """Return the server's Retry-After delay in seconds (0 if absent or not
    numeric). NCBI sends seconds; an HTTP-date value falls back to normal backoff."""
    try:
        return max(0.0, float(response.headers.get("Retry-After", 0)))
    except ValueError:
        return 0.0


# ---------------------------------------------------------------------------
# NCBI E-utilities requests
# ---------------------------------------------------------------------------

######## # Source: SRT request helper, with type hints and three attempts from the ref.
######## # New fix: pause before EVERY request, including successful requests.
######## # Retry only temporary read failures, using bounded exponential backoff.
# Source: SRT request helper; type hints and three attempts from the reference script.
def ncbi_get(
    endpoint: str,
    params: Dict[str, Any],
    max_retries: int = 3,
    sleep_between: float = SLEEPBETWEEN,
) -> Optional[requests.Response]:
    """Send a GET request to an E-utilities endpoint and return the response,
    or None on failure.

    Pauses before every request to respect NCBI rate limits. Only temporary
    failures (timeouts, connection errors, HTTP 408/429/5xx) are retried, with
    exponential backoff. max_retries is the total number of attempts."""
    #### # All callers are internal and pass fixed, valid arguments, so these
    #### # checks cannot fail in practice.
    if endpoint not in {"esearch.fcgi", "esummary.fcgi", "efetch.fcgi"}:
        raise AnnotationError("Unsupported NCBI endpoint.")
    #### if not isinstance(params, dict):
    ####     raise AnnotationError("NCBI request parameters must be a dictionary.")
    if not isinstance(max_retries, int) or not 1 <= max_retries <= 5:
        raise AnnotationError("The number of request attempts must be between 1 and 5.")
    if not isinstance(sleep_between, (int, float)) or not math.isfinite(sleep_between) or sleep_between < 0:
        raise AnnotationError("The request delay must be finite and nonnegative.")
    url = f"{NCBI_BASE}/{endpoint}"
    params = dict(params)
    params.setdefault("tool", NCBI_TOOL)

    if NCBI_EMAIL:
        if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", NCBI_EMAIL):
            raise AnnotationError("NCBI_EMAIL must contain a valid email address.")
        params.setdefault("email", NCBI_EMAIL)
    if NCBI_API_KEY:
        params.setdefault("api_key", NCBI_API_KEY)

    for attempt in range(1, max_retries + 1):
        time.sleep(sleep_between)
        retry_after = 0.0

        try:
            resp = requests.get(url, params=params, timeout=30, allow_redirects=False)
            if resp.status_code == 200:
                return resp

            status = resp.status_code
            retry_after = retry_after_seconds(resp)
            resp.close()
            if status not in RETRYABLE_STATUS_CODES:
                report(f"[ERROR] {endpoint} returned HTTP {status}; this response will not be retried.", logging.ERROR)
                return None
            report(
                f"[WARN] Attempt {attempt}/{max_retries} for {endpoint} "
                f"returned temporary HTTP {status}."
            )
        except requests.exceptions.SSLError:
            report("[ERROR] HTTPS certificate verification failed; the request will not be retried.", logging.ERROR)
            return None
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as error:
            report(
                f"[WARN] Attempt {attempt}/{max_retries} for {endpoint} "
                f"failed: {type(error).__name__}."
            )
        except requests.exceptions.RequestException as error:
            report(f"[ERROR] {endpoint} failed: {type(error).__name__}; this error will not be retried.", logging.ERROR)
            return None

        if attempt < max_retries:
            ######## # Extra waits grow 1, 2, 4, ... seconds; the five-second request
            ######## # delay above still applies before each attempt by default.
            # Backoff grows 1, 2, 4, ... seconds (or the server's Retry-After if
            # longer), in addition to the standard pause before each attempt.
            backoff = max(float(2 ** (attempt - 1)), retry_after)
            if backoff > 60:
                report("[ERROR] NCBI requested a longer wait; stop and try again later.", logging.ERROR)
                return None
            time.sleep(backoff)

    report(
        f"[ERROR] Failed to retrieve {endpoint} after {max_retries} attempts.",
        logging.ERROR,
    )
    return None


######## # Source: SRT ESearch function, with type hints from p1_basic_anno_ref.py.
######## # Resolve a gene symbol and organism, checking exact official-symbol matches.
######## # New additions: handle invalid JSON with ValueError and warn about ambiguity.
######## (This block described esearch_gene and has been moved above that function.)


def normalize_ensembl_gene_id(identifier: str) -> Optional[str]:
    """Return an Ensembl human gene ID (ENSG...) in upper case without its
    version suffix, or None if the input is not an ENSG ID."""
    if not isinstance(identifier, str):
        return None

    match = re.fullmatch(
        r"(ENSG[0-9]+)(?:\.[0-9]+)?",
        identifier.strip(),
        re.IGNORECASE,
    )
    return match.group(1).upper() if match else None


# Source: SRT ESearch function, with type hints from p1_basic_anno_ref.py.
# Extended to accept Ensembl gene IDs, to choose between several hits by exact
# official symbol or alias, and to report why an input could not be resolved.
def esearch_gene(
    symbol: str,
    organism: str = "Homo sapiens",
    resolution_details: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Resolve a gene symbol, alias or ENSG ID to one current NCBI Gene ID.

    Returns the Gene ID, or None if the input is not found, ambiguous or the
    request failed. If resolution_details is a dict, it is filled with
    'status' (resolved / not found / ambiguous / failed), 'note' and
    'candidate_gene_ids'."""
    symbol = validate_gene_symbol(symbol)

    if (
        not isinstance(organism, str)
        or not organism.strip()
        or any(c in organism for c in '\"[]\r\n')
    ):
        raise AnnotationError("Supply a valid organism name.")

    candidates = []

    # Record the outcome in resolution_details and return the Gene ID (or None).
    def finish(value, status, note):
        if status == "ambiguous" and candidates:
            note += f" Candidate Gene IDs: {', '.join(candidates)}."

        if resolution_details is not None:
            resolution_details.update({
                "status": status,
                "note": note,
                "candidate_gene_ids": list(candidates),
            })

        return value

    # ENSG inputs are searched in all fields; symbols in the Gene Name field.
    ensembl_id = normalize_ensembl_gene_id(symbol)

    if ensembl_id:
        lookup = f'"{ensembl_id}"[All Fields]'
    else:
        lookup = f"{symbol}[Gene Name]"

    term = (
        f'{lookup} AND "{organism.strip()}"[Organism] '
        "AND alive[prop]"
    )

    response = ncbi_get("esearch.fcgi", {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "retmax": 100,
    })

    if response is None:
        return finish(None, "failed", "ESearch request failed.")

    try:
        data = response.json()
    except ValueError:
        return finish(None, "failed", "ESearch returned invalid JSON.")

    # Validate the response structure before reading any values from it.
    if (
        not isinstance(data, dict)
        or not isinstance(data.get("esearchresult"), dict)
    ):
        return finish(
            None, "failed",
            "ESearch returned an unexpected response structure.",
        )

    search = data["esearchresult"]
    errors = search.get("errorlist") or {}

    if not isinstance(errors, dict):
        return finish(None, "failed", "Invalid ESearch error structure.")

    # 'phrasesnotfound' can accompany a valid search with zero results, so it
    # is handled separately from genuine query errors.
    other_errors = any(
        value for key, value in errors.items()
        if key != "phrasesnotfound"
    )

    if (
        data.get("error") or data.get("ERROR")
        or search.get("error") or other_errors
    ):
        return finish(None, "failed", "ESearch reported a query error.")

    ids = search.get("idlist")

    if not isinstance(ids, list) or any(
        not isinstance(gid, str)
        or not re.fullmatch(r"[1-9][0-9]*", gid)
        for gid in ids
    ):
        return finish(None, "failed", "Invalid Gene ID list.")

    candidates = list(dict.fromkeys(ids))   # Remove duplicates, keep order.

    try:
        count = int(search.get("count", len(candidates)))
    except (TypeError, ValueError):
        return finish(None, "failed", "Invalid search match count.")

    if count < len(candidates) or count < 0:
        return finish(None, "failed", "Inconsistent search match count.")

    if count == 0:
        return finish(
            None, "not found",
            "No current NCBI Gene record matched this input "
            "in the selected organism.",
        )

    if errors.get("phrasesnotfound"):
        return finish(
            None, "failed",
            "ESearch only partially interpreted the query.",
        )

    if count > len(candidates):
        return finish(
            None, "ambiguous",
            "Not all candidates were retrieved; no Gene ID was selected.",
        )

    # ENSG input: accept only candidates whose gene-level Ensembl
    # cross-reference matches exactly. The XML is cached by efetch_gene_xml,
    # so it is not downloaded again when the Ensembl IDs are extracted later.
    if ensembl_id:
        matches = []

        for gene_id in candidates:
            xml_text = efetch_gene_xml(gene_id)

            if xml_text is None:
                return finish(
                    None, "failed",
                    "Could not retrieve every candidate's XML "
                    "to verify the Ensembl ID.",
                )

            try:
                root = ET.fromstring(xml_text)
            except ET.ParseError:
                return finish(None, "failed", "Candidate Gene XML was invalid.")

            # Strip any XML namespaces so the plain element paths below match.
            for element in root.iter():
                element.tag = element.tag.rsplit("}", 1)[-1]

            if root.find(".//Entrezgene_gene/Gene-ref") is None:
                return finish(
                    None, "failed",
                    "Candidate XML did not contain the expected Gene record.",
                )

            path = ".//Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag"

            for dbtag in root.findall(path):
                database = (dbtag.findtext("Dbtag_db") or "").strip()

                if database.casefold() != "ensembl":
                    continue

                linked_id = dbtag.findtext(
                    "Dbtag_tag/Object-id/Object-id_str",
                    default="",
                )

                if normalize_ensembl_gene_id(linked_id) == ensembl_id:
                    matches.append(gene_id)
                    break

        if len(matches) == 1:
            return finish(
                matches[0], "resolved",
                "Verified the exact Ensembl gene cross-reference.",
            )

        if not matches:
            return finish(
                None, "not found",
                "No candidate contained the exact Ensembl gene cross-reference.",
            )

        return finish(
            None, "ambiguous",
            "Several current records contained the exact Ensembl cross-reference.",
        )

    # Symbol input with a single hit: accept it.
    if len(candidates) == 1:
        return finish(
            candidates[0], "resolved",
            "Resolved from one current symbol/alias search match.",
        )

    # Several hits: prefer a unique exact official symbol, then a unique alias.
    summaries = esummary_genes(candidates)

    if any(gid not in summaries for gid in candidates):
        return finish(
            None, "failed",
            "Could not obtain usable summaries for every candidate.",
        )

    exact = []
    alias_matches = []

    for gene_id in candidates:
        record = summaries[gene_id]
        official = record.get("nomenclaturesymbol") or record.get("name") or ""

        if str(official).casefold() == symbol.casefold():
            exact.append(gene_id)

        aliases = str(record.get("otheraliases") or "").split(",")

        if any(
            alias.strip().casefold() == symbol.casefold()
            for alias in aliases
        ):
            alias_matches.append(gene_id)

    if len(exact) == 1:
        return finish(
            exact[0], "resolved",
            "Selected the unique exact official-symbol match.",
        )

    if not exact and len(alias_matches) == 1:
        return finish(
            alias_matches[0], "resolved",
            "Selected the unique exact alias match.",
        )

    return finish(
        None, "ambiguous",
        "Several candidates remained without a unique "
        "official-symbol or alias match.",
    )


######## # Source: both original scripts; retained from the SRT base.
######## # Retrieve structured annotations for candidate Gene IDs, with JSON validation.
######## # NEW (improvement): a batch lookup is used to compare ambiguous search hits.
# Source: SRT base (both original scripts). Extended to a batch lookup so that
# several ESearch candidates can be compared in a single request.
def esummary_genes(gene_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    """Return {Gene ID: ESummary record} for the requested IDs, keeping only
    valid, current (live, not replaced) records. Returns {} on failure."""
    if not isinstance(gene_ids, (list, tuple)) or not gene_ids:
        raise AnnotationError("Supply a nonempty list of NCBI Gene IDs.")
    gene_ids = [validate_gene_id(gid) for gid in gene_ids]
    params = {"db": "gene", "id": ",".join(gene_ids), "retmode": "json"}
    response = ncbi_get("esummary.fcgi", params)
    if response is None:
        return {}
    try:
        data = response.json()
    except ValueError:
        report("[ERROR] Invalid NCBI ESummary JSON.", logging.ERROR)
        return {}
    if not isinstance(data, dict) or not isinstance(data.get("result"), dict):
        report("[ERROR] Unexpected NCBI ESummary response structure.", logging.ERROR)
        return {}
    result = data["result"]
    return {
        gid: record for gid in gene_ids
        if isinstance(record := result.get(gid), dict) and record and not record.get("error")
        # Status '' or '0' means live. A nonzero currentid marks a record that
        # has been replaced, so it is not treated as a current match.
        and str(record.get("status") or "") in {"", "0"}
        and str(record.get("currentid") or "") in {"", "0"}
    }


def esummary_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    """Return the ESummary record for one Gene ID, or None if unavailable."""
    gene_id = validate_gene_id(gene_id)
    record = esummary_genes([gene_id]).get(gene_id)
    if record is None:
        report(f"[WARN] No usable ESummary for '{gene_id}'.")
    return record


######## # Source: both original scripts; retained from the SRT base.
######## # Download the full NCBI Gene XML record so its database links can be inspected.
######## def efetch_gene_xml(gene_id: str) -> Optional[str]:
########     gene_id = validate_gene_id(gene_id)
########     # NEW (fix): Gene XML uses retmode=xml; rettype=gp is for protein records.
########     # Reference: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
########     params = {"db": "gene", "id": gene_id, "retmode": "xml"}
########     response = ncbi_get("efetch.fcgi", params)
########     return None if response is None else response.text
# Source: SRT base (both original scripts). Now cached, so a record fetched to
# verify an ENSG input is reused rather than downloaded a second time.
def efetch_gene_xml(gene_id: str) -> Optional[str]:
    """Return the full NCBI Gene XML record (used for its database
    cross-references), or None if the request failed."""
    gene_id = validate_gene_id(gene_id)
    if gene_id not in _GENE_XML_CACHE:
        # Gene records use retmode=xml with no rettype ('gp' applies to proteins).
        # Reference: https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
        response = ncbi_get("efetch.fcgi", {"db": "gene", "id": gene_id, "retmode": "xml"})
        if response is None:
            return None
        _GENE_XML_CACHE[gene_id] = response.text
    return _GENE_XML_CACHE[gene_id]


######## # Source: SRT XML parser, with its return type hint from the reference script.
######## # Collect unique ENSG identifiers from Ensembl database links in NCBI XML.
######## # ENSG identifiers are human gene IDs; other species use different prefixes.
# Source: SRT XML parser; return type hint from the reference script.
# Note: the ENSG prefix is specific to human genes (mouse genes use ENSMUSG).
def extract_ensembl_gene_ids_from_xml(
    xml_text: str,
    #### *,
    #### strict: bool = False,    # Never set by any caller.
) -> List[str]:
    """Return the unique ENSG IDs found in Ensembl Dbtag cross-references of an
    NCBI Gene XML record, in order of appearance ([] if the XML is invalid)."""
    ensg_ids = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        report(
            "[ERROR] Invalid NCBI Gene XML; Ensembl IDs left empty.",
            logging.ERROR,
        )
        #### if strict:
        ####     raise
        return ensg_ids

    ######## for dbtag in root.iter():
    ########     if not dbtag.tag.endswith("Dbtag"):
    ########         continue
    ########
    ########     db_name = None
    ########
    ########     for child in dbtag:
    ########         if child.tag.endswith("Dbtag_db"):
    ########             db_name = (child.text or "").strip()
    ########
    ########     if not db_name or "ensembl" not in db_name.lower():
    ########         continue
    ########
    ########     for child in dbtag:
    ########         if child.tag.endswith("Dbtag_tag"):
    ########             for sub in child.iter():
    ########                 ensembl_id = normalize_ensembl_gene_id(sub.text or "")
    ########
    ########                 if ensembl_id and ensembl_id not in ensg_ids:
    ########                     ensg_ids.append(ensembl_id)
    for dbtag in root.iter():
        if not dbtag.tag.endswith("Dbtag"):
            continue
        children = list(dbtag)
        # Keep only Dbtags whose database name mentions Ensembl.
        if not any(c.tag.endswith("Dbtag_db") and "ensembl" in (c.text or "").lower() for c in children):
            continue
        for child in (c for c in children if c.tag.endswith("Dbtag_tag")):
            for sub in child.iter():
                ensembl_id = normalize_ensembl_gene_id(sub.text or "")
                if ensembl_id and ensembl_id not in ensg_ids:
                    ensg_ids.append(ensembl_id)

    return ensg_ids


# ---------------------------------------------------------------------------
# Processing a gene list
# ---------------------------------------------------------------------------

######## # Source: both scripts; keep the SRT sequence of ESearch, ESummary and EFetch.
######## # Process each gene and retain the same result keys for later project stages.
######## # Validate the complete input list before starting any requests.
######## # NEW (improvement): validate once, trim spaces and remove repeated symbols.
def clean_gene_list(gene_symbols: List[str]) -> List[str]:
    """Validate the whole list before any request is made; trim whitespace and
    drop repeated symbols (case-insensitive), keeping the first occurrence."""
    if not isinstance(gene_symbols, (list, tuple)) or not gene_symbols:
        raise AnnotationError("Supply a nonempty list of gene symbols.")
    cleaned = []
    seen = set()
    for symbol in gene_symbols:
        symbol = validate_gene_symbol(symbol)
        if symbol.casefold() not in seen:
            seen.add(symbol.casefold())
            cleaned.append(symbol)
    return cleaned


# Source: both scripts; keeps the SRT sequence of ESearch, ESummary and EFetch
# and the same result keys for later project stages.
def fetch_ncbi_info_for_gene_list(
    gene_symbols: List[str],
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    """Annotate each input and return one result row per input, including
    inputs that could not be resolved (with the reason recorded)."""
    results = []

    for symbol in clean_gene_list(gene_symbols):
        print(f"\n=== Processing gene symbol: {symbol} ===")

        # 1. Resolve the input; esearch_gene records the status and reason.
        details: Dict[str, Any] = {}
        gene_id = esearch_gene(symbol, organism=organism, resolution_details=details)

        # 2. Every input gets a row with the same set of keys.
        entry = {
            "query_symbol": symbol,
            "gene_id": gene_id or "",
            "summary": {},
            "ensembl_gene_ids": [],
            "resolution_status": details.get("status", "failed"),
            "resolution_note": details.get("note", ""),
            "same_gene_note": "",
        }

        # 3. Unresolved input: keep the row with its reason and continue.
        if gene_id is None:
            report(f"[WARN] '{symbol}': {entry['resolution_status']}. {entry['resolution_note']}")
            results.append(entry)
            continue

        print(f"  -> Resolved to NCBI Gene ID: {gene_id}")

        # 4. Gene ID found but no usable summary: keep the row, marked as failed.
        summary = esummary_gene(gene_id)
        if summary is None:
            entry["resolution_status"] = "failed"
            entry["resolution_note"] = "Gene ID found but ESummary failed."
            results.append(entry)
            continue
        entry["summary"] = summary

        official = str(summary.get("nomenclaturesymbol") or summary.get("name") or "")
        if official and official.casefold() != symbol.casefold():
            report(f"[WARN] Queried '{symbol}'; the resolved official symbol is '{official}'.")

        # 5. Ensembl gene IDs from the Gene XML cross-references.
        xml_text = efetch_gene_xml(gene_id)
        if xml_text is None:
            report("[WARN] XML retrieval failed; ENSG IDs left empty.")
        else:
            entry["ensembl_gene_ids"] = extract_ensembl_gene_ids_from_xml(xml_text)
        print(f"  -> Ensembl Gene ID(s): {';'.join(entry['ensembl_gene_ids']) or 'none retrieved'}")

        results.append(entry)

    # 6. Flag inputs that resolve to a gene already listed (e.g. MYC and c-Myc).
    ######## first_query_by_gene_id: Dict[str, str] = {}
    ######## for entry in results:
    ########     gid = str(entry["gene_id"])
    ########     if not gid:
    ########         continue
    ########     if gid in first_query_by_gene_id:
    ########         entry["same_gene_note"] = (
    ########             f"Resolves to the same gene as '{first_query_by_gene_id[gid]}' (Gene ID {gid})."
    ########         )
    ########     else:
    ########         first_query_by_gene_id[gid] = entry["query_symbol"]
    first_query_by_gene_id: Dict[str, str] = {}
    for entry in results:
        gid = entry["gene_id"]
        first = first_query_by_gene_id.setdefault(gid, entry["query_symbol"]) if gid else None
        if first and first != entry["query_symbol"]:
            entry["same_gene_note"] = f"Resolves to the same gene as '{first}' (Gene ID {gid})."

    return results


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

# FIX: this function had been removed but is still called by
# connect_ncbi_gene_info, which raised a NameError. Restored from the SRT base,
# with the same name fallbacks as the TSV, the resolution status, and handling
# of unresolved rows (empty summary).
def print_human_readable_summary(results: List[Dict[str, Any]]) -> None:
    """Print a concise, human-readable summary of every result row."""
    print("\n\n===== NCBI GENE SUMMARY (human-readable) =====")
    for entry in results:
        s = entry.get("summary") or {}
        print(f"Symbol queried   : {entry['query_symbol']}")
        print(f"  Resolution     : {entry.get('resolution_status', 'NA')} - {entry.get('resolution_note') or 'NA'}")
        if not s:
            print()
            continue
        print(f"  NCBI Gene ID   : {entry['gene_id']}")
        print(f"  Official Symbol: {s.get('nomenclaturesymbol') or s.get('name') or 'NA'}")
        print(f"  Official Name  : {s.get('nomenclaturename') or s.get('description') or 'NA'}")
        print(f"  Chromosome     : {s.get('chromosome') or 'NA'}")
        print(f"  Map Location   : {s.get('maplocation') or 'NA'}")
        print(f"  Ensembl Gene ID(s): {';'.join(entry.get('ensembl_gene_ids', [])) or 'NA'}")
        if entry.get("same_gene_note"):
            print(f"  Note           : {entry['same_gene_note']}")
        print(f"  Summary        : {s.get('summary') or 'NA'}")
        print()


######## # Source: both scripts; retain the SRT table columns and their order.
######## # Write a TSV table for later database loading; multiple ENSG IDs use semicolons.
######## # New fix: apply the same name fallbacks used in the console summary.
# Source: both scripts; keeps the SRT column order, with the resolution columns
# added at the end. Name fields use the same fallbacks as the console summary.
def write_table_to_file(
    results: List[Dict[str, Any]],
    filename: str,
    delimiter: str = "\t",
) -> None:
    """Write one row per input to a delimited file (tab by default) for loading
    into MySQL. Multiple ENSG IDs are separated by semicolons."""
    header_cols = [
        "Query_Symbol", "Gene_ID", "ENSG_IDs", "Official_Symbol",
        "Official_Name", "Chromosome", "Map_Location",
        "Other_Aliases", "Summary", "Same_Gene_Note", "Resolution_Status", "Resolution_Note",
    ]

    ######## #Source: both scripts. Remove embedded line breaks and delimiters from cells.
    # Remove line breaks and delimiter characters inside a cell so that each
    # record stays on one line with the correct number of columns.
    def clean(text: Any) -> str:
        if text is None:
            return ""
        return str(text).replace("\n", " ").replace("\r", " ").replace(delimiter, " ")

    with open(filename, "w", encoding="utf-8") as f:
        f.write(delimiter.join(header_cols) + "\n")

        for entry in results:
            s = entry["summary"]
            row = [
                entry.get("query_symbol"),
                entry.get("gene_id"),
                ";".join(entry.get("ensembl_gene_ids", [])),
                s.get("nomenclaturesymbol") or s.get("name"),
                s.get("nomenclaturename") or s.get("description"),
                s.get("chromosome"),
                s.get("maplocation"),
                s.get("otheraliases"),
                s.get("summary"),
                ######## entry.get("Same_Gene_Note")
                # FIX: the key is lower case ("Same_Gene_Note" always gave an
                # empty cell), and the two resolution columns were missing, so
                # rows had fewer columns than the header.
                entry.get("same_gene_note"),
                entry.get("resolution_status"),
                entry.get("resolution_note"),
            ]
            f.write(delimiter.join(clean(value) for value in row) + "\n")

    print(f"[INFO] Wrote table to {os.path.abspath(filename)}")


######## # Source: both original scripts; retained from the SRT base.
######## # Save each complete ESummary record and its extracted Ensembl IDs as JSON.
# Source: SRT base (both original scripts).
def save_results_to_json(results: List[Dict[str, Any]], filename: str) -> None:
    """Save the full result rows, including complete ESummary records, as JSON."""
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"[INFO] Saved JSON results to {os.path.abspath(filename)}")


######## # Source: p1_basic_anno.py (SRT base); this reusable function is retained.
######## # Fetch annotations, print them, save TSV/JSON files and return the records.
######## # Report missing annotations and missing ENSG IDs, then save successful records.
# Source: p1_basic_anno.py (SRT base); reusable entry point for other scripts.
def connect_ncbi_gene_info(
    gene_list: List[str],
    output_filename: str,
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    """Annotate gene_list, print a summary, write <output_filename>.tsv and
    .json, and return the result rows.

    Raises AnnotationError, leaving existing output files untouched, if no
    input could be resolved."""
    requested = clean_gene_list(gene_list)
    results = fetch_ncbi_info_for_gene_list(requested, organism=organism)
    ######## print(f"[INFO] Retrieved {len(results)} of {len(requested)} requested gene entries.")
    ######## retrieved = {entry["query_symbol"] for entry in results}
    ######## missing = [symbol for symbol in requested if symbol not in retrieved]
    # FIX: unresolved inputs are now kept as rows, so len(results) always equals
    # len(requested) and the old 'missing' check could never find anything.
    # Count rows by resolution status instead.
    missing = [entry["query_symbol"] for entry in results if entry["resolution_status"] != "resolved"]
    print(f"[INFO] Retrieved {len(results) - len(missing)} of {len(requested)} requested gene entries.")
    if missing:
        report(f"[WARN] No annotation retrieved for: {', '.join(missing)}")
    ######## no_ensg = [entry["query_symbol"] for entry in results if not entry["ensembl_gene_ids"]]
    # FIX: only resolved genes can have Ensembl IDs, so list only those here.
    no_ensg = [
        entry["query_symbol"] for entry in results
        if entry["resolution_status"] == "resolved" and not entry["ensembl_gene_ids"]
    ]
    if no_ensg:
        report(f"[WARN] No Ensembl Gene ID retrieved for: {', '.join(no_ensg)}")
    ######## if not results:
    # FIX: results is never empty now; stop only if no input was resolved.
    if len(missing) == len(results):
        raise AnnotationError("No gene annotations were retrieved; existing output files were not replaced.")
    ######## # NEW (fix): create the destination even when this wrapper is imported.
    # Create the output folder here so the function also works when imported.
    output_filename = os.path.expanduser(os.fspath(output_filename))
    parent = os.path.dirname(os.path.abspath(output_filename))
    os.makedirs(parent, exist_ok=True)
    print_human_readable_summary(results)
    write_table_to_file(results, output_filename + ".tsv")
    save_results_to_json(results, output_filename + ".json")
    return results


######## # Source: main() adapted from p1_basic_anno_ref.py.
######## # Adaptation: call the SRT wrapper instead of repeating its fetch/save steps.
######## # Edit gene_list for your assignment; output_filename is a prefix, without .tsv.
######## # Data files are saved under ~/BDB_ICA/OP1_data/development.
######## # Source: main() adapted from p1_basic_anno_ref.py.
######## # The ten development genes come from the assignment brief.
######## # Use the SRT wrapper to retrieve human annotations and save TSV/JSON files.
# Source: main() adapted from p1_basic_anno_ref.py; calls the SRT wrapper
# rather than repeating its fetch and save steps.
# Runs the development list first and checks it; the challenging list is only
# run if every development gene was fully annotated.
# Configuration: edit the gene lists, organism and output folder below.
def main() -> None:
    # 1. Input lists, both taken from the assignment brief.
    development_genes = [
        "BRCA1", "BRCA2", "TP53", "EGFR", "MYC",
        "APOE", "CFTR", "ACE2", "TNF", "VEGFA",
    ]

    challenge_genes = [
        "MYC", "c-Myc", "Blah", "ENSG00000136997", "AADACL2",
        "Plet1", "HYDIN", "ENSG00000261846", "XIST", "MALAT1",
    ]

    organism = "Homo sapiens"

    # 2. Output prefixes (without .tsv/.json); each list has its own folder.
    base_dir = os.path.expanduser(
        "~/BDB_ICA/OP1_data/final_script_run_for_both_lists"
    )

    development_prefix = os.path.join(
        base_dir, "development", "ncbi_gene_summary2"
    )

    challenge_prefix = os.path.join(
        base_dir, "challenge", "ncbi_gene_summary2"
    )

    # 3. Development run: annotate and save the TSV and JSON files.
    print("[INFO] Starting development run.")

    development_results = connect_ncbi_gene_info(
        development_genes,
        development_prefix,
        organism=organism,)

    #### # 4. Check that every expected input has exactly one row, in order.
    #### # Not needed: fetch_ncbi_info_for_gene_list builds exactly one row per
    #### # cleaned input, in input order, so this check can never fail.
    # development completeness checks,
    expected_inputs = clean_gene_list(development_genes)
    returned_inputs = [
        entry.get("query_symbol")
        for entry in development_results]
    ####
    if returned_inputs != expected_inputs:
        raise AnnotationError(
        "Development results do not contain the expected input rows "
        "in order; challenge run was not started.")

    ######## # 5. Check the contents of the development rows.
    ######## problems = []
    ########
    ######## for entry in development_results:
    ########     reasons = []
    ########
    ########     if entry.get("resolution_status") != "resolved":
    ########         reasons.append(
    ########             f"status={entry.get('resolution_status', 'missing')}"
    ########         )
    ########
    ########     if not entry.get("gene_id"):
    ########         reasons.append("missing NCBI Gene ID")
    ########
    ########     if not entry.get("summary"):
    ########         reasons.append("missing annotation summary")
    ########
    ########     if not entry.get("ensembl_gene_ids"):
    ########         reasons.append("missing Ensembl Gene ID")
    ########
    ########     if reasons:
    ########         problems.append(
    ########             f"{entry['query_symbol']}: {', '.join(reasons)}"
    ########         )
    # 4. Check every development row. A row is only marked "resolved" once it
    # has both a Gene ID and a summary, so the status check covers both.
    problems = []
    for entry in development_results:
        reasons = [] if entry["resolution_status"] == "resolved" else [f"status={entry['resolution_status']}"]
        if not entry["ensembl_gene_ids"]:
            reasons.append("missing Ensembl Gene ID")
        if reasons:
            problems.append(f"{entry['query_symbol']}: {', '.join(reasons)}")

    # 5. Stop if any development check failed.
    if problems:
        raise AnnotationError(
            "Development checks failed; challenge run was not started. "
            + " | ".join(problems)
        )

    # 6. Challenge run, started automatically once development passes.
    print("[INFO] Development checks passed. Starting challenge run.")

    connect_ncbi_gene_info(
        challenge_genes,
        challenge_prefix,
        organism=organism,
    )

    print("[INFO] Challenge results saved.")


if __name__ == "__main__":
    try:
        # Only safe messages are logged. Exception text can contain a request
        # URL (and so the API key), so only the exception type and code
        # locations are recorded.
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
        message = str(error) if isinstance(error, AnnotationError) else "Run stopped. Check the input, connection and ncbi_annotation.log."
        report(f"[ERROR] {message}", logging.ERROR)
        sys.exit(1)
