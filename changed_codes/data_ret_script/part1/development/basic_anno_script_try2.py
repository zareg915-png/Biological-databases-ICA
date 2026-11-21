#!/usr/bin/env python3
# Source: p1_basic_anno.py (SRT base), combined with p1_basic_anno_ref.py.
# Purpose: retrieve human gene annotations and Ensembl gene IDs from NCBI. This version runs as a single script.
# GenAI assistance: ChatGPT helped edit. 
# Load tools for HTTP requests, JSON files, XML parsing, delays and type hints.
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


# Source: p1_basic_anno.py (SRT base), with runtime contact settings added.
# Define the NCBI service, contact details and a five-second request delay.
# Set NCBI_EMAIL in your environment before running; NCBI_API_KEY is optional.
SLEEPBETWEEN = 5.0
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
NCBI_API_KEY: Optional[str] = os.environ.get("NCBI_API_KEY", "").strip() or None
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}
LOG_FILE = "ncbi_annotation.log"
LOGGER = logging.getLogger("ncbi_annotation")
LOGGER.addHandler(logging.NullHandler())
LOGGER.propagate = False


class AnnotationError(ValueError):
    """A safe explanation of invalid input or an unsuccessful annotation run."""


# New: validate inputs before making requests. Do not accept query expressions.
def validate_gene_symbol(symbol: str) -> str:
    if not isinstance(symbol, str):
        raise AnnotationError("Each gene symbol must be text.")
    symbol = symbol.strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", symbol):
        raise AnnotationError("Gene symbols must be nonempty and use letters, digits, dots, underscores or hyphens.")
    return symbol


def validate_gene_id(gene_id: str) -> str:
    gene_id = str(gene_id)
    if not re.fullmatch(r"[1-9][0-9]*", gene_id):
        raise AnnotationError("An NCBI Gene ID must be a positive whole number.")
    return gene_id


# New: send the same safe message to the console and diagnostic log.
# Never pass raw exception text, response bodies, URLs or credentials here.
def report(message: str, level: int = logging.WARNING) -> None:
    LOGGER.log(level, message)
    print(message, file=sys.stderr)


def retry_after_seconds(response: requests.Response) -> float:
    value = response.headers.get("Retry-After")
    if value is None:
        return 0.0
    try:
        seconds = float(value)
        return max(0.0, seconds) if math.isfinite(seconds) else 0.0
    except ValueError:
        try:
            when = parsedate_to_datetime(value)
            if when.tzinfo is None:
                when = when.replace(tzinfo=timezone.utc)
            return max(0.0, (when - datetime.now(timezone.utc)).total_seconds())
        except (TypeError, ValueError, OverflowError):
            return 0.0


# Source: SRT request helper, with type hints and three attempts from the ref.
# New fix: pause before EVERY request, including successful requests.
# Retry only temporary read failures, using bounded exponential backoff.
def ncbi_get(
    endpoint: str,
    params: Dict[str, Any],
    max_retries: int = 3,
    sleep_between: float = SLEEPBETWEEN,
) -> Optional[requests.Response]:
    """Make an NCBI request; max_retries counts total attempts, including the first."""
    if endpoint not in {"esearch.fcgi", "esummary.fcgi", "efetch.fcgi"}:
        raise AnnotationError("Unsupported NCBI endpoint.")
    if not isinstance(params, dict):
        raise AnnotationError("NCBI request parameters must be a dictionary.")
    if not isinstance(max_retries, int) or not 1 <= max_retries <= 5:
        raise AnnotationError("The number of request attempts must be between 1 and 5.")
    if not isinstance(sleep_between, (int, float)) or not math.isfinite(sleep_between) or sleep_between < 0:
        raise AnnotationError("The request delay must be finite and nonnegative.")
    url = f"{NCBI_BASE}/{endpoint}"
    params = dict(params)

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
            # Extra waits grow 1, 2, 4, ... seconds; the five-second request
            # delay above still applies before each attempt by default.
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


# Source: SRT ESearch function, with type hints from p1_basic_anno_ref.py.
# Find the first NCBI Gene ID matching a gene symbol and organism.
# New additions: handle invalid JSON with ValueError and warn about ambiguity.
def esearch_gene(
    symbol: str,
    organism: str = "Homo sapiens",
) -> Optional[str]:
    symbol = validate_gene_symbol(symbol)
    if not isinstance(organism, str) or not organism.strip() or any(c in organism for c in '"[]\r\n'):
        raise AnnotationError("Supply a valid organism name.")
    organism = organism.strip()
    term = f'{symbol}[Gene Name] AND "{organism}"[Organism]'
    params = {"db": "gene", "term": term, "retmode": "json"}
    resp = ncbi_get("esearch.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError:
        report(f"[ERROR] Invalid ESearch JSON for '{symbol}'.", logging.ERROR)
        return None

    if not isinstance(data, dict) or not isinstance(data.get("esearchresult"), dict):
        report("[ERROR] Unexpected NCBI ESearch response structure.", logging.ERROR)
        return None
    search_result = data["esearchresult"]
    if data.get("error") or search_result.get("errorlist"):
        report("[ERROR] NCBI ESearch reported a query error.", logging.ERROR)
        return None
    id_list = search_result.get("idlist")
    if not isinstance(id_list, list) or any(not isinstance(gid, str) or not re.fullmatch(r"[1-9][0-9]*", gid) for gid in id_list):
        report("[ERROR] NCBI ESearch returned an invalid Gene ID list.", logging.ERROR)
        return None
    if not id_list:
        report(
            f"[WARN] No Gene ID found for '{symbol}' in '{organism}'."
        )
        return None

    if len(id_list) > 1:
        report(
            f"[WARN] Multiple matches for '{symbol}'; using Gene ID "
            f"{id_list[0]}. Check this matches your intended gene."
        )

    return id_list[0]


# Source: both original scripts; retained from the SRT base.
# Retrieve structured annotation for a Gene ID, including names and chromosome.
# New fix: catch ValueError to handle invalid JSON returned by requests.
def esummary_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    gene_id = validate_gene_id(gene_id)
    params = {"db": "gene", "id": gene_id, "retmode": "json"}
    resp = ncbi_get("esummary.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError:
        report(f"[ERROR] Invalid ESummary JSON for '{gene_id}'.", logging.ERROR)
        return None

    if not isinstance(data, dict) or not isinstance(data.get("result"), dict):
        report("[ERROR] Unexpected NCBI ESummary response structure.", logging.ERROR)
        return None
    result = data["result"]
    summary = result.get(gene_id)
    if not isinstance(summary, dict) or not summary or summary.get("error"):
        report(f"[WARN] No usable ESummary for '{gene_id}'.")
        return None

    return summary


# Source: both original scripts; retained from the SRT base.
# Download the full NCBI Gene XML record so its database links can be inspected.
def efetch_gene_xml(gene_id: str) -> Optional[str]:
    gene_id = validate_gene_id(gene_id)
    params = {
        "db": "gene",
        "id": gene_id,
        "rettype": "gp",
        "retmode": "xml",
    }
    resp = ncbi_get("efetch.fcgi", params)
    if resp is None:
        return None
    return resp.text


# Source: SRT XML parser, with its return type hint from the reference script.
# Collect unique ENSG identifiers from Ensembl database links in NCBI XML.
# ENSG identifiers are human gene IDs; other species use different prefixes.
def extract_ensembl_gene_ids_from_xml(xml_text: str) -> List[str]:
    ensg_ids: List[str] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        report("[ERROR] Invalid NCBI Gene XML; Ensembl IDs left empty.", logging.ERROR)
        return ensg_ids

    # Source: both scripts. Recognise Dbtag elements, including namespaced tags.
    def iter_dbtag_elements(elem):
        for dbtag in elem.iter():
            if dbtag.tag.endswith("Dbtag"):
                yield dbtag

    for dbtag in iter_dbtag_elements(root):
        db_name = None
        for child in dbtag:
            if child.tag.endswith("Dbtag_db"):
                db_name = (child.text or "").strip()

        if not db_name or "ensembl" not in db_name.lower():
            continue

        for child in dbtag:
            if child.tag.endswith("Dbtag_tag"):
                for sub in child.iter():
                    text = (sub.text or "").strip()
                    if text.startswith("ENSG") and text not in ensg_ids:
                        ensg_ids.append(text)

    return ensg_ids


# Source: both scripts; keep the SRT sequence of ESearch, ESummary and EFetch.
# Process each gene and retain the same result keys for later project stages.
# Validate the complete input list before starting any requests.
def fetch_ncbi_info_for_gene_list(
    gene_symbols: List[str],
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    if not isinstance(gene_symbols, (list, tuple)) or not gene_symbols:
        raise AnnotationError("Supply a nonempty list of gene symbols.")
    gene_symbols = [validate_gene_symbol(symbol) for symbol in gene_symbols]
    results = []

    for symbol in gene_symbols:

        print(f"\n=== Processing gene symbol: {symbol} ===")
        gene_id = esearch_gene(symbol, organism=organism)
        if gene_id is None:
            continue

        print(f"  -> Resolved to NCBI Gene ID: {gene_id}")
        summary = esummary_gene(gene_id)
        if summary is None:
            continue

        xml_text = efetch_gene_xml(gene_id)
        if xml_text is None:
            report("[WARN] XML retrieval failed; ENSG IDs left empty.")
            ensg_ids = []
        else:
            ensg_ids = extract_ensembl_gene_ids_from_xml(xml_text)

        print(f"  -> Ensembl Gene ID(s): {';'.join(ensg_ids) or 'none retrieved'}")
        results.append({
            "query_symbol": symbol,
            "gene_id": gene_id,
            "summary": summary,
            "ensembl_gene_ids": ensg_ids,
        })

    return results


# Source: both scripts; retain the SRT console summary.
# Display the annotations in a readable form.
# New fix: use NCBI name/description when official nomenclature fields are empty.
def print_human_readable_summary(results: List[Dict[str, Any]]) -> None:
    print("\n===== NCBI GENE SUMMARY =====")

    for entry in results:
        s = entry["summary"]
        official_symbol = s.get("nomenclaturesymbol") or s.get("name") or "NA"
        official_name = s.get("nomenclaturename") or s.get("description") or "NA"
        ensg_str = ";".join(entry.get("ensembl_gene_ids", [])) or "NA"

        print(f"Symbol queried   : {entry['query_symbol']}")
        print(f"  NCBI Gene ID   : {entry['gene_id']}")
        print(f"  Official Symbol: {official_symbol}")
        print(f"  Official Name  : {official_name}")
        print(f"  Chromosome     : {s.get('chromosome') or 'NA'}")
        print(f"  Map Location   : {s.get('maplocation') or 'NA'}")
        print(f"  Ensembl Gene ID(s): {ensg_str}")
        print(f"  Summary        : {s.get('summary') or 'NA'}")
        print()


# Source: both scripts; retain the SRT table columns and their order.
# Write a TSV table for later database loading; multiple ENSG IDs use semicolons.
# New fix: apply the same name fallbacks used in the console summary.
def write_table_to_file(
    results: List[Dict[str, Any]],
    filename: str,
    delimiter: str = "\t",
) -> None:
    header_cols = [
        "Query_Symbol", "Gene_ID", "ENSG_IDs", "Official_Symbol",
        "Official_Name", "Status", "Chromosome", "Map_Location",
        "Other_Aliases", "Summary",
    ]

    # Source: both scripts. Remove embedded line breaks and delimiters from cells.
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
                s.get("status"),
                s.get("chromosome"),
                s.get("maplocation"),
                s.get("otheraliases"),
                s.get("summary"),
            ]
            f.write(delimiter.join(clean(value) for value in row) + "\n")

    print(f"[INFO] Wrote table to {filename}")


# Source: both original scripts; retained from the SRT base.
# Save each complete ESummary record and its extracted Ensembl IDs as JSON.
def save_results_to_json(results: List[Dict[str, Any]], filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"[INFO] Saved JSON results to {filename}")


# Source: p1_basic_anno.py (SRT base); this reusable function is retained.
# Fetch annotations, print them, save TSV/JSON files and return the records.
# New additions: type hints, optional organism argument and a retrieval count.
def connect_ncbi_gene_info(
    gene_list: List[str],
    output_filename: str,
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    results = fetch_ncbi_info_for_gene_list(gene_list, organism=organism)
    if not results:
        raise AnnotationError("No gene annotations were retrieved; existing output files were not replaced.")
    print(f"[INFO] Retrieved {len(results)} of {len(gene_list)} requested gene entries.")
    if len(results) < len(gene_list):
        report("[WARN] Some requested genes have no retrieved annotation; output contains the successful entries.")
    print_human_readable_summary(results)
    write_table_to_file(results, output_filename + ".tsv")
    save_results_to_json(results, output_filename + ".json")
    return results


# Source: main() adapted from p1_basic_anno_ref.py.
# Adaptation: call the SRT wrapper instead of repeating its fetch/save steps.
# Edit gene_list for your assignment; output_filename is a prefix, without .tsv.
# Files are saved in the current working directory, as in the SRT version.
# Source: main() adapted from p1_basic_anno_ref.py.
# The ten development genes were supplied by you.
# Use the SRT wrapper to retrieve human annotations and save TSV/JSON files.
def main() -> None:
    gene_list = [
        "BRCA1",
        "BRCA2",
        "TP53",
        "EGFR",
        "MYC",
        "APOE",
        "CFTR",
        "ACE2",
        "TNF",
        "VEGFA",
    ]

    organism = "Homo sapiens"
    output_dir = os.path.expanduser("~/BDB_ICA/OP1_data/development")
   output_filename = os.path.join(output_dir, "ncbi_gene_summary2_try2")

    connect_ncbi_gene_info(
        gene_list,
        output_filename,
        organism=organism,
    )


# Source: p1_basic_anno_ref.py.
# Run main() when this file is executed; importing it does not run the workflow.
if __name__ == "__main__":
    try:
        # Only safe messages are logged. Exception text can contain a request
        # URL or API key, so log its type and code locations instead.
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
