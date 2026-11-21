#!/usr/bin/env python3
# Source: p1_basic_anno.py (SRT base), combined with p1_basic_anno_ref.py.
# Purpose: retrieve human gene annotations and Ensembl gene IDs from NCBI.
# Load tools for HTTP requests, JSON files, XML parsing, delays and type hints.
import sys
import time
import json
from typing import List, Dict, Any, Optional
import xml.etree.ElementTree as ET

import requests


# Source: p1_basic_anno.py (SRT base), with your email substituted.
# Define the NCBI service, contact details and a five-second request delay.
# An API key is optional; leaving it as None is fine for this script.
SLEEPBETWEEN = 5.0
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "s2924525@ed.ac.uk"
NCBI_API_KEY: Optional[str] = None


# Source: SRT request helper, with type hints and three attempts from the ref.
# New fix: pause before EVERY request, including successful requests.
# Send a request, retry failures and return None if every attempt fails.
def ncbi_get(
    endpoint: str,
    params: Dict[str, Any],
    max_retries: int = 3,
    sleep_between: float = SLEEPBETWEEN,
) -> Optional[requests.Response]:
    """Make an NCBI request; max_retries counts total attempts, including the first."""
    url = f"{NCBI_BASE}/{endpoint}"
    params = dict(params)

    if NCBI_EMAIL:
        params.setdefault("email", NCBI_EMAIL)
    if NCBI_API_KEY:
        params.setdefault("api_key", NCBI_API_KEY)

    for attempt in range(1, max_retries + 1):
        time.sleep(sleep_between)

        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.status_code == 200:
                return resp

            print(
                f"[WARN] Attempt {attempt}/{max_retries} for {endpoint} "
                f"returned HTTP {resp.status_code}.",
                file=sys.stderr,
            )
        except requests.exceptions.RequestException as error:
            # New fix: report the error type without exposing URLs or API keys.
            print(
                f"[ERROR] Attempt {attempt}/{max_retries} for {endpoint} "
                f"failed: {type(error).__name__}.",
                file=sys.stderr,
            )

    print(
        f"[ERROR] Failed to retrieve {endpoint} after {max_retries} attempts.",
        file=sys.stderr,
    )
    return None


# Source: SRT ESearch function, with type hints from p1_basic_anno_ref.py.
# Find the first NCBI Gene ID matching a gene symbol and organism.
# New additions: handle invalid JSON with ValueError and warn about ambiguity.
def esearch_gene(
    symbol: str,
    organism: str = "Homo sapiens",
) -> Optional[str]:
    term = f'{symbol}[Gene Name] AND "{organism}"[Organism]'
    params = {"db": "gene", "term": term, "retmode": "json"}
    resp = ncbi_get("esearch.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError:
        print(f"[ERROR] Invalid ESearch JSON for '{symbol}'.", file=sys.stderr)
        return None

    id_list = data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        print(
            f"[WARN] No Gene ID found for '{symbol}' in '{organism}'.",
            file=sys.stderr,
        )
        return None

    if len(id_list) > 1:
        print(
            f"[WARN] Multiple matches for '{symbol}'; using Gene ID "
            f"{id_list[0]}. Check this matches your intended gene.",
            file=sys.stderr,
        )

    return id_list[0]


# Source: both original scripts; retained from the SRT base.
# Retrieve structured annotation for a Gene ID, including names and chromosome.
# New fix: catch ValueError to handle invalid JSON returned by requests.
def esummary_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    params = {"db": "gene", "id": gene_id, "retmode": "json"}
    resp = ncbi_get("esummary.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError:
        print(f"[ERROR] Invalid ESummary JSON for '{gene_id}'.", file=sys.stderr)
        return None

    result = data.get("result", {})
    if gene_id not in result or result[gene_id].get("error"):
        print(f"[WARN] No usable ESummary for '{gene_id}'.", file=sys.stderr)
        return None

    return result[gene_id]


# Source: both original scripts; retained from the SRT base.
# Download the full NCBI Gene XML record so its database links can be inspected.
def efetch_gene_xml(gene_id: str) -> Optional[str]:
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
    except ET.ParseError as error:
        print(f"[ERROR] XML parse error: {error}", file=sys.stderr)
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
# New addition: remove surrounding spaces and skip empty gene symbols.
def fetch_ncbi_info_for_gene_list(
    gene_symbols: List[str],
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    results = []

    for symbol in gene_symbols:
        symbol = symbol.strip()
        if not symbol:
            continue

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
            print("[WARN] XML retrieval failed; ENSG IDs left empty.", file=sys.stderr)
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
    print(f"[INFO] Retrieved {len(results)} of {len(gene_list)} requested gene entries.")
    print_human_readable_summary(results)
    write_table_to_file(results, output_filename + ".tsv")
    save_results_to_json(results, output_filename + ".json")
    return results


# Source: main() and the three example genes from p1_basic_anno_ref.py.
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
        "VEGFA",]

    organism = "Homo sapiens"
   output_dir = os.path.expanduser("~/BDB_ICA/OP1_data/development")
      output_filename = os.path.join(output_dir, "ncbi_gene_summary2_try1")

    connect_ncbi_gene_info(
        gene_list,
        output_filename,
        organism=organism,)


# Source: p1_basic_anno_ref.py.
# Run main() when this file is executed; importing it does not run the workflow.
if __name__ == "__main__":
    main()