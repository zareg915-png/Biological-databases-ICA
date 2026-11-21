#!/usr/bin/env python3
"""
ncbi_gene_fetch_efetch_ensg.py

Use NCBI E-utilities (ESearch + ESummary + EFetch) to download information
for a list of genes from the NCBI Gene database, and include Ensembl Gene IDs
(ENSG...) extracted from dbXrefs in the full XML record.

- No calls to Ensembl APIs.
- No use of Plumber.
- Single-file script suitable for running in PyCharm.
"""

import sys
import time
import json
from typing import List, Dict, Any, Optional

import requests
import xml.etree.ElementTree as ET

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI recommends providing an email; replace with a real one if possible.
NCBI_EMAIL = "your_email@example.com"  # <- change this
NCBI_API_KEY: Optional[str] = None     # <- optional: "your_ncbi_api_key"


def ncbi_get(endpoint: str,
             params: Dict[str, Any],
             max_retries: int = 3,
             sleep_between: float = 0.5) -> Optional[requests.Response]:
    """
    Generic helper for calling NCBI E-utilities with basic retry logic.

    endpoint: e.g. "esearch.fcgi", "esummary.fcgi", "efetch.fcgi"
    params: dict of query parameters
    """
    url = f"{NCBI_BASE}/{endpoint}"

    params = dict(params)  # shallow copy
    if NCBI_EMAIL:
        params.setdefault("email", NCBI_EMAIL)
    if NCBI_API_KEY:
        params.setdefault("api_key", NCBI_API_KEY)

    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.status_code == 200:
                return resp
            else:
                print(f"[WARN] Attempt {attempt} for {url} returned "
                      f"status {resp.status_code}. Retrying...",
                      file=sys.stderr)
        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Attempt {attempt} for {url} raised {e}",
                  file=sys.stderr)

        time.sleep(sleep_between)

    print(f"[ERROR] Failed to retrieve {url} after {max_retries} attempts.",
          file=sys.stderr)
    return None


def esearch_gene(symbol: str,
                 organism: str = "Homo sapiens") -> Optional[str]:
    """
    Use ESearch to find a gene's NCBI Gene ID (UID) by symbol and organism.

    Query format:
      SYMBOL[Gene Name] AND "Homo sapiens"[Organism]

    Returns:
        First matching Gene ID as a string, or None if not found.
    """
    term = f"{symbol}[Gene Name] AND \"{organism}\"[Organism]"
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json"
    }
    resp = ncbi_get("esearch.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except json.JSONDecodeError:
        print(f"[ERROR] Failed to parse JSON for ESearch term '{term}'",
              file=sys.stderr)
        return None

    id_list = data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        print(f"[WARN] No Gene ID found for symbol '{symbol}' in '{organism}'.",
              file=sys.stderr)
        return None

    return id_list[0]


def esummary_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    """
    Use ESummary to get structured information about a gene by NCBI Gene ID.

    Returns:
        Dictionary with the ESummary record for that gene, or None.
    """
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json"
    }
    resp = ncbi_get("esummary.fcgi", params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except json.JSONDecodeError:
        print(f"[ERROR] Failed to parse JSON for ESummary id '{gene_id}'",
              file=sys.stderr)
        return None

    result = data.get("result", {})
    if gene_id not in result:
        print(f"[WARN] ESummary record for '{gene_id}' not in 'result'.",
              file=sys.stderr)
        return None

    return result[gene_id]


def efetch_gene_xml(gene_id: str) -> Optional[str]:
    """
    Use EFetch to get the full Gene XML record for a given NCBI Gene ID.

    According to NCBI Gene docs, use:
      db=gene, rettype=gp, retmode=xml

    Returns:
        XML string, or None if retrieval fails.
    """
    params = {
        "db": "gene",
        "id": gene_id,
        "rettype": "gp",
        "retmode": "xml"
    }
    resp = ncbi_get("efetch.fcgi", params)
    if resp is None:
        return None

    return resp.text


def extract_ensembl_gene_ids_from_xml(xml_text: str) -> List[str]:
    """
    Parse the NCBI Gene XML and extract Ensembl Gene IDs (ENSG...)
    from dbXrefs.

    In the XML, external references typically appear as:
      <Gene-commentary_dbxref>
        <Dbtag>
          <Dbtag_db>Ensembl</Dbtag_db>
          <Dbtag_tag>
            <Object-id>
              <Object-id_id>ENSG00000...</Object-id_id>
            </Object-id>
          </Dbtag_tag>
        </Dbtag>
      </Gene-commentary_dbxref>

    or sometimes as text like "Ensembl:ENSG00000...".
    This function tries to be robust to common patterns, but focuses on
    extracting strings that start with 'ENSG'.
    """
    ensg_ids: List[str] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"[ERROR] XML parse error: {e}", file=sys.stderr)
        return ensg_ids

    # Search for any Dbtag elements that might reference Ensembl
    # Depending on schema, tag names may be fully qualified or plain.
    # We'll search by local-name to be robust against namespaces.
    def iter_dbtag_elements(elem):
        for dbtag in elem.iter():
            if dbtag.tag.endswith("Dbtag"):
                yield dbtag

    for dbtag in iter_dbtag_elements(root):
        db_name = None
        # Find Dbtag_db child (database name)
        for child in dbtag:
            if child.tag.endswith("Dbtag_db"):
                db_name = (child.text or "").strip()
        if not db_name:
            continue

        # We are interested in Ensembl links
        if "ensembl" not in db_name.lower():
            continue

        # Now extract IDs from Dbtag_tag subtree
        for child in dbtag:
            if child.tag.endswith("Dbtag_tag"):
                # inside we may find Object-id or text entries
                texts = []

                for sub in child.iter():
                    if sub.text:
                        texts.append(sub.text.strip())

                for t in texts:
                    if t.startswith("ENSG"):
                        if t not in ensg_ids:
                            ensg_ids.append(t)

    return ensg_ids


def fetch_ncbi_info_for_gene_list(gene_symbols: List[str],
                                  organism: str = "Homo sapiens") -> List[Dict[str, Any]]:
    """
    For each gene symbol:
      1. ESearch (db=gene) to get Gene ID.
      2. ESummary on that Gene ID for general info.
      3. EFetch (XML) on that Gene ID, parse dbXrefs to get ENSG IDs.

    Returns:
        List of records:
        {
          "query_symbol": ...,
          "gene_id": ...,
          "summary": {...},
          "ensembl_gene_ids": ["ENSG...", ...],
        }
    """
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
            print("  -> Failed to retrieve XML via EFetch; ENSG IDs will be empty.",
                  file=sys.stderr)
            ensg_ids = []
        else:
            ensg_ids = extract_ensembl_gene_ids_from_xml(xml_text)

        if ensg_ids:
            print(f"  -> Extracted Ensembl Gene ID(s) from dbXrefs: {', '.join(ensg_ids)}")
        else:
            print("  -> No Ensembl Gene IDs found in dbXrefs of the XML record.")

        results.append({
            "query_symbol": symbol,
            "gene_id": gene_id,
            "summary": summary,
            "ensembl_gene_ids": ensg_ids,
        })

    return results


def print_human_readable_summary(results: List[Dict[str, Any]]) -> None:
    """
    Print a concise, human-readable summary to the console.
    """
    print("\n\n===== NCBI GENE SUMMARY (human-readable) =====")

    for entry in results:
        s = entry["summary"]
        ensg_str = ";".join(entry.get("ensembl_gene_ids", [])) or "NA"
        print(f"Symbol queried   : {entry['query_symbol']}")
        print(f"  NCBI Gene ID   : {entry['gene_id']}")
        print(f"  Official Symbol: {s.get('nomenclaturesymbol', s.get('name', 'NA'))}")
        print(f"  Official Name  : {s.get('nomenclaturename', 'NA')}")
        print(f"  Chromosome     : {s.get('chromosome', 'NA')}")
        print(f"  Map Location   : {s.get('maplocation', 'NA')}")
        print(f"  Ensembl Gene ID(s) from dbXrefs: {ensg_str}")
        print(f"  Summary        : {s.get('summary', 'NA')}")
        print()


def write_table_to_file(results: List[Dict[str, Any]],
                        filename: str,
                        delimiter: str = "\t") -> None:
    """
    Write a delimited summary table (default: tab-delimited) to a file.

    Columns include Ensembl gene IDs (ENSG), taken from dbXrefs in NCBI Gene XML:
      Query_Symbol, Gene_ID, ENSG_IDs, Official_Symbol, Official_Name,
      Status, Chromosome, Map_Location, Other_Aliases, Summary
    """
    header_cols = [
        "Query_Symbol",
        "Gene_ID",
        "ENSG_IDs",         # Ensembl Gene IDs from NCBI Gene XML dbXrefs
        "Official_Symbol",
        "Official_Name",
        "Status",
        "Chromosome",
        "Map_Location",
        "Other_Aliases",
        "Summary",
    ]

    def clean(text: Any) -> str:
        if text is None:
            return ""
        s = str(text)
        s = s.replace("\n", " ").replace("\r", " ")
        s = s.replace(delimiter, " ")
        return s

    with open(filename, "w", encoding="utf-8") as f:
        f.write(delimiter.join(header_cols) + "\n")

        for entry in results:
            s = entry["summary"]
            ensg_str = ";".join(entry.get("ensembl_gene_ids", []))
            row = [
                clean(entry.get("query_symbol")),
                clean(entry.get("gene_id")),
                clean(ensg_str),
                clean(s.get("nomenclaturesymbol", s.get("name"))),
                clean(s.get("nomenclaturename")),
                clean(s.get("status")),
                clean(s.get("chromosome")),
                clean(s.get("maplocation")),
                clean(s.get("otheraliases")),
                clean(s.get("summary")),
            ]
            f.write(delimiter.join(row) + "\n")

    print(f"[INFO] Wrote tab-delimited table to {filename}")


def save_results_to_json(results: List[Dict[str, Any]], filename: str) -> None:
    """
    Save full records (ESummary data + extracted ENSG IDs) for later use.
    """
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"[INFO] Saved full JSON results to {filename}")


def main():
    # Example list of 10 human gene symbols. Edit as needed.
    gene_list = [
        "BRCA1",
        "BRCA2",
        "TP53",
    ]

    # Change organism if required (e.g., "Mus musculus")
    organism = "Homo sapiens"

    # Fetch NCBI Gene info via ESearch + ESummary + EFetch XML
    results = fetch_ncbi_info_for_gene_list(
        gene_symbols=gene_list,
        organism=organism
    )

    # Human-readable summary
    print_human_readable_summary(results)

    # Tab-delimited table including ENSG IDs from dbXrefs
    table_output_file = "ncbi_gene_summary.tsv"
    write_table_to_file(results, filename=table_output_file, delimiter="\t")

    # Save raw JSON (ESummary + ENSG list)
    json_output_file = "ncbi_gene_info.json"
    save_results_to_json(results, json_output_file)


if __name__ == "__main__":
    main()
