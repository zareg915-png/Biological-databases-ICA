#!/usr/bin/env python3
# Source: p1_basic_anno.py (SRT base), combined with p1_basic_anno_ref.py.
# Purpose: retrieve human gene annotations and Ensembl gene IDs from NCBI.
#
# Markers used in this checked version:
#   "# NEW (fix)"          = a mistake in the earlier version that was corrected
#   "# NEW (improvement)"  = an addition that makes the script safer or clearer
# Everything without a NEW marker is unchanged from your version.

# Load tools for HTTP requests, JSON files, XML parsing, delays and type hints.
import os
import sys
import time
import json
from typing import List, Dict, Any, Optional
import xml.etree.ElementTree as ET

import requests


# Define the NCBI service, contact details and a five-second request delay.
# An API key is optional; leaving it as None is fine for this script.
# Note: NCBI allows up to 3 requests per second without a key, so 5 seconds is
# very cautious. It is a sensible choice on a shared university server, where
# many students may share one internet address. Lower it only if runs feel slow.
SLEEPBETWEEN = 5.0
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "s2924525@ed.ac.uk"
NCBI_API_KEY: Optional[str] = None

# NEW (improvement): NCBI asks every script to say which tool is calling it
# (the "tool" parameter) as well as giving an email. This lets NCBI contact you
# instead of blocking you if something goes wrong.
NCBI_TOOL = "bd_t3_gene_annotation"

# NEW (improvement): HTTP status codes that mean "temporary problem, try again"
# (429 = too many requests, 5xx = NCBI server trouble). Other errors such as
# 400 or 404 mean our request itself is wrong, so retrying would be pointless.
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


# Send a request, retry failures and return None if every attempt fails.
# The script pauses before EVERY request, including successful ones.
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
    # NEW (improvement): send the tool name (see NCBI_TOOL above).
    params.setdefault("tool", NCBI_TOOL)
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

            # NEW (improvement): stop at once for errors that retrying cannot fix.
            if resp.status_code not in RETRYABLE_STATUS:
                break

            # NEW (improvement): wait a little longer after each failed try
            # (5 s, then 10 s, ...) so a rate-limited request has time to clear.
            if attempt < max_retries:
                time.sleep(sleep_between * attempt)
        except requests.exceptions.RequestException as error:
            # Report the error type without exposing URLs or API keys.
            print(
                f"[ERROR] Attempt {attempt}/{max_retries} for {endpoint} "
                f"failed: {type(error).__name__}.",
                file=sys.stderr,
            )

    # NEW (fix): the old message always said "after 3 attempts", which is wrong
    # now that non-retryable errors stop early.
    print(f"[ERROR] Giving up on {endpoint}.", file=sys.stderr)
    return None


# NEW (improvement): ESummary can look up several Gene IDs in one request.
# This is used below to tell apart several search hits for one symbol.
# Retrieve structured annotation for Gene IDs, including names and chromosome.
# Returns a dictionary {gene_id: record}; IDs NCBI could not describe are left out.
def esummary_genes(gene_ids: List[str]) -> Dict[str, Dict[str, Any]]:
    params = {"db": "gene", "id": ",".join(gene_ids), "retmode": "json"}
    resp = ncbi_get("esummary.fcgi", params)
    if resp is None:
        return {}

    try:
        data = resp.json()
    except ValueError:
        print(f"[ERROR] Invalid ESummary JSON for '{','.join(gene_ids)}'.", file=sys.stderr)
        return {}

    result = data.get("result", {})
    usable: Dict[str, Dict[str, Any]] = {}
    for gene_id in gene_ids:
        record = result.get(gene_id)
        if isinstance(record, dict) and not record.get("error"):
            usable[gene_id] = record
    return usable


# Retrieve structured annotation for ONE Gene ID (a thin wrapper around the
# function above, so the rest of the script works exactly as before).
def esummary_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    record = esummary_genes([gene_id]).get(gene_id)
    if record is None:
        print(f"[WARN] No usable ESummary for '{gene_id}'.", file=sys.stderr)
    return record


# Find the NCBI Gene ID matching a gene symbol and organism.
# NEW (fix): the old version simply took the FIRST search hit. NCBI's
# [Gene Name] search also matches aliases, so for a symbol like TNF the first
# hit is not guaranteed to be the gene you meant. Now, when there are several
# hits, the script looks at each one's official symbol and picks the exact match.
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

    if len(id_list) == 1:
        return id_list[0]

    # NEW (fix): several hits -> choose the one whose official symbol matches.
    summaries = esummary_genes(id_list)
    for gene_id in id_list:
        official = str(summaries.get(gene_id, {}).get("name", ""))
        if official.upper() == symbol.upper():
            if gene_id != id_list[0]:
                print(
                    f"  -> {len(id_list)} matches for '{symbol}'; chose Gene ID "
                    f"{gene_id} because its official symbol matches exactly."
                )
            return gene_id

    print(
        f"[WARN] {len(id_list)} matches for '{symbol}' and none has that exact "
        f"official symbol; using Gene ID {id_list[0]}. Check it is the gene you want.",
        file=sys.stderr,
    )
    return id_list[0]


# Download the full NCBI Gene XML record so its database links can be inspected.
# NEW (fix): the old version sent rettype=gp. NCBI's table of valid EFetch
# options lists "gp" only for the protein database; for db=gene the XML record
# is requested with an empty rettype and retmode=xml. So rettype is now left out.
def efetch_gene_xml(gene_id: str) -> Optional[str]:
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "xml",
    }
    resp = ncbi_get("efetch.fcgi", params)
    if resp is None:
        return None
    return resp.text


# Collect unique ENSG identifiers from Ensembl database links in NCBI XML.
# ENSG identifiers are human gene IDs; other species use different prefixes.
def extract_ensembl_gene_ids_from_xml(xml_text: str) -> List[str]:
    ensg_ids: List[str] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as error:
        print(f"[ERROR] XML parse error: {error}", file=sys.stderr)
        return ensg_ids

    # Recognise Dbtag elements, including namespaced tags.
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


# NEW (improvement): tidy the gene list before using it. Removes spaces and
# empty entries, and removes repeats (so a gene typed twice, e.g. "TP53" and
# "tp53", is only looked up once). The original order is kept.
def clean_gene_list(gene_symbols: List[str]) -> List[str]:
    cleaned: List[str] = []
    seen = set()
    for symbol in gene_symbols:
        symbol = symbol.strip()
        if symbol and symbol.upper() not in seen:
            seen.add(symbol.upper())
            cleaned.append(symbol)
    return cleaned


# Process each gene and retain the same result keys for later project stages.
def fetch_ncbi_info_for_gene_list(
    gene_symbols: List[str],
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    results = []

    for symbol in clean_gene_list(gene_symbols):
        print(f"\n=== Processing gene symbol: {symbol} ===")
        gene_id = esearch_gene(symbol, organism=organism)
        if gene_id is None:
            continue

        print(f"  -> Resolved to NCBI Gene ID: {gene_id}")
        summary = esummary_gene(gene_id)
        if summary is None:
            continue

        # NEW (improvement): warn if NCBI's official symbol differs from what
        # you typed (for example you typed an old alias). Data is still kept.
        official = str(summary.get("name", ""))
        if official and official.upper() != symbol.upper():
            print(
                f"[WARN] You asked for '{symbol}' but NCBI's official symbol "
                f"is '{official}'.",
                file=sys.stderr,
            )

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


# Display the annotations in a readable form.
# Uses NCBI name/description when official nomenclature fields are empty.
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


# Write a TSV table for later database loading; multiple ENSG IDs use semicolons.
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

    # Remove embedded line breaks and delimiters from cells.
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

    # NEW (improvement): show the full path so you know where the file landed.
    print(f"[INFO] Wrote table to {os.path.abspath(filename)}")


# Save each complete ESummary record and its extracted Ensembl IDs as JSON.
def save_results_to_json(results: List[Dict[str, Any]], filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    # NEW (improvement): show the full path so you know where the file landed.
    print(f"[INFO] Saved JSON results to {os.path.abspath(filename)}")


# Fetch annotations, print them, save TSV/JSON files and return the records.
def connect_ncbi_gene_info(
    gene_list: List[str],
    output_filename: str,
    organism: str = "Homo sapiens",
) -> List[Dict[str, Any]]:
    # NEW (improvement): tidy the list first, then tell the user roughly how
    # long the run will take. Each gene needs 3 requests (search, summary, XML)
    # and the script waits SLEEPBETWEEN seconds before each one, so the quiet
    # pauses are normal and the script has not frozen.
    requested = clean_gene_list(gene_list)
    estimate = int(len(requested) * 3 * SLEEPBETWEEN)
    print(f"[INFO] Looking up {len(requested)} genes; expect about {estimate} seconds.")

    results = fetch_ncbi_info_for_gene_list(requested, organism=organism)
    print(f"[INFO] Retrieved {len(results)} of {len(requested)} requested gene entries.")

    # NEW (improvement): say exactly which genes failed, so a gene that is
    # missing from the output files is never silently overlooked.
    retrieved = {entry["query_symbol"] for entry in results}
    missing = [g for g in requested if g not in retrieved]
    if missing:
        print(f"[WARN] No annotation retrieved for: {', '.join(missing)}", file=sys.stderr)

    # NEW (improvement): also list genes that were found but have no Ensembl ID.
    no_ensg = [e["query_symbol"] for e in results if not e["ensembl_gene_ids"]]
    if no_ensg:
        print(f"[WARN] No Ensembl gene ID found for: {', '.join(no_ensg)}", file=sys.stderr)

    # NEW (improvement): if nothing at all came back, say so before writing
    # files that would otherwise contain only a header line.
    if not results:
        print("[WARN] Nothing retrieved; the output files will hold only headers.",
              file=sys.stderr)

    print_human_readable_summary(results)
    write_table_to_file(results, output_filename + ".tsv")
    save_results_to_json(results, output_filename + ".json")
    return results


# Edit gene_list for your assignment; output_filename is a prefix, without .tsv.
# Files are saved in the current working directory.
# The ten development genes were supplied by you.
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
    os.makedirs(output_dir, exist_ok=True)
    
    output_filename = os.path.join(output_dir, "ncbi_gene_summary2_try3")

    connect_ncbi_gene_info(
        gene_list,
        output_filename,
        organism=organism,
    )


# Run main() when this file is executed; importing it does not run the workflow.
if __name__ == "__main__":
    main()
