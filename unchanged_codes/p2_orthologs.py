#!/usr/bin/env python3
"""
combined_ncbi_orthologs_from_mixed_human_ids.py

Given a small example list of *mixed* human gene identifiers
(gene symbols and/or NCBI Gene IDs), this script:

1) Resolves each identifier to a unique human NCBI Gene ID and
   official symbol using NCBI ESearch + ESummary (db=gene).

2) Uses NCBI Datasets v2 Gene orthologs API to find mouse orthologs
   for those human NCBI Gene IDs.

3) Uses NCBI EFetch (db=gene, XML) to retrieve Ensembl Gene IDs
   (via dbXrefs) for both human and mouse NCBI Gene IDs.

This script was generated using ELM (elm.edina.ed.ac.uk) and edited.
Simon Tomlinson Nov 2026
"""

import sys
import csv
from typing import List, Dict, Any, Optional, Set

import requests
import xml.etree.ElementTree as ET

import time


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


# Mouse taxonomy ID (NCBI)
MOUSE_TAX_ID = 10090
SLEEPTIME=1

# NCBI endpoints
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_DATASETS_ORTHOLOGS_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/orthologs"
)


# ---------------------------------------------------------------------------
# Generic HTTP helper
# ---------------------------------------------------------------------------

def http_get(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    timeout: int = 30,
) -> Optional[requests.Response]:
    """
    Simple GET wrapper with basic error handling.
    """
    try:
        resp = requests.get(url, params=params, timeout=timeout)
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] GET {url} failed: {e}", file=sys.stderr)
        return None

    if resp.status_code != 200:
        print(
            f"[ERROR] GET {url} returned status {resp.status_code}: "
            f"{resp.text[:200]}",
            file=sys.stderr,
        )
        return None

    return resp


# ---------------------------------------------------------------------------
# Step 0: Resolve mixed human identifiers to NCBI Gene IDs + symbols
# ---------------------------------------------------------------------------

def resolve_human_identifier_to_gene(
    identifier: str,
) -> Optional[Dict[str, str]]:
    """
    Resolve a human gene identifier (symbol or NCBI Gene ID) to a
    unique NCBI Gene ID and official symbol using:

      1) A fast heuristic if identifier looks like a numeric Gene ID.
      2) ESearch (db=gene) restricted to human [9606].
      3) ESummary to obtain official symbol.

    Returns:
        {
          "input_id": <original identifier>,
          "gene_id":  <NCBI Gene ID as string>,
          "symbol":   <official gene symbol>
        }
    or None if resolution fails.
    """
    identifier = identifier.strip()
    if not identifier:
        return None

    # Heuristic: if it's purely digits, treat as candidate Gene ID.
    looks_numeric = identifier.isdigit()

    if looks_numeric:
        # Try ESummary directly on this numeric ID, restricted to human
        # by checking the taxid in the summary record.
        print(f"[INFO] Treating '{identifier}' as a candidate NCBI Gene ID")
        summary = fetch_gene_summary_by_id(identifier)
        if summary is not None and summary.get("taxid") == "9606":
            return {
                "input_id": identifier,
                "gene_id": summary["gene_id"],
                "symbol": summary["symbol"],
            }
        else:
            print(
                f"[WARN] '{identifier}' did not resolve as a human Gene ID; "
                f"falling back to symbol search.",
                file=sys.stderr,
            )

    # General path: ESearch for the identifier as a query term, restricted
    # to human [9606].
    return search_human_gene_by_term(identifier)


def search_human_gene_by_term(term: str) -> Optional[Dict[str, str]]:
    """
    Use ESearch (db=gene) to find a human gene record matching 'term',
    then ESummary to get the official symbol.

    Returns:
        {
          "input_id": <original term>,
          "gene_id":  <NCBI Gene ID as string>,
          "symbol":   <official gene symbol>
        }
    or None if resolution fails.
    """
    print(f"[INFO] Searching human gene for term: {term}")

    # ESearch restricted to human genes
    url = f"{NCBI_EUTILS_BASE}/esearch.fcgi"
    params = {
        "db": "gene",
        "term": f"{term}[Gene Name] AND 9606[TaxID]",
        "retmode": "json",
        "retmax": "1",
    }
    resp = http_get(url, params=params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError as e:
        print(f"[ERROR] Failed to parse ESearch JSON for term '{term}': {e}",
              file=sys.stderr)
        return None

    idlist = data.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        print(f"[WARN] No human Gene ID found for term '{term}'", file=sys.stderr)
        return None

    gene_id = idlist[0]
    summary = fetch_gene_summary_by_id(gene_id)
    if summary is None:
        return None

    return {
        "input_id": term,
        "gene_id": summary["gene_id"],
        "symbol": summary["symbol"],
    }


def fetch_gene_summary_by_id(gene_id: str) -> Optional[Dict[str, str]]:
    """
    Use ESummary (db=gene) to retrieve a small summary including
    taxid and official symbol for a single Gene ID.

    Returns:
        {
          "gene_id": <NCBI Gene ID>,
          "symbol":  <official gene symbol>,
          "taxid":   <taxon ID as string>,
        }
    or None on failure.
    """
    url = f"{NCBI_EUTILS_BASE}/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
    }
    resp = http_get(url, params=params)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError as e:
        print(f"[ERROR] Failed to parse ESummary JSON for Gene ID '{gene_id}': {e}",
              file=sys.stderr)
        return None

    result = data.get("result", {})
    doc = result.get(gene_id)
    if not doc:
        print(f"[WARN] No ESummary document for Gene ID '{gene_id}'", file=sys.stderr)
        return None

    taxid = str(doc.get("taxid", ""))
    symbol = (doc.get("name") or "").strip()
    gid = str(doc.get("uid", gene_id))

    if not symbol:
        print(
            f"[WARN] No official symbol found in ESummary for Gene ID '{gene_id}'",
            file=sys.stderr,
        )

    return {
        "gene_id": gid,
        "symbol": symbol,
        "taxid": taxid,
    }


# ---------------------------------------------------------------------------
# Step 1: Get mouse orthologs via NCBI Datasets API
# ---------------------------------------------------------------------------

def fetch_orthologs_for_gene(gene_id: str) -> Optional[Dict[str, Any]]:
    """
    Call NCBI Datasets gene orthologs web service for a given NCBI Gene ID.

    Endpoint:
      GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/orthologs

    Returns:
        Parsed JSON dict or None on failure.
    """
    url = NCBI_DATASETS_ORTHOLOGS_URL.format(gene_id=gene_id)
    resp = http_get(url)
    if resp is None:
        return None

    try:
        data = resp.json()
    except ValueError as e:
        print(f"[ERROR] Failed to parse JSON for gene_id {gene_id}: {e}",
              file=sys.stderr)
        return None

    return data


def extract_mouse_orthologs_from_json(
    orthologs_json: Dict[str, Any],
    mouse_tax_id: int = MOUSE_TAX_ID,
) -> List[Dict[str, Any]]:
    """
    From the orthologs JSON returned by NCBI Datasets, extract mouse orthologs.

    Expected structure (simplified / may evolve):

      {
        "reports": [
          {
            "gene": {
              "tax_id": 10090,
              "gene_id": "12189",
              "symbol": "Brca1",
              "name": "BRCA1, DNA repair associated",
              ...
            }
          },
          ...
        ]
      }

    Returns:
        List of dicts with keys:
          {
            "mouse_gene_id",
            "mouse_symbol",
            "mouse_name",
            "mouse_tax_id"
          }
    """
    mouse_orthos: List[Dict[str, Any]] = []

    orthologs = orthologs_json.get("reports", [])
    if not isinstance(orthologs, list):
        return mouse_orthos

    for o in orthologs:
        o2 = o.get("gene")
        if not isinstance(o2, dict):
            continue

        try:
            tax_id = int(o2.get("tax_id", -1))
        except (TypeError, ValueError):
            continue

        if tax_id != mouse_tax_id:
            continue

        mouse_orthos.append(
            {
                "mouse_gene_id": str(o2.get("gene_id", "")),
                "mouse_symbol": o2.get("symbol", "") or "",
                "mouse_name": o2.get("name", "") or "",
                "mouse_tax_id": tax_id,
            }
        )
        time.sleep(SLEEPTIME)

    return mouse_orthos


def build_human_mouse_ortholog_rows(
    resolved_genes: List[Dict[str, str]]
) -> List[Dict[str, str]]:
    """
    For each resolved human gene, call the Datasets orthologs service and
    collect mouse orthologs into row dicts.

    Each input dict must contain:
      "input_id", "gene_id", "symbol"

    Each output row dict has keys:
      "Human_Input_ID", "Human_Gene_ID", "Human_Symbol",
      "Mouse_Gene_ID",  "Mouse_Symbol",  "Mouse_Tax_ID"
    """
    rows: List[Dict[str, str]] = []

    for rec in resolved_genes:
        human_input = rec["input_id"]
        human_id = rec["gene_id"]
        human_symbol = rec["symbol"]

        print(
            f"\n=== Processing human gene {human_symbol} "
            f"(input '{human_input}', ID {human_id}) ==="
        )
        orthologs_json = fetch_orthologs_for_gene(human_id)
        if orthologs_json is None:
            print("  -> Failed to retrieve ortholog data.")
            continue

        mouse_orthos = extract_mouse_orthologs_from_json(orthologs_json)
        if not mouse_orthos:
            print("  -> No mouse orthologs found in web service response.")
            continue

        for m in mouse_orthos:
            row = {
                "Human_Input_ID": human_input,
                "Human_Gene_ID": human_id,
                "Human_Symbol": human_symbol,
                "Mouse_Gene_ID": m["mouse_gene_id"],
                "Mouse_Symbol": m["mouse_symbol"],
                "Mouse_Tax_ID": str(m["mouse_tax_id"]),
            }
            rows.append(row)
            time.sleep(SLEEPTIME)

    print(f"[INFO] Collected {len(rows)} human–mouse ortholog row(s).")
    return rows


# ---------------------------------------------------------------------------
# Step 2: EFetch NCBI Gene XML and extract Ensembl IDs
# ---------------------------------------------------------------------------

def efetch_gene_xml(gene_ids: List[str]) -> Dict[str, str]:
    """
    Use NCBI EFetch to retrieve XML records for a list of Gene IDs.

    Args:
        gene_ids: list of NCBI Gene IDs as strings.

    Returns:
        Dict mapping gene_id -> XML text (single <Entrezgene> snippet).
    """
    xml_by_id: Dict[str, str] = {}
    if not gene_ids:
        return xml_by_id

    url = f"{NCBI_EUTILS_BASE}/efetch.fcgi"
    params = {
        "db": "gene",
        "id": ",".join(gene_ids),
        "retmode": "xml",
    }
    resp = http_get(url, params=params)
    if resp is None:
        return xml_by_id

    text = resp.text

    try:
        root = ET.fromstring(text)
    except ET.ParseError as e:
        print(f"[ERROR] Failed to parse combined Gene XML: {e}", file=sys.stderr)
        return xml_by_id

    # Typical structure:
    # <Entrezgene-Set>
    #   <Entrezgene> ... </Entrezgene>
    #   <Entrezgene> ... </Entrezgene>
    # </Entrezgene-Set>
    for entrezgene in root.findall(".//Entrezgene"):
        gene_id = None
        track_info = entrezgene.find(".//Gene-track_geneid")
        if track_info is not None and track_info.text:
            gene_id = track_info.text.strip()

        if gene_id:
            xml_snippet = ET.tostring(entrezgene, encoding="unicode")
            xml_by_id[gene_id] = xml_snippet

    return xml_by_id


def extract_ensembl_ids_from_gene_xml(xml_text: str) -> List[str]:
    """
    Parse a single Gene XML entry (<Entrezgene>) and extract Ensembl gene IDs
    from dbXrefs (Dbtag).

    Strategy:
    - Look for Dbtag where Dbtag_db contains "Ensembl".
    - Collect text values that look like Ensembl gene IDs:
        - Start with "ENSG" (human) or "ENSMUSG" (mouse).
    """
    ensembl_ids: List[str] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"[ERROR] XML parse error in extract_ensembl_ids_from_gene_xml: {e}",
              file=sys.stderr)
        return ensembl_ids

    for dbtag in root.iter():
        if not dbtag.tag.endswith("Dbtag"):
            continue

        db_name = None
        for child in dbtag:
            if child.tag.endswith("Dbtag_db"):
                db_name = (child.text or "").strip()
        if not db_name:
            continue

        if "ensembl" not in db_name.lower():
            continue

        # Within Dbtag_tag, collect all text snippets and pick those starting
        # with ENSG or ENSMUSG (gene-level Ensembl IDs).
        for child in dbtag:
            if not child.tag.endswith("Dbtag_tag"):
                continue

            texts: List[str] = []
            for sub in child.iter():
                if sub.text:
                    texts.append(sub.text.strip())

            for t in texts:
                if t.startswith("ENSG") or t.startswith("ENSMUSG"):
                    if t not in ensembl_ids:
                        ensembl_ids.append(t)

    return ensembl_ids


def add_ensembl_ids_to_rows(
    rows: List[Dict[str, str]]
) -> List[Dict[str, str]]:
    """
    Given human–mouse ortholog rows with:
      Human_Input_ID, Human_Gene_ID, Human_Symbol,
      Mouse_Gene_ID, Mouse_Symbol, Mouse_Tax_ID

    Use EFetch to retrieve Ensembl gene IDs for all unique
    human and mouse gene IDs, and return an updated list of rows
    with added fields:
      Human_Ensembl_Gene_ID, Mouse_Ensembl_Gene_ID
    """
    if not rows:
        return []

    # Collect unique human and mouse Gene IDs
    human_ids: Set[str] = set()
    mouse_ids: Set[str] = set()

    for row in rows:
        h_id = (row.get("Human_Gene_ID") or "").strip()
        m_id = (row.get("Mouse_Gene_ID") or "").strip()
        if h_id:
            human_ids.add(h_id)
        if m_id:
            mouse_ids.add(m_id)

    print(f"[INFO] Unique human Gene IDs: {len(human_ids)}")
    print(f"[INFO] Unique mouse Gene IDs: {len(mouse_ids)}")

    # Fetch XML for human and mouse Gene IDs
    human_xml_by_id = efetch_gene_xml(sorted(human_ids))
    mouse_xml_by_id = efetch_gene_xml(sorted(mouse_ids))

    # Extract Ensembl gene IDs from XML and build lookup dicts
    human_ens_by_id: Dict[str, str] = {}
    for gid, xml_text in human_xml_by_id.items():
        ids = extract_ensembl_ids_from_gene_xml(xml_text)
        ens_h = [i for i in ids if i.startswith("ENSG")]
        if ens_h:
            human_ens_by_id[gid] = ";".join(ens_h)
        elif ids:
            human_ens_by_id[gid] = ";".join(ids)

    mouse_ens_by_id: Dict[str, str] = {}
    for gid, xml_text in mouse_xml_by_id.items():
        ids = extract_ensembl_ids_from_gene_xml(xml_text)
        ens_m = [i for i in ids if i.startswith("ENSMUSG")]
        if ens_m:
            mouse_ens_by_id[gid] = ";".join(ens_m)
        elif ids:
            mouse_ens_by_id[gid] = ";".join(ids)

    # Attach Ensembl IDs to each row
    updated_rows: List[Dict[str, str]] = []
    for row in rows:
        h_id = (row.get("Human_Gene_ID") or "").strip()
        m_id = (row.get("Mouse_Gene_ID") or "").strip()

        new_row = dict(row)  # shallow copy
        new_row["Human_Ensembl_Gene_ID"] = human_ens_by_id.get(h_id, "")
        new_row["Mouse_Ensembl_Gene_ID"] = mouse_ens_by_id.get(m_id, "")
        updated_rows.append(new_row)

    return updated_rows


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def connect_ncbi_ortholog(EXAMPLE_HUMAN_INPUTS, ncbi_ortholog_file):
    # Step 0: Resolve mixed human identifiers to NCBI Gene IDs and symbols
    resolved_genes: List[Dict[str, str]] = []
    for ident in EXAMPLE_HUMAN_INPUTS:
        res = resolve_human_identifier_to_gene(ident)
        if res is None:
            print(
                f"[WARN] Could not resolve identifier '{ident}' to a human Gene ID",
                file=sys.stderr,
            )
        else:
            print(
                f"[INFO] Resolved '{ident}' -> Gene ID {res['gene_id']} "
                f"({res['symbol']})"
            )
            resolved_genes.append(res)
        time.sleep(SLEEPTIME)
    if not resolved_genes:
        print("[ERROR] No human genes resolved; exiting.", file=sys.stderr)
        sys.exit(1)

    # Step 1: Use Datasets API to build ortholog rows
    ortholog_rows = build_human_mouse_ortholog_rows(resolved_genes)
    if not ortholog_rows:
        print("[ERROR] No ortholog rows generated; exiting.", file=sys.stderr)
        sys.exit(1)

    # Step 2: Use EFetch to attach Ensembl IDs
    rows_with_ensembl = add_ensembl_ids_to_rows(ortholog_rows)

    # Final output file
    output_file = ncbi_ortholog_file +".tsv"

    header = [
        "Human_Input_ID",
        "Human_Gene_ID",
        "Human_Symbol",
        "Human_Ensembl_Gene_ID",
        "Mouse_Gene_ID",
        "Mouse_Symbol",
        "Mouse_Tax_ID",
        "Mouse_Ensembl_Gene_ID",
    ]

    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)

        for row in rows_with_ensembl:
            writer.writerow(
                [
                    row.get("Human_Input_ID", ""),
                    row.get("Human_Gene_ID", ""),
                    row.get("Human_Symbol", ""),
                    row.get("Human_Ensembl_Gene_ID", ""),
                    row.get("Mouse_Gene_ID", ""),
                    row.get("Mouse_Symbol", ""),
                    row.get("Mouse_Tax_ID", ""),
                    row.get("Mouse_Ensembl_Gene_ID", ""),
                ]
            )

    print(f"[INFO] Wrote TAB-delimited table to {output_file}")
    return(rows_with_ensembl)

