# =============================================================================
# GenAI declaration: Tool(s): Claude (Anthropic)
# Use: feedback, review and comments on code I adapted from simon tomlinson
# =============================================================================
"""Part 2: map each Part 1 human gene to its mouse ortholog(s) in NCBI.

Adapted from Simon Tomlinson's ELM-assisted OrthologMap3_SRT.py (2026).

Input:  OP1_data/part1_data/final_data/<dataset>/ncbi_gene_summary2.tsv (Part 1)
Output: OP1_data/part2_data/<dataset>/ncbi_ortholog.tsv

Every Part 1 input row produces at least one output row, in input order and
numbered by Input_Row, so unresolved inputs are kept. A human gene with
several mouse orthologs produces one row per mouse gene; the human-side
columns are repeated on each of those rows. Mapping_Status records the outcome:

  ortholog_found / multiple_mouse_orthologs   one or more mouse genes found
  no_mouse_ortholog_returned                  lookup worked, NCBI lists none
  input_not_resolved                          Part 1 did not resolve the input
  invalid_human_gene_id, input_gene_is_not_human, human_record_not_current,
  human_lookup_failed, ortholog_lookup_failed problems that need attention

Shared helpers (HTTP, NCBI lookups, file I/O, the run() entry point and its
exit codes) come from the Part 1 module p1_basic_anno.py, which must be in
the same folder.
"""

import re
from functools import cache  # Unused at present (see review note on mouse_genes).

from p1_basic_anno import (API_ERRORS, NCBI_GENES, ServiceError, current_human, gene_ensembl_ids, gene_summaries, get_json, part_dir, read_table, run, write_tsv)

# Output column order; the MySQL loader requires exactly this header.
COLUMNS = """Human_Input_ID Human_Gene_ID Human_Symbol Human_Ensembl_Gene_ID Mouse_Gene_ID
             Mouse_Symbol Mouse_Tax_ID Mouse_Ensembl_Gene_ID Mouse_Official_Name Input_Row
             Same_Gene_Note Resolution_Status Resolution_Note Mouse_Ortholog_Count Mapping_Status
             Mapping_Note Human_Ensembl_Status Mouse_Ensembl_Status""".split()
# Mapping statuses that are valid results rather than lookup problems.
MAPPING_OUTCOMES = {"ortholog_found", "multiple_mouse_orthologs", "no_mouse_ortholog_returned", "input_not_resolved"}

def mouse_genes(gid):
    """Return the distinct mouse (taxon 10090) ortholog gene records for a human Gene ID.

    Queries the NCBI Datasets v2 orthologs endpoint, following
    next_page_token until all pages are read. Each record is the Datasets
    "gene" object (gene_id, symbol, description, ensembl_gene_ids, ...).
    Raises ServiceError if a page is malformed, a report is flagged as
    incomplete, a mouse Gene ID is invalid or conflicts with an earlier one,
    or the API repeats a page token (to avoid an endless loop).
    """
    url = f"{NCBI_GENES}{gid}/orthologs"
    params = {"taxon_filter": "10090", "page_size": 100}
    mice, seen = {}, set()
    while True:
        page = get_json(url, params)
        reports = page.get("reports")
        # Datasets omits "reports" (or returns {}) when there are no orthologs.
        if reports is None and (not page or str(page.get("total_count")) == "0"):
            reports = []
        if not isinstance(reports, list):
            raise ServiceError("NCBI did not return an ortholog report list.")
        for item in reports:
            gene = item["gene"]
            if item.get("errors"):
                raise ServiceError("NCBI marked an ortholog report as incomplete.")
            # Keep mouse genes only; other taxa in the report are ignored.
            if str(gene["tax_id"]) == "10090":
                mouse_id = str(gene["gene_id"])
                if not re.fullmatch(r"[1-9]\d*", mouse_id) or (mouse_id in mice and mice[mouse_id] != gene):
                    raise ServiceError("Invalid or conflicting mouse Gene ID.")
                mice[mouse_id] = gene  # Keyed by Gene ID, so repeats across pages collapse.
        token = page.get("next_page_token")
        if not token:
            # Last page: flag any difference between NCBI's count and the
            # mouse genes kept.
            if "total_count" in page and int(page["total_count"]) != len(mice):
                print(f"[WARN] Gene {gid}: NCBI counted {page['total_count']} reports. Kept {len(mice)} distinct mouse genes.")
            break
        if token in seen:
            raise ServiceError("NCBI repeated an ortholog page.")
        seen.add(token)
        params["page_token"] = token
    return list(mice.values())


def ensembl_field(gid, values, prefix, source):
    """Return (semicolon-joined Ensembl IDs, status) for one gene.

    Uses the IDs already supplied in `values` when any are valid for this
    species (`prefix`, e.g. "ENSG" or "ENSMUSG"), and reports `source` as
    the status. Otherwise falls back to the NCBI Gene XML cross-references.
    The status is then "from_ncbi_xml", "no_id_returned" or "lookup_failed".
    Version suffixes are removed.
    """
    ids = {v.split(".")[0] for v in values if isinstance(v, str) and re.fullmatch(prefix + r"\d+(\.\d+)?", v)}
    if ids:
        return ";".join(sorted(ids)), source
    try:
        found = gene_ensembl_ids(gid, prefix)
        return ";".join(found), "from_ncbi_xml" if found else "no_id_returned"
    except API_ERRORS:
        return "", "lookup_failed"


def map_input(number, source, known):
    """Return the output row(s) for one Part 1 input row.

    `number` is the 1-based Input_Row, `source` the Part 1 row and `known` the
    per-dataset ESummary cache shared with gene_summaries(). The human Gene ID
    is checked again (live, human, current) before the ortholog lookup, so a
    Part 1 file that has gone stale is not mapped silently. Every return path
    gives a Mapping_Status; a failure is recorded in the row, not raised.
    """
    gid = source["Gene_ID"]
    # Carry the Part 1 resolution columns through so each row explains itself.
    row = {key: source.get(key, "") for key in ("Resolution_Status", "Resolution_Note", "Same_Gene_Note")}
    row.update(Input_Row=number, Human_Input_ID=source["Query_Symbol"], Human_Gene_ID=gid,
               Human_Symbol=source["Official_Symbol"], Human_Ensembl_Gene_ID=source["ENSG_IDs"],
               Human_Ensembl_Status="from_part1" if source["ENSG_IDs"] else "not_requested",
               Mouse_Ensembl_Status="not_requested")
    if row["Resolution_Status"] != "resolved":
        return [dict(row, Mapping_Status="input_not_resolved", Mapping_Note="See Part 1 resolution.")]
    if not re.fullmatch(r"[1-9]\d*", gid):
        return [dict(row, Mapping_Status="invalid_human_gene_id", Mapping_Note="Invalid Part 1 Gene ID.")]
    # Re-validate the human record before asking for its orthologs.
    try:
        record = gene_summaries([gid], known)[0]
        taxon = (record.get("organism") or {}).get("taxid") or record.get("taxid")
        if not taxon:
            raise ServiceError("Missing taxonomy ID in the NCBI summary.")
        if str(taxon) != "9606":
            return [dict(row, Mapping_Status="input_gene_is_not_human", Mapping_Note="NCBI returned another species.")]
        if not current_human(record):
            return [dict(row, Mapping_Status="human_record_not_current", Mapping_Note="NCBI record is no longer current.")]
    except API_ERRORS as error:
        return [dict(row, Mapping_Status="human_lookup_failed", Mapping_Note=str(error))]
    try:
        mice = mouse_genes(gid)
    except API_ERRORS as error:
        return [dict(row, Mapping_Status="ortholog_lookup_failed", Mapping_Note=str(error))]
    # Human Ensembl IDs: reuse Part 1's, or look them up if Part 1 had none.
    row["Human_Ensembl_Gene_ID"], row["Human_Ensembl_Status"] = ensembl_field(
        gid, source["ENSG_IDs"].split(";"), "ENSG", "from_part1")
    row["Mouse_Ortholog_Count"] = len(mice)  # Left blank on failure rows: count not assessed.
    if not mice:
        return [dict(row, Mapping_Status="no_mouse_ortholog_returned", Mapping_Note="NCBI returned no mouse ortholog.")]
    row["Mapping_Status"] = "ortholog_found" if len(mice) == 1 else "multiple_mouse_orthologs"
    # One output row per mouse ortholog; the human columns are repeated.
    rows = []
    for gene in mice:
        mouse_id = str(gene["gene_id"])
        ids, status = ensembl_field(mouse_id, gene.get("ensembl_gene_ids") or [], "ENSMUSG", "from_ncbi_datasets")
        rows.append(dict(row, Mouse_Gene_ID=mouse_id, Mouse_Tax_ID="10090", Mouse_Symbol=gene.get("symbol"),
                         Mouse_Official_Name=gene.get("description"), Mouse_Ensembl_Gene_ID=ids,
                         Mouse_Ensembl_Status=status))
    return rows


def process(data_root, dataset):
    """Map one dataset's Part 1 output to mouse orthologs and write the TSV.

    Returns True (the run exits with code 2) if any row has a lookup problem
    or carries a failed Part 1 lookup, and, for the development list, if any
    input is unresolved, since every development gene is expected to resolve.
    """
    # plain=True: Part 1 writes unquoted TSV, so quotes are literal text.
    inputs = read_table(part_dir(data_root, 1, dataset) / "ncbi_gene_summary2.tsv", plain=True)
    known = {}
    rows = [row for number, source in enumerate(inputs, 1) for row in map_input(number, source, known)]
    # Default (quoted) TSV: Part 3 and the MySQL loader read Part 2 with CSV quoting.
    write_tsv(part_dir(data_root, 2, dataset) / "ncbi_ortholog.tsv", COLUMNS, rows)
    print(f"Part 2 {dataset}: {len(inputs)} inputs, {len(rows)} mapping rows.")
    failed = any(r["Mapping_Status"] not in MAPPING_OUTCOMES or r["Resolution_Status"] == "failed" or
                 "lookup_failed" in (r["Human_Ensembl_Status"], r["Mouse_Ensembl_Status"]) for r in rows)
    return failed or (dataset == "development" and any(r["Resolution_Status"] != "resolved" for r in rows))


if __name__ == "__main__":
    run(process)
