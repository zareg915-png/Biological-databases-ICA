#!/usr/bin/env python3
"""Part 3: compare each human gene with its mouse ortholog using UniProt and Bgee.

Input:  OP1_data/part2_data/<dataset>/ncbi_ortholog.tsv (Part 2)
Output: OP1_data/part3_data/<dataset>/
  p3_pair_summary.tsv  one row per input and mouse ortholog: identifiers, each
                       species' UniProt and Bgee fields (accessions, protein
                       names, GO terms, function, locations), Jaccard scores
                       and the expression comparison
  p3_results.json      {"genes": one record per distinct gene looked up,
                        "proteins": one record per UniProt entry used};
                       loaded by all_sql_loads.py as p3_gene_status and p3_uniprot

For every distinct resolved human and mouse gene (looked up once, even if
it appears in several rows) the script:
  1. confirms the gene's identity, type and Ensembl IDs with NCBI Datasets;
  2. finds its UniProtKB entries, verified against the NCBI Gene ID or
     Ensembl IDs, preferring reviewed (Swiss-Prot) entries;
  3. collects protein names, GO terms (MF, BP, CC), UniProt subcellular
     locations, function text and location evidence codes;
  4. retrieves Bgee present/absent expression calls per anatomical entity.

Each Part 2 row (one per input and mouse ortholog) then gets one summary row
with Jaccard similarities for the GO and location sets and a tissue-level
expression comparison. Rows that cannot be compared keep a status and a note
(Pair_Note) explaining why.

Shared helpers and the run() entry point come from p1_basic_anno.py, and
MAPPING_OUTCOMES from p2_orthologs.py; both files must be in the same folder.
"""
import json
import re
import time
from collections import defaultdict
import importlib.util, sys
from pathlib import Path
CODE = Path(__file__).resolve().parents[1]   # the data_ret_script folder

def _load(name, relative_path):
    """Import another part's script under the module name this code expects."""
    spec = importlib.util.spec_from_file_location(name, CODE / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)

_load("p1_basic_anno", "part1/final_script_run_for_both_lists/final_script_for_both_end_res.py")
_load("p2_orthologs", "part2/final_ortho_script.py")
from p1_basic_anno import (API_ERRORS, NCBI_GENES, ServiceError, decode, fetch, get_json, part_dir, read_table, run, write_tsv)
from p2_orthologs import MAPPING_OUTCOMES

UNIPROT = "https://rest.uniprot.org"
BGEE = "https://www.bgee.org/api/"
TAXA = {"Human": "9606", "Mouse": "10090"}
# NCBI gene types (normalised to upper case with underscores) that are not expected to encode a protein, so no UniProt entry is sought for them.
NONCODING = set("NCRNA LNCRNA RRNA TRNA SNRNA SNORNA SCRNA MIRNA PSEUDO PSEUDOGENE NONCODING BIOLOGICAL_REGION".split())
# Annotation sets compared between species with the Jaccard index.
SETS = ("GO_MF", "GO_BP", "GO_CC", "Location")
# Per-species columns, written once with a Human_ prefix and once with Mouse_.
SIDE_FIELDS = "UniProt_Status UniProt_Accessions Protein_Names GO_Terms Function Subcellular_Locations Location_Evidence Bgee_Status".split()
# Output column order for p3_pair_summary.tsv.
COLUMNS = """Dataset Input_Row Query_Symbol Human_Gene_ID Human_Symbol Mouse_Gene_ID Mouse_Symbol
             Resolution_Status Mapping_Status""".split()
COLUMNS += [f"{side}_{field}" for side in TAXA for field in SIDE_FIELDS]
COLUMNS += [kind + "_Jaccard" for kind in SETS]
COLUMNS += """Expression_Comparison_Status Comparable_Tissue_Count Both_Present_Count Both_Absent_Count
             Shared_Present_Tissues Expression_Call_Agreement Present_Tissue_Jaccard Pair_Note""".split()


# UniProt
def uniprot_pages(url, params):
    """Yield every result from a paginated UniProt REST query.

    UniProt gives the next page as a URL in the HTTP Link header; that URL
    already contains the query, so `params` are only sent with the first
    request. A page URL seen twice raises ServiceError to prevent a loop.
    """
    seen = set()
    while url:
        if url in seen:
            raise ServiceError("UniProt repeated a results page.")
        seen.add(url)
        body, headers = fetch(url, params)
        yield from decode(body)["results"]
        link = re.search(r'<([^>]+)>;\s*rel="?next"?', headers.get("Link", ""))
        url, params = (link[1] if link else None), None


def mapped_entries(ids, database):
    """Yield UniProtKB entries mapped from Ensembl IDs via the UniProt ID-mapping service.

    `database` is the UniProt source database name ("Ensembl" for gene IDs or
    "Ensembl_Transcript"). The job is submitted, polled every 2 s for up to
    about a minute, and its results are paged. Only rows whose source ID
    (ignoring any version suffix) was one of the requested IDs are kept.
    """
    base = UNIPROT + "/idmapping"
    form = {"from": database, "to": "UniProtKB", "ids": ",".join(ids)}
    job = get_json(base + "/run", form=form)["jobId"]  # POST; not retried by fetch().
    for _ in range(30):
        # A finished job may redirect to its results, which carry no jobStatus.
        state = get_json(f"{base}/status/{job}").get("jobStatus")
        if state in (None, "FINISHED"):
            break
        if state not in ("RUNNING", "NEW"):
            raise ServiceError(f"UniProt mapping job returned {state}.")
        time.sleep(2)
    else:
        raise ServiceError("UniProt mapping job timed out.")
    for row in uniprot_pages(f"{base}/uniprotkb/results/{job}", {"format": "json", "size": 500}):
        if row["from"].split(".")[0] in ids:
            yield row["to"]


def transcripts(gene):
    """Return the unversioned Ensembl transcript IDs of a gene's Ensembl gene IDs.

    Used as a last-resort route to UniProt. Ensembl must confirm both the
    gene ID and the expected species; only transcripts whose parent is that
    gene are kept.
    """
    species = "homo_sapiens" if gene["Tax_ID"] == "9606" else "mus_musculus"
    ids = set()
    for gid in gene["Ensembl_IDs"]:
        record = get_json("https://rest.ensembl.org/lookup/id/" + gid, {"expand": 1, "content-type": "application/json"})
        if record["id"] != gid or record["species"] != species:
            raise ServiceError("Ensembl returned a different gene or species.")
        ids.update(t["id"].split(".")[0] for t in record.get("Transcript", []) if t.get("Parent") == gid)
    return ids


def verified(entries, gene, mapped=False):
    """Keep only UniProt entries that belong to this gene, without duplicates.

    An entry must be from the gene's species and must not point to a
    different NCBI Gene ID. Entries from a Gene ID search (mapped=False) must
    name this Gene ID. Entries reached through Ensembl ID mapping
    (mapped=True) are also rejected if they list only other Ensembl genes.
    """
    result = {}
    for entry in entries:
        if str(entry["organism"]["taxonId"]) != gene["Tax_ID"]:
            continue
        xrefs = entry.get("uniProtKBCrossReferences", [])
        ids = {x["id"] for x in xrefs if x["database"] == "GeneID"}
        # Ensembl cross-references are transcripts; the gene ID is a property.
        genes = {p["value"].split(".")[0] for x in xrefs if x["database"] == "Ensembl"
                 for p in x.get("properties", []) if p["key"] == "GeneId"}
        if ids and gene["Gene_ID"] not in ids:
            continue
        if mapped and genes and genes.isdisjoint(gene["Ensembl_IDs"]):
            continue
        if mapped or gene["Gene_ID"] in ids:
            result[entry["primaryAccession"]] = entry
    return list(result.values())


def protein_data(gene):
    """Return (UniProt status, entries) for one gene.

    Search order, stopping at the first route that finds verified entries:
    NCBI Gene ID cross-reference search, then Ensembl gene ID mapping, then
    Ensembl transcript ID mapping. Reviewed (Swiss-Prot) entries are used
    when any exist; otherwise unreviewed (TrEMBL) entries are used and the
    status says so. Non-coding genes are skipped.
    """
    if gene["Gene_Type"].upper().replace("-", "_").replace(" ", "_") in NONCODING:
        return "not_applicable_gene_type", []
    query = f"organism_id:{gene['Tax_ID']} AND xref:geneid-{gene['Gene_ID']}"
    params = {"query": query, "format": "json", "size": 500}
    entries = verified(uniprot_pages(UNIPROT + "/uniprotkb/search", params), gene)
    if not entries and gene["Ensembl_IDs"]:
        entries = verified(mapped_entries(gene["Ensembl_IDs"], "Ensembl"), gene, mapped=True)
        if not entries:
            ids = transcripts(gene)
            if ids:
                entries = verified(mapped_entries(sorted(ids), "Ensembl_Transcript"), gene, mapped=True)
    reviewed = [e for e in entries if e.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"]
    if reviewed:
        return "reviewed_found", reviewed
    return ("unreviewed_only" if entries else "no_verified_entry"), entries


def protein_name(entry):
    """Return the entry's protein name.

    Reviewed entries have a recommended name; unreviewed (TrEMBL) entries
    often have only a submitted name. Returns "" if neither is present.
    """
    description = entry.get("proteinDescription", {})
    name = description.get("recommendedName") or (description.get("submissionNames") or [{}])[0]
    return name.get("fullName", {}).get("value", "")


def protein_annotations(entries):
    """Summarise a gene's UniProt entries for the output table.

    Returns the accessions, protein names, GO terms ("GO:ID name"), function
    text, subcellular location names and location evidence (ECO codes) as
    text, plus `sets`: the GO MF/BP/CC term IDs and UniProt subcellular
    location IDs used for the Jaccard comparison. All evidence levels are
    included.
    """
    sets = {kind: set() for kind in SETS}
    functions, locations, evidence, go_terms = set(), set(), set(), set()
    for entry in entries:
        for c in entry.get("comments", []):
            if c["commentType"] == "FUNCTION":
                functions.update(t["value"] for t in c.get("texts", []) if t.get("value"))
            if c["commentType"] == "SUBCELLULAR LOCATION":
                for item in c.get("subcellularLocations", []):
                    loc = item.get("location") or {}
                    if loc.get("id"):
                        sets["Location"].add(loc["id"])  # Controlled ID, e.g. SL-0191.
                    if loc.get("value"):
                        locations.add(loc["value"])
                        evidence.update(f"{loc['value']} ({e['evidenceCode']})" for e in loc.get("evidences", [])
                                        if e.get("evidenceCode"))
        for ref in entry.get("uniProtKBCrossReferences", []):
            if ref["database"] == "GO":
                # GoTerm is "F:...", "P:..." or "C:...": the letter gives the GO aspect.
                props = {p["key"]: p.get("value", "") for p in ref.get("properties", [])}
                kind = {"F": "GO_MF", "P": "GO_BP", "C": "GO_CC"}.get(props.get("GoTerm", "")[:1])
                if kind:
                    sets[kind].add(ref["id"])
                    go_terms.add(f'{ref["id"]} {props["GoTerm"][2:]}')  # e.g. "GO:0003677 DNA binding".
    return dict(UniProt_Accessions="; ".join(sorted(e["primaryAccession"] for e in entries)),
                Protein_Names="; ".join(sorted({protein_name(e) for e in entries} - {""})),
                GO_Terms="; ".join(sorted(go_terms)),
                Function=" | ".join(sorted(functions)), Subcellular_Locations="; ".join(sorted(locations)),
                Location_Evidence="; ".join(sorted(evidence)), sets=sets)


# Bgee
def expression_data(gene):
    """Return (Bgee status, calls) for one gene.

    `calls` maps an UBERON anatomical entity ID to ("present" or "absent",
    tissue name). Only gold- and silver-quality calls on UBERON anatomy are
    kept. A tissue with more than one call (for example both present and
    absent) is treated as conflicting and dropped. Bgee is only queried when
    the gene has exactly one Ensembl ID, so the calls clearly belong to it.
    A 404 means Bgee does not have the gene.
    """
    ids = gene["Ensembl_IDs"]
    if len(ids) != 1:
        return ("missing_ensembl_id" if not ids else "multiple_ensembl_ids"), {}
    calls = {}
    total = 0    # All calls returned, before filtering.
    missing = 0  # Requests answered with 404.
    for requested in ("EXPRESSED", "NOT_EXPRESSED"):
        params = dict(page="gene", action="expression", gene_id=ids[0], species_id=gene["Tax_ID"],
                      display_type="json", cond_param="anat_entity", expr_type=requested)
        try:
            records = get_json(BGEE, params)["data"]["calls"]
        except ServiceError as error:
            if error.status != 404:
                raise
            missing += 1
            continue
        if not isinstance(records, list):
            raise ServiceError("Bgee returned an invalid call list.")
        total += len(records)
        for record in records:
            tissue = record["condition"]["anatEntity"]
            tissue_id = tissue["id"]
            if tissue_id.startswith("UBERON:") and (record.get("expressionQuality") or "").lower() in ("gold", "silver"):
                state = {"expressed": "present", "not_expressed": "absent"}[record["expressionState"]]
                # Second call for the same tissue: mark it as conflicting (None).
                calls[tissue_id] = None if tissue_id in calls else (state, tissue.get("name") or tissue_id)
    if missing:
        return ("not_found" if missing == 2 else "incomplete"), {}
    return ("retrieved" if total else "no_calls_returned"), {t: call for t, call in calls.items() if call}


# Per-gene annotation and pair comparison
def annotate(gid, side, ensembl):
    """Collect NCBI, UniProt and Bgee data for one gene ("Human" or "Mouse" side).

    `ensembl` holds the Ensembl IDs carried from Parts 1-2; NCBI's current IDs
    replace them when available. Each source is tried independently, so a
    failure in one is recorded in `notes` and a "lookup_failed" status
    without stopping the others. If NCBI reports a different gene or species,
    nothing else is looked up.
    """
    taxon = TAXA[side]
    gene = dict(Gene_ID=gid, Tax_ID=taxon, Gene_Type="", Ensembl_IDs=ensembl, notes=[])
    try:
        rec = get_json(f"{NCBI_GENES}{gid}/dataset_report")["reports"][0]["gene"]
        if str(rec["gene_id"]) != gid or str(rec["tax_id"]) != taxon:
            gene.update(UniProt_Status="invalid_gene_identity", Bgee_Status="invalid_gene_identity")
            gene["notes"].append("NCBI returned a different gene or species.")
            return gene
        gene["Gene_Type"] = rec.get("type") or ""
        prefix = "ENSG" if side == "Human" else "ENSMUSG"
        current_ids = {v.split(".")[0] for v in rec.get("ensembl_gene_ids", [])
                       if re.fullmatch(prefix + r"\d+(\.\d+)?", v)}
        gene["Ensembl_IDs"] = sorted(current_ids) or ensembl
    except API_ERRORS as error:
        # Carry on with the Part 2 Ensembl IDs; gene type stays unknown.
        gene["notes"].append(f"NCBI metadata: {error}")
    try:
        gene["UniProt_Status"], entries = protein_data(gene)
        gene["entries"] = entries  # Kept for the per-protein JSON table.
        gene.update(protein_annotations(entries))
    except API_ERRORS as error:
        gene["UniProt_Status"] = "lookup_failed"
        gene["notes"].append(f"UniProt: {error}")
    try:
        gene["Bgee_Status"], gene["calls"] = expression_data(gene)
    except API_ERRORS as error:
        gene["Bgee_Status"] = "lookup_failed"
        gene["notes"].append(f"Bgee: {error}")
    return gene


def compare_expression(human, mouse):
    """Compare Bgee calls over the tissues that have a usable call in both species.

    Expression_Call_Agreement is the share of those tissues with the same call
    (both present or both absent). Present_Tissue_Jaccard is the number of
    tissues present in both divided by the number present in at least one;
    it is left empty when no tissue is present in either species.
    """
    finished = {"retrieved", "no_calls_returned"}
    if any(g.get("Bgee_Status") not in finished for g in (human, mouse)):
        return {"Expression_Comparison_Status": "not_comparable"}
    left, right = human["calls"], mouse["calls"]
    both = left.keys() & right.keys()
    if not both:
        return dict(Expression_Comparison_Status="no_shared_comparable_tissues", Comparable_Tissue_Count=0)
    present = sorted(t for t in both if left[t][0] == right[t][0] == "present")
    absent = sum(left[t][0] == right[t][0] == "absent" for t in both)
    any_present = len(both) - absent
    return dict(Expression_Comparison_Status="comparable", Comparable_Tissue_Count=len(both),
                Both_Present_Count=len(present), Both_Absent_Count=absent,
                Shared_Present_Tissues="; ".join(f"{left[t][1]} ({t})" for t in present),
                Expression_Call_Agreement=(len(present) + absent) / len(both),
                Present_Tissue_Jaccard=len(present) / any_present if any_present else None)


def pair_summary(dataset, row, genes):
    """Build the output row for one Part 2 row (one input and one mouse ortholog).

    Copies the identifiers and statuses from Part 2, adds each species'
    UniProt and Bgee fields ("no_resolved_gene" where there is no gene on that
    side), the Jaccard index for each annotation set, the expression
    comparison and a Pair_Note that gathers every resolution, mapping and
    lookup note. A Jaccard index is left empty unless both species have
    at least one term of that kind.
    """
    result = {k: row[k] for k in ("Input_Row", "Human_Gene_ID", "Human_Symbol", "Mouse_Gene_ID", "Mouse_Symbol",
                                  "Resolution_Status", "Mapping_Status")}
    result.update(Dataset=dataset, Query_Symbol=row["Human_Input_ID"])
    notes = [row["Resolution_Note"]] if row["Resolution_Status"] != "resolved" else []
    if row["Mapping_Status"] not in ("ortholog_found", "multiple_mouse_orthologs") and row["Mapping_Note"]:
        notes.append(row["Mapping_Note"])
    human = genes.get(("Human", row["Human_Gene_ID"]), {})
    mouse = genes.get(("Mouse", row["Mouse_Gene_ID"]), {})
    for side, gene in (("Human", human), ("Mouse", mouse)):
        for field in SIDE_FIELDS:
            result[f"{side}_{field}"] = gene.get(field, "no_resolved_gene" if field.endswith("_Status") else "")
        notes.extend(f"{side} {note}" for note in gene.get("notes", []))
    for kind in SETS:
        left = human.get("sets", {}).get(kind, set())
        right = mouse.get("sets", {}).get(kind, set())
        result[kind + "_Jaccard"] = len(left & right) / len(left | right) if left and right else None
    result.update(compare_expression(human, mouse), Pair_Note="; ".join(notes))
    return result


def process(data_root, dataset):
    """Annotate every gene in one dataset's Part 2 output and write both output files.

    Returns True (the run exits with code 2) if any gene lookup left a note
    or incomplete Bgee data, or if Part 2 carried a failed Part 1 lookup or a
    mapping problem. p3_results.json is rewritten on every run, so it always
    matches the TSV. After a clean run, Part 3 TSV files from earlier script
    versions are deleted so they cannot be mistaken for current output.
    """
    # Part 2 writes standard quoted TSV.
    mappings = read_table(part_dir(data_root, 2, dataset) / "ncbi_ortholog.tsv")
    # Collect each distinct (side, Gene ID) once, merging its Ensembl IDs from every row it appears in.
    ids = defaultdict(set)
    for row in mappings:
        for side in TAXA:
            gid = row[side + "_Gene_ID"]
            if gid and row["Resolution_Status"] == "resolved":
                ids[side, gid].update(v for v in row[side + "_Ensembl_Gene_ID"].split(";") if v)
    genes = {(side, gid): annotate(gid, side, sorted(values)) for (side, gid), values in ids.items()}
    summaries = [pair_summary(dataset, row, genes) for row in mappings]
    folder = part_dir(data_root, 3, dataset)
    write_tsv(folder / "p3_pair_summary.tsv", COLUMNS, summaries)
    # Per-gene and per-protein detail for the loader's p3_gene_status and p3_uniprot tables.
    gene_rows = [{"Species": side, "Tax_ID": g["Tax_ID"], "Gene_ID": gid, "Gene_Type": g["Gene_Type"],
                  "Ensembl_IDs": ";".join(g["Ensembl_IDs"]), "UniProt_Status": g.get("UniProt_Status", ""),
                  "Bgee_Status": g.get("Bgee_Status", ""), "Notes": "; ".join(g["notes"])}
                 for (side, gid), g in genes.items()]
    protein_rows = [{"Species": side, "Tax_ID": g["Tax_ID"], "Gene_ID": gid,
                     "UniProt_Accession": e["primaryAccession"], "Protein_Name": protein_name(e),
                     "Reviewed": e.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"}
                    for (side, gid), g in genes.items() for e in g.get("entries", [])]
    (folder / "p3_results.json").write_text(
        json.dumps({"genes": gene_rows, "proteins": protein_rows}, indent=2, ensure_ascii=False), encoding="utf-8")
    failed = (any(g["notes"] or g["Bgee_Status"] == "incomplete" for g in genes.values())
              or any(r["Resolution_Status"] == "failed" or r["Mapping_Status"] not in MAPPING_OUTCOMES for r in mappings))
    if not failed:
        for old in ("p3_tissue_comparison.tsv", "p3_gene_status.tsv", "p3_uniprot.tsv", "p3_bgee_calls.tsv"):
            (folder / old).unlink(missing_ok=True)
    print(f"Part 3 {dataset}: {len(summaries)} input/ortholog rows saved.")
    return failed


if __name__ == "__main__":
    run(process)
