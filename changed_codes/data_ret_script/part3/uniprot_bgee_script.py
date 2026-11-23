#!/usr/bin/env python3
#Part 3: compare human-mouse orthologs by UniProt annotation and Bgee expression.
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

# Editable configuration: these paths match the combined Part 1/2 MySQL loader.
PROJECT_ROOT = Path(os.environ.get("BDB_PROJECT_ROOT", "~/BDB_ICA")).expanduser()
PART1_ROOT = PROJECT_ROOT / "OP1_data/part1_data/final_data"
PART2_ROOT = PROJECT_ROOT / "OP1_data/part2_data"
OUTPUT_ROOT = PROJECT_ROOT / "OP1_data/part3_data"
DATASETS = ("development", "challenge")
NCBI_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{}/dataset_report"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
ID_MAPPING_URL = "https://rest.uniprot.org/idmapping"
BGEE_URL = "https://www.bgee.org/api/"
TAXA = {"human": "9606", "mouse": "10090"}
TIMEOUT, MAX_ATTEMPTS, MIN_INTERVAL, CACHE_DAYS, MAX_POLLS = 30, 3, 0.4, 7, 30
RETRYABLE = {408, 429, 500, 502, 503, 504}
NONCODING_TYPES = {"NCRNA", "RRNA", "TRNA", "SNRNA", "SNORNA", "SCRNA", "NONCODING",
                   "LNCRNA", "MIRNA", "PSEUDO", "PSEUDOGENE", "BIOLOGICAL_REGION"}

# Output columns.
GENE_KEYS = ["Dataset", "Species", "Tax_ID", "Gene_ID"]
PAIR_KEYS = ["Dataset", "Input_Row", "Query_Symbol", "Human_Gene_ID", "Mouse_Gene_ID"]
GENE_COLUMNS = GENE_KEYS + [
    "Symbol", "Ensembl_IDs", "Gene_Type", "Gene_Type_Source", "Metadata_Status", "Metadata_Note",
    "UniProt_Status", "UniProt_Note", "Bgee_Status", "Bgee_Note"]
PROTEIN_COLUMNS = GENE_KEYS + [
    "UniProt_Accession", "Protein_Name", "Reviewed", "Selected", "Match_Method",
    "Function_Summary", "GO_MF_IDs", "GO_BP_IDs", "GO_CC_IDs", "Location_IDs",
    "Subcellular_Locations", "Function_Annotations_JSON", "GO_Annotations_JSON",
    "Location_Annotations_JSON", "Protein_Existence", "Annotation_Score", "Entry_Audit_JSON",
    "Source_URL"]
EXPRESSION_COLUMNS = GENE_KEYS + [
    "Bgee_Gene_ID", "Tissue_ID", "Tissue_Name", "Expression_Call", "Call_Quality",
    "Expression_Score", "Score_Confidence", "FDR", "FDR_Qualifier", "FDR_Reported", "Data_Types",
    "Condition_Scope", "Eligible_For_Comparison", "Source_URL"]
CALL_DETAIL = ["Call_Quality", "Expression_Score", "FDR", "FDR_Qualifier"]
TISSUE_COLUMNS = PAIR_KEYS + [
    "Tissue_ID", "Tissue_Name", "Human_Call", "Mouse_Call", "Comparison_Status"] + [
    f"{side}_{key}" for side in ("Human", "Mouse") for key in CALL_DETAIL]
ANNOTATION_SETS = ["GO_MF", "GO_BP", "GO_CC", "Location"]
SET_FIELDS = ["Human_Count", "Mouse_Count", "Shared_IDs", "Human_Only_IDs", "Mouse_Only_IDs",
              "Jaccard", "Status"]
SUMMARY_COLUMNS = PAIR_KEYS + [
    "Human_Symbol", "Mouse_Symbol", "Resolution_Status", "Resolution_Note", "Same_Gene_Note",
    "Mapping_Status", "Mapping_Note", "Human_Gene_Type", "Mouse_Gene_Type",
    "Human_UniProt_Status", "Mouse_UniProt_Status", "Human_UniProt_Accessions",
    "Mouse_UniProt_Accessions", "Protein_Comparison_Status"] + [
    f"{prefix}_{field}" for prefix in ANNOTATION_SETS for field in SET_FIELDS] + [
    "Human_Bgee_Status", "Mouse_Bgee_Status", "Expression_Comparison_Status",
    "Comparable_Tissue_Count", "Concordant_Tissue_Count", "Expression_Call_Agreement",
    "Shared_Present_Tissues", "Shared_Absent_Tissues", "Different_Call_Tissues",
    "Present_Tissue_Jaccard"]


class Part3Error(ValueError):
    """An input problem that must be fixed before any API request."""


class APIError(RuntimeError):
    def __init__(self, message, status=None):
        super().__init__(message)
        self.status = status


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_text(path, text):
    """Write to a temporary file, then rename, so no output or cache file is half-written."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=path.parent)
    with os.fdopen(handle, "w", encoding="utf-8", newline="") as stream:
        stream.write(text)
    os.replace(temporary, path)


def to_number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class JSONClient:
    """REST client with throttling, retries and a dated JSON cache (offline mode reads it only)."""

    def __init__(self, cache_dir, refresh=False, offline=False):
        self.cache_dir, self.refresh, self.offline = Path(cache_dir), refresh, offline
        self.last_request, self.audit = {}, []

    def _path(self, key):
        return self.cache_dir / (hashlib.sha256(key.encode()).hexdigest() + ".json")

    def load(self, key):
        path = self._path(key)
        if self.refresh or not path.exists():
            return None
        try:
            saved = json.loads(path.read_text(encoding="utf-8"))
            age = datetime.now(timezone.utc) - datetime.fromisoformat(saved["retrieved_utc"])
        except (ValueError, KeyError):
            return None  # Unreadable cache entry: fetch again.
        if not self.offline and age.days >= CACHE_DAYS:
            return None
        self.audit.append({"url": key, "retrieved_utc": saved["retrieved_utc"], "cached": True})
        return saved["response"], saved["headers"]

    def store(self, key, data, headers):
        retrieved = utc_now()
        atomic_text(self._path(key), json.dumps({"url": key, "retrieved_utc": retrieved,
                                                 "headers": headers, "response": data}))
        self.audit.append({"url": key, "retrieved_utc": retrieved, "cached": False})

    def get(self, base_url, params=None, form=None, cache=True):
        """GET, or POST when form is given. POSTs and temporary job URLs are not cached."""
        url = base_url + ("?" + urlencode(params) if params else "")
        cache = cache and not form
        saved = self.load(url) if cache else None
        if saved:
            return saved
        if self.offline:
            raise APIError("Offline mode: no cached response for " + url)
        host, problem = urlsplit(url).netloc, ""
        for attempt in range(MAX_ATTEMPTS):
            time.sleep(max(0, MIN_INTERVAL - (time.monotonic() - self.last_request.get(host, 0))))
            self.last_request[host] = time.monotonic()
            request = Request(url, data=urlencode(form).encode() if form else None,
                              headers={"Accept": "application/json", "User-Agent": "BDB_ICA"})
            try:
                with urlopen(request, timeout=TIMEOUT) as response:
                    data = json.loads(response.read())
                    headers = {k.lower(): v for k, v in response.headers.items()}
            except HTTPError as error:
                if error.code not in RETRYABLE:
                    raise APIError(f"HTTP {error.code}.", error.code) from None
                problem = f"HTTP {error.code}"
            except OSError as error:  # Connection failures and timeouts.
                problem = type(error).__name__
            except ValueError:
                raise APIError("The service returned invalid JSON.") from None
            else:
                if not isinstance(data, dict) or data.get("error") or data.get("errors") or str(
                        data.get("code", 200)) != "200":
                    raise APIError("The service returned an error response.")
                if cache:
                    self.store(url, data, headers)
                return data, headers
            if attempt < MAX_ATTEMPTS - 1:
                time.sleep(2 ** attempt)
        raise APIError(f"{problem} after {MAX_ATTEMPTS} attempts.")


def read_tsv(path, plain=False):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t",
                                quoting=csv.QUOTE_NONE if plain else csv.QUOTE_MINIMAL)
        rows = [{k: (v or "").strip() for k, v in row.items() if k} for row in reader]
        return reader.fieldnames or [], rows


def ensembl_ids(value, species):
    """Validate and de-version Ensembl gene IDs given as a list or delimited string."""
    prefix = "ENSG" if species == "human" else "ENSMUSG"
    tokens = value if isinstance(value, list) else re.split(r"[;,\s]+", value or "")
    result = set()
    for token in filter(None, tokens):
        if not re.fullmatch(prefix + r"\d+(\.\d+)?", token):
            raise Part3Error(f"Invalid {species} Ensembl gene ID: {token!r}.")
        result.add(token.split(".")[0])
    return sorted(result)


def load_inputs(part1_path, part2_path, dataset):
    """One pair row per Part 1 input and Part 2 mouse ortholog (joined on Input_Row)."""
    f1, part1 = read_tsv(part1_path, plain=True)
    f2, part2 = read_tsv(part2_path)
    for name, fields, needed in [("Part 1", f1, {"Query_Symbol", "Gene_ID"}),
                                 ("Part 2", f2, {"Input_Row", "Human_Gene_ID", "Mouse_Gene_ID"})]:
        if not needed <= set(fields):
            raise Part3Error(f"{name} is missing columns {sorted(needed - set(fields))}.")
    humans = {}
    for number, row in enumerate(part1, 1):
        row["Input_Row"] = row.get("Input_Row") or str(number)
        row["Resolution_Status"] = row.get("Resolution_Status") or (
            "resolved" if row["Gene_ID"] else "unresolved")
        humans[row["Input_Row"]] = row
    if not humans or len(humans) != len(part1):
        raise Part3Error("Part 1 is empty or has duplicate Input_Row values.")
    orthologs = defaultdict(list)
    for row in part2:
        human = humans.get(row["Input_Row"])
        if human is None or human["Gene_ID"] != row["Human_Gene_ID"]:
            raise Part3Error(f"Part 2 row {row['Input_Row']} does not match Part 1; rerun Part 2.")
        orthologs[row["Input_Row"]].append(row)
    pairs = []
    for key, human in humans.items():
        for mouse in orthologs[key] or [{"Mapping_Status": "part2_row_missing"}]:
            mouse_id = mouse.get("Mouse_Gene_ID", "")
            pairs.append({
                "Dataset": dataset, "Input_Row": key, "Query_Symbol": human["Query_Symbol"],
                "Human_Gene_ID": human["Gene_ID"], "Human_Symbol": human.get("Official_Symbol", ""),
                "Mouse_Gene_ID": mouse_id, "Mouse_Symbol": mouse.get("Mouse_Symbol", ""),
                "Mapping_Status": mouse.get("Mapping_Status")
                or ("ortholog_found" if mouse_id else "unspecified_no_mouse_row"),
                "Mapping_Note": mouse.get("Mapping_Note", ""),
                **{k: human.get(k, "") for k in ("Resolution_Status", "Resolution_Note",
                                                  "Same_Gene_Note")},
                "Human_Ensembl_IDs": ensembl_ids(
                    human.get("ENSG_IDs") or mouse.get("Human_Ensembl_Gene_ID"), "human"),
                "Mouse_Ensembl_IDs": ensembl_ids(mouse.get("Mouse_Ensembl_Gene_ID"), "mouse"),
                "Human_Gene_Type": human.get("Gene_Type", "")})
    return pairs, len(part1)


def ncbi_metadata(client, gene):
    """Confirm gene ID and species with NCBI; take its current gene type and Ensembl IDs."""
    reports = client.get(NCBI_URL.format(gene["Gene_ID"]))[0].get("reports")
    if not reports:
        raise APIError("NCBI returned no gene record.")
    record = reports[0].get("gene", {})
    if str(record.get("gene_id")) != gene["Gene_ID"] or str(record.get("tax_id")) != gene["Tax_ID"]:
        raise APIError("NCBI returned a different gene ID or species.", "identity_error")
    gene["Gene_Type"], gene["Gene_Type_Source"] = record.get("type") or "", "NCBI Datasets Gene"
    current = ensembl_ids(record.get("ensembl_gene_ids", []), gene["Species"])
    gene["Ensembl_IDs"] = current or gene["Ensembl_IDs"]
    gene["Symbol"] = record.get("symbol") or gene["Symbol"]
    gene["Metadata_Status"] = "retrieved"


def uniprot_pages(client, url, params, cache=True):
    """Yield every page of a UniProt result set by following next-page Link headers."""
    for _ in range(100):
        data, headers = client.get(url, params, cache=cache)
        yield data
        link = re.search(r'<([^>]+)>;\s*rel="?next"?', headers.get("link", ""))
        if not link:
            return
        url, params = link.group(1), None
    raise APIError("UniProt exceeded the pagination limit.")


def map_ensembl_ids(client, ids):
    """UniProt ID mapping (Ensembl gene -> UniProtKB); only the finished result is cached."""
    form = {"from": "Ensembl", "to": "UniProtKB", "ids": ",".join(ids)}
    key = ID_MAPPING_URL + "/run?" + urlencode(form)
    saved = client.load(key)
    if saved:
        return saved[0]["results"]
    job = client.get(ID_MAPPING_URL + "/run", form=form)[0]["jobId"]
    for _ in range(MAX_POLLS):
        status = client.get(f"{ID_MAPPING_URL}/status/{job}", cache=False)[0].get("jobStatus")
        if status not in ("NEW", "RUNNING"):
            break
        time.sleep(2)
    # A finished job reports FINISHED or redirects straight to its results (no jobStatus).
    if status not in (None, "FINISHED"):
        raise APIError(f"UniProt ID mapping job status: {status}.")
    results_url = f"{ID_MAPPING_URL}/uniprotkb/results/{job}"
    rows = [row for page in uniprot_pages(client, results_url, {"format": "json", "size": 500},
                                          cache=False) for row in page.get("results", [])]
    client.store(key, {"results": rows}, {})
    return rows


def gene_links(entry):
    """NCBI GeneIDs and Ensembl gene IDs cross-referenced by a UniProt entry."""
    xrefs = entry.get("uniProtKBCrossReferences", [])
    ncbi = {x.get("id") for x in xrefs if x.get("database") == "GeneID"}
    ensembl = {p.get("value", "").split(".")[0] for x in xrefs if x.get("database") == "Ensembl"
               for p in x.get("properties", []) if p.get("key") == "GeneId"}
    return ncbi, ensembl


def parse_protein(entry, gene, method):
    """Keep annotation IDs plus the evidence attached to each annotation."""
    accession, description = entry["primaryAccession"], entry.get("proteinDescription", {})
    names = [description.get("recommendedName", {})] + description.get("submissionNames", [])
    name = next((n["fullName"]["value"] for n in names if n.get("fullName")), "")
    comments = entry.get("comments", [])
    functions = [t for c in comments if c.get("commentType") == "FUNCTION"
                 for t in c.get("texts", [])]
    locations = [{"location": item["location"], "topology": item.get("topology"),
                  "orientation": item.get("orientation"), "molecule": c.get("molecule"),
                  "note": c.get("note")}
                 for c in comments if c.get("commentType") == "SUBCELLULAR LOCATION"
                 for item in c.get("subcellularLocations", []) if item.get("location")]
    go = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            props = {p["key"]: p.get("value", "") for p in xref.get("properties", [])}
            term = props.get("GoTerm", "")  # e.g. "F:ATP binding"
            go.append({"id": xref["id"], "aspect": term[:1], "term": term[2:],
                       "evidence_type": props.get("GoEvidenceType", ""),
                       "evidences": xref.get("evidences", [])})
    go_ids = {aspect: sorted({t["id"] for t in go if t["aspect"] == letter})
              for aspect, letter in (("MF", "F"), ("BP", "P"), ("CC", "C"))}
    return {
        **{k: gene[k] for k in GENE_KEYS}, "UniProt_Accession": accession, "Protein_Name": name,
        "Reviewed": entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "Selected": False, "Match_Method": method,
        "Function_Summary": " ".join(t.get("value", "") for t in functions),
        "GO_MF_IDs": go_ids["MF"], "GO_BP_IDs": go_ids["BP"], "GO_CC_IDs": go_ids["CC"],
        "Location_IDs": sorted({l["location"]["id"] for l in locations if "id" in l["location"]}),
        "Subcellular_Locations": sorted({l["location"].get("value", "") for l in locations}),
        "Function_Annotations_JSON": functions, "GO_Annotations_JSON": go,
        "Location_Annotations_JSON": locations,
        "Protein_Existence": entry.get("proteinExistence", ""),
        "Annotation_Score": entry.get("annotationScore", ""),
        "Entry_Audit_JSON": entry.get("entryAudit", {}),
        "Source_URL": "https://www.uniprot.org/uniprotkb/" + accession}


def fetch_uniprot(client, gene):
    """Exact GeneID cross-reference first; Ensembl gene ID mapping only if that finds nothing."""
    if gene["Gene_Type"].upper().replace("-", "_").replace(" ", "_") in NONCODING_TYPES:
        gene["UniProt_Status"] = "not_applicable_gene_type"
        gene["UniProt_Note"] = f"NCBI gene type {gene['Gene_Type']}; no protein expected."
        return []
    taxon, gene_id = gene["Tax_ID"], gene["Gene_ID"]
    method = "exact_GeneID_cross_reference"
    query = {"query": f"organism_id:{taxon} AND xref:geneid-{gene_id}", "format": "json",
             "size": 500}
    candidates = [(gene_id, e) for page in uniprot_pages(client, UNIPROT_URL, query)
                  for e in page["results"]]
    if not candidates and gene["Ensembl_IDs"]:
        method = "Ensembl_gene_ID_mapping"
        candidates = [(r["from"], r["to"]) for r in map_ensembl_ids(client, gene["Ensembl_IDs"])]
    proteins = {}
    for source, entry in candidates:
        ncbi, ensembl = gene_links(entry)
        # Keep entries of the right species that do not point to a different gene.
        verified = str(entry.get("organism", {}).get("taxonId")) == taxon and (
            gene_id in ncbi if source == gene_id else
            (not ncbi or gene_id in ncbi) and (not ensembl or source in ensembl))
        if verified:
            proteins.setdefault(entry["primaryAccession"], parse_protein(entry, gene, method))
    proteins = list(proteins.values())
    reviewed = any(p["Reviewed"] for p in proteins)
    for protein in proteins:  # All reviewed entries, else all verified unreviewed entries.
        protein["Selected"] = protein["Reviewed"] or not reviewed
    gene["UniProt_Status"] = ("reviewed_found" if reviewed else "unreviewed_only" if proteins
                              else "no_verified_entry")
    gene["UniProt_Note"] = f"Match route: {method}." if proteins else (
        f"{len(candidates)} candidate entries; none verified.")
    return proteins


def bgee_row(call, gene, ensembl):
    anatomy = call["condition"]["anatEntity"]
    quality = str(call.get("expressionQuality") or "").lower()
    score = call.get("expressionScore") or {}
    reported = str(call.get("fdr") or "").strip()  # e.g. "<= 1.00e-14": keep the inequality
    qualifier, value = re.fullmatch(r"(<=|>=|<|>|≤|≥)?\s*(.*)", reported).groups()
    fdr = to_number(value)
    return {
        **{k: gene[k] for k in GENE_KEYS}, "Bgee_Gene_ID": ensembl,
        "Tissue_ID": anatomy["id"], "Tissue_Name": anatomy.get("name", ""),
        "Expression_Call": {"expressed": "present", "not_expressed": "absent"}[
            call["expressionState"]],
        "Call_Quality": quality, "Expression_Score": to_number(score.get("expressionScore")),
        "Score_Confidence": score.get("expressionScoreConfidence", ""),
        "FDR": fdr, "FDR_Reported": reported,
        "FDR_Qualifier": "" if fdr is None else {"≤": "<=", "≥": ">="}.get(qualifier,
                                                                           qualifier or "="),
        "Data_Types": call.get("dataTypesWithData", []),
        "Condition_Scope": "anatomy only; stages/sexes/strains aggregated",
        "Eligible_For_Comparison": anatomy["id"].startswith("UBERON:")
        and quality in {"gold", "silver"},
        "Source_URL": "https://www.bgee.org/gene/" + ensembl}


def fetch_bgee(client, gene):
    """Present and absent calls per anatomical entity, for every Ensembl ID of the gene."""
    if not gene["Ensembl_IDs"]:
        gene["Bgee_Status"], gene["Bgee_Note"] = "missing_ensembl_id", "No Ensembl gene ID."
        return []
    rows, outcomes, notes = [], [], []
    for ensembl in gene["Ensembl_IDs"]:
        for expr_type in ("EXPRESSED", "NOT_EXPRESSED"):
            params = {"page": "gene", "action": "expression", "gene_id": ensembl,
                      "species_id": gene["Tax_ID"], "display_type": "json",
                      "cond_param": "anat_entity", "expr_type": expr_type}
            try:
                calls = client.get(BGEE_URL, params)[0]["data"]["calls"]
                rows.extend([bgee_row(call, gene, ensembl) for call in calls])
                outcomes.append("retrieved")
            except (APIError, KeyError, TypeError, AttributeError, ValueError) as error:
                outcomes.append("not_found" if getattr(error, "status", None) == 404
                                else "lookup_failed")
                notes.append(f"{ensembl} {expr_type}: {error or type(error).__name__}")
    found = set(outcomes)
    gene["Bgee_Status"] = (("retrieved" if rows else "no_calls_returned")
                           if found == {"retrieved"} else "incomplete" if "retrieved" in found
                           else "not_found" if found == {"not_found"} else "lookup_failed")
    if len(gene["Ensembl_IDs"]) > 1:
        notes.append("Several Ensembl IDs; tissue comparison skipped.")
    gene["Bgee_Note"] = "; ".join(notes)
    return rows


def collect_genes(client, pairs, dataset):
    """Look up each distinct resolved human and mouse gene once."""
    genes = {}
    for pair in pairs:
        for species, tax_id in TAXA.items():
            side = species.title()
            gene_id = pair[side + "_Gene_ID"]
            if not gene_id or pair["Resolution_Status"] != "resolved":
                continue
            gene_type = pair.get(side + "_Gene_Type", "")
            gene = genes.setdefault((species, gene_id), {
                "Dataset": dataset, "Species": species, "Tax_ID": tax_id, "Gene_ID": gene_id,
                "Symbol": pair[side + "_Symbol"], "Ensembl_IDs": [], "Gene_Type": gene_type,
                "Gene_Type_Source": "Part 1" if gene_type else "",
                **{f"{step}_{field}": value for step in ("Metadata", "UniProt", "Bgee")
                   for field, value in (("Status", "not_requested"), ("Note", ""))}})
            gene["Ensembl_IDs"] = sorted({*gene["Ensembl_IDs"], *pair[side + "_Ensembl_IDs"]})
    proteins, expression = [], []
    for gene in genes.values():
        print(f"[INFO] {dataset}: {gene['Species']} {gene['Symbol'] or gene['Gene_ID']}",
              flush=True)
        try:
            ncbi_metadata(client, gene)
        except (APIError, Part3Error, KeyError, TypeError, AttributeError) as error:
            gene["Metadata_Status"], gene["Metadata_Note"] = "lookup_failed", str(error)
            if getattr(error, "status", None) == "identity_error":
                gene["UniProt_Status"] = gene["Bgee_Status"] = "invalid_gene_identity"
                continue
        try:
            proteins.extend(fetch_uniprot(client, gene))
        except (APIError, KeyError, TypeError, AttributeError, ValueError) as error:
            gene["UniProt_Status"], gene["UniProt_Note"] = "lookup_failed", str(error)
        expression.extend(fetch_bgee(client, gene))
    return genes, proteins, expression


def set_comparison(left, right):
    """Jaccard of two ID sets; blank (not zero) when either side has no annotation."""
    shared, both = left & right, bool(left and right)
    return {"Human_Count": len(left), "Mouse_Count": len(right), "Shared_IDs": sorted(shared),
            "Human_Only_IDs": sorted(left - right), "Mouse_Only_IDs": sorted(right - left),
            "Jaccard": len(shared) / len(left | right) if both else None,
            "Status": "comparable" if both else "missing_annotation"}


def compare_proteins(human, mouse, proteins):
    sides = [("Human", human), ("Mouse", mouse)]
    output = {f"{side}_UniProt_Status": g["UniProt_Status"] if g else "no_resolved_gene"
              for side, g in sides}
    found = {"reviewed_found", "unreviewed_only"}
    if not (human and mouse and {g["UniProt_Status"] for _, g in sides} <= found):
        output["Protein_Comparison_Status"] = "not_comparable"
        return output
    chosen = {side: [p for p in proteins if p["Selected"] and p["Species"] == g["Species"]
                     and p["Gene_ID"] == g["Gene_ID"]] for side, g in sides}
    for side, _ in sides:
        output[side + "_UniProt_Accessions"] = [p["UniProt_Accession"] for p in chosen[side]]
    for prefix in ANNOTATION_SETS:  # Gene level: union over all selected proteins.
        sets = [{i for p in chosen[side] for i in p[prefix + "_IDs"]} for side, _ in sides]
        output.update({f"{prefix}_{k}": v for k, v in set_comparison(*sets).items()})
    output["Protein_Comparison_Status"] = "annotation_overlap_only"
    return output


def index_expression(rows):
    """Eligible calls by (species, gene) and tissue; a tissue seen twice is marked None."""
    index = defaultdict(dict)
    for row in rows:
        if row["Eligible_For_Comparison"]:
            calls = index[(row["Species"], row["Gene_ID"])]
            calls[row["Tissue_ID"]] = None if row["Tissue_ID"] in calls else row
    return index


def compare_expression(pair, human, mouse, index):
    summary = {"Comparable_Tissue_Count": 0, "Concordant_Tissue_Count": 0,
               "Expression_Call_Agreement": None, "Shared_Present_Tissues": [],
               "Shared_Absent_Tissues": [], "Different_Call_Tissues": [],
               "Present_Tissue_Jaccard": None,
               "Human_Bgee_Status": human["Bgee_Status"] if human else "no_resolved_gene",
               "Mouse_Bgee_Status": mouse["Bgee_Status"] if mouse else "no_resolved_gene"}
    detail = []

    def stop(status):
        summary["Expression_Comparison_Status"] = status
        return summary, detail

    if not human or not mouse:
        return stop("not_comparable")
    if len(human["Ensembl_IDs"]) != 1 or len(mouse["Ensembl_IDs"]) != 1:
        return stop("missing_or_multiple_ensembl_ids")
    if any(g["Bgee_Status"] not in {"retrieved", "no_calls_returned"} for g in (human, mouse)):
        return stop("incomplete_expression_lookup")
    hcalls = index.get(("human", human["Gene_ID"]), {})
    mcalls = index.get(("mouse", mouse["Gene_ID"]), {})
    if not hcalls and not mcalls:
        return stop("no_eligible_tissues")

    def label(calls, tissue):
        return (calls[tissue]["Expression_Call"] if calls.get(tissue)
                else "conflicting_calls" if tissue in calls else "no_data")

    both = []  # Tissues with one eligible call in each species: (tissue, human, mouse).
    for tissue in sorted({*hcalls, *mcalls}):
        h, m = hcalls.get(tissue), mcalls.get(tissue)
        if (tissue in hcalls and h is None) or (tissue in mcalls and m is None):
            status = "conflicting_calls"
        elif h and m:
            both.append((tissue, h["Expression_Call"], m["Expression_Call"]))
            status = "same_call" if both[-1][1] == both[-1][2] else "different_call"
        else:
            status = "no_comparable_call_in_one_species"
        row = {k: pair[k] for k in PAIR_KEYS}
        row.update({"Tissue_ID": tissue, "Tissue_Name": (h or m or {}).get("Tissue_Name", ""),
                    "Human_Call": label(hcalls, tissue), "Mouse_Call": label(mcalls, tissue),
                    "Comparison_Status": status})
        for side, call in (("Human", h), ("Mouse", m)):
            row.update({f"{side}_{k}": call.get(k) if call else None for k in CALL_DETAIL})
        detail.append(row)
    if not both:
        return stop("no_shared_comparable_tissues")
    present = [t for t, h, m in both if h == m == "present"]
    absent = [t for t, h, m in both if h == m == "absent"]
    any_present = [t for t, h, m in both if "present" in (h, m)]
    summary.update({
        "Comparable_Tissue_Count": len(both), "Concordant_Tissue_Count": len(present + absent),
        "Expression_Call_Agreement": len(present + absent) / len(both),
        "Shared_Present_Tissues": present, "Shared_Absent_Tissues": absent,
        "Different_Call_Tissues": [t for t, h, m in both if h != m],
        "Present_Tissue_Jaccard": len(present) / len(any_present) if any_present else None})
    return stop("comparable")


def write_tsv(path, rows, columns):
    """Lists and dicts are written as compact JSON text, booleans as 0/1."""
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({
            k: json.dumps(v, ensure_ascii=False, separators=(",", ":"))
            if isinstance(v, (dict, list)) else int(v) if isinstance(v, bool) else v
            for k, v in row.items()})
    atomic_text(path, buffer.getvalue())


def run_dataset(pairs, input_count, dataset, output_dir, source_paths, client):
    genes, proteins, expression = collect_genes(client, pairs, dataset)
    index = index_expression(expression)
    summaries, tissues = [], []
    for pair in pairs:
        resolved = pair["Resolution_Status"] == "resolved"
        human, mouse = (genes.get((s, pair[s.title() + "_Gene_ID"])) if resolved else None
                        for s in TAXA)
        row = {k: v for k, v in pair.items() if k in SUMMARY_COLUMNS}
        row["Human_Gene_Type"] = human["Gene_Type"] if human else ""
        row["Mouse_Gene_Type"] = mouse["Gene_Type"] if mouse else ""
        row.update(compare_proteins(human, mouse, proteins))
        expression_summary, expression_detail = compare_expression(pair, human, mouse, index)
        summaries.append({**row, **expression_summary})
        tissues.extend(expression_detail)
    tables = {"p3_gene_status.tsv": (list(genes.values()), GENE_COLUMNS),
              "p3_uniprot.tsv": (proteins, PROTEIN_COLUMNS),
              "p3_bgee_calls.tsv": (expression, EXPRESSION_COLUMNS),
              "p3_tissue_comparison.tsv": (tissues, TISSUE_COLUMNS),
              "p3_pair_summary.tsv": (summaries, SUMMARY_COLUMNS)}
    for filename, (rows, columns) in tables.items():
        write_tsv(Path(output_dir) / filename, rows, columns)
    problems = [g for g in genes.values() if "lookup_failed" in (g["Metadata_Status"],
                g["UniProt_Status"], g["Bgee_Status"]) or g["Bgee_Status"] == "incomplete"
                or g["UniProt_Status"] == "invalid_gene_identity"]
    manifest = {
        "dataset": dataset, "generated_utc": utc_now(), "inputs": [os.path.relpath(p, PROJECT_ROOT) for p in source_paths],
        "input_count": input_count, "pair_count": len(pairs), "unique_gene_count": len(genes),
        "genes_with_lookup_problems": len(problems),
        "files": {name: len(rows) for name, (rows, _) in tables.items()},
        "protein_comparison_status_counts":
            dict(Counter(r["Protein_Comparison_Status"] for r in summaries)),
        "expression_comparison_status_counts":
            dict(Counter(r["Expression_Comparison_Status"] for r in summaries))}
    results = {"manifest": manifest, "genes": list(genes.values()), "proteins": proteins,
               "expression_calls": expression, "tissue_comparisons": tissues,
               "pair_summaries": summaries, "request_provenance": client.audit}
    atomic_text(Path(output_dir) / "p3_results.json",
                json.dumps(results, ensure_ascii=False, indent=2))
    print(f"[SAVED] {dataset}: {input_count} inputs, {len(pairs)} pairs -> {output_dir}")
    if problems:
        print(f"[WARN] {len(problems)} genes had lookup problems; see p3_gene_status.tsv.")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", choices=["both", *DATASETS], default="both")
    parser.add_argument("--part1-root", type=Path, default=PART1_ROOT)
    parser.add_argument("--part2-root", type=Path, default=PART2_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--check-inputs", action="store_true", help="Validate inputs only")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--refresh", action="store_true", help="Ignore the dated API cache")
    mode.add_argument("--offline", action="store_true", help="Use cached API responses only")
    args = parser.parse_args(argv)
    output_root = args.output_root.expanduser().resolve()
    jobs = []
    for dataset in DATASETS if args.run == "both" else (args.run,):
        sources = [(args.part1_root / dataset / "ncbi_gene_summary2.tsv").expanduser(),
                   (args.part2_root / dataset / "ncbi_ortholog.tsv").expanduser()]
        pairs, count = load_inputs(*sources, dataset)
        jobs.append((pairs, count, dataset, output_root / dataset, sources))
        print(f"[CHECKED] {dataset}: {count} inputs -> {len(pairs)} pair rows.")
    if args.check_inputs:
        return 0
    client = JSONClient(output_root / "api_cache", args.refresh, args.offline)
    failed = [run_dataset(*job, client)["genes_with_lookup_problems"] for job in jobs]
    return 2 if any(failed) else 0  # 2 = outputs saved, but some lookups need attention.


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("[ERROR] Cancelled; finished outputs and cached responses are kept.",
              file=sys.stderr)
        sys.exit(130)
    except Exception as error:  # Input, file or unexpected errors: no traceback.
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
