#!/usr/bin/env python3
import os
from getpass import getpass, getuser
from html import escape
from pathlib import Path
from urllib.parse import quote

import mysql.connector

OUT_DIR = Path(os.environ.get("BDB_HTML_DIR", "~/BDB_ICA/OP3_HTML")).expanduser()
SOCKET = os.environ.get("MYSQL_SOCKET", "/var/run/mysqld/mysqld.sock")
MY_CNF = Path("~/.my.cnf").expanduser()
DATASETS = ("development", "challenge")
TISSUES_SHOWN = 4

# One row per human input and mouse ortholog: Part 1 (NCBI), Part 2 (orthologs),
# Part 3 (UniProt + Bgee pair summary) and Part 4 (local API check).
# "<=>" is NULL-safe equals, so inputs with no mouse ortholog keep their Part 3 row.
SUMMARY_SQL = """
SELECT a.Input_Row, a.Query_Symbol, a.Resolution_Status, a.Resolution_Note,
       a.Gene_ID, a.Official_Symbol, a.Official_Name, a.Map_Location, a.ENSG_IDs,
       (SELECT MIN(d.Input_Row) FROM p1_annotations d
         WHERE d.Dataset = a.Dataset AND d.Gene_ID = a.Gene_ID) AS First_Input,
       o.Mouse_Gene_ID, o.Mouse_Symbol, o.Mouse_Official_Name, o.Mapping_Status,
       s.Human_UniProt_Accessions, s.Mouse_UniProt_Accessions, s.Protein_Comparison_Status,
       s.GO_MF_Jaccard, s.GO_BP_Jaccard, s.GO_CC_Jaccard, s.Location_Jaccard,
       s.Expression_Comparison_Status, s.Comparable_Tissue_Count, s.Expression_Call_Agreement,
       s.Present_Tissue_Jaccard, s.Shared_Present_Tissues, s.Pair_Note,
       l.Comparison_Status AS Local_Comparison, l.Local_Note, l.Differences
FROM p1_annotations a
LEFT JOIN p2_orthologs o    ON o.Dataset = a.Dataset AND o.Input_Row = a.Input_Row
LEFT JOIN p3_pair_summary s ON s.Dataset = a.Dataset AND s.Input_Row = a.Input_Row
                           AND s.Mouse_Gene_ID <=> o.Mouse_Gene_ID
LEFT JOIN p4_local_ncbi_comparisons l ON l.Dataset = a.Dataset AND l.Input_Row = a.Input_Row
WHERE a.Dataset = %s
ORDER BY a.Input_Row, o.Mouse_Gene_ID
"""


def e(value):
    return escape("" if value is None else str(value))


def link(url, text):
    return f'<a href="{e(url)}" target="_blank">{e(text)}</a>'


def na(reason):
    return f'<span class="na">N/A: {e(reason)}</span>'


def score(label, value):
    if value is None:
        return f"{label} {na('no terms')}<br>"
    return f"{label} <b>{float(value):.2f}</b><br>"


def reason(r):
    """Short reason a Part 3 comparison is missing."""
    if not r["Gene_ID"]:
        return "input not resolved"
    if not r["Mouse_Gene_ID"]:
        return "no mouse ortholog"
    if "not_applicable_gene_type" in (r["Pair_Note"] or ""):
        return "non-coding gene: no protein expected"
    return r["Pair_Note"] or "not compared"


def input_cell(r):
    html = f'<b>{e(r["Query_Symbol"])}</b> #{r["Input_Row"]}<br>'
    if r["Resolution_Status"] != "resolved":
        return html + na(r["Resolution_Note"] or "not resolved")
    if r["First_Input"] != r["Input_Row"]:
        html += f'<span class="small">Same gene as input #{r["First_Input"]}</span>'
    return html


def gene_cell(gene_id, symbol, name, extra=""):
    url = "https://www.ncbi.nlm.nih.gov/gene/" + quote(str(gene_id))
    return f'{link(url, symbol)}<br><span class="small">{e(name)}<br>{e(extra)}</span>'


def uniprot_cell(r):
    if r["Protein_Comparison_Status"] != "annotation_overlap_only":
        return na(reason(r))
    lines = []
    for species in ("Human", "Mouse"):
        accessions = json.loads(r[species + "_UniProt_Accessions"] or "[]")
        lines.append(f"{species}: " + ", ".join(
            link("https://www.uniprot.org/uniprotkb/" + quote(a), a) for a in accessions))
    return "<br>".join(lines)


def go_cell(r):
    if r["Protein_Comparison_Status"] != "annotation_overlap_only":
        return na(reason(r))
    return (score("MF", r["GO_MF_Jaccard"]) + score("BP", r["GO_BP_Jaccard"])
            + score("CC", r["GO_CC_Jaccard"]) + score("Location", r["Location_Jaccard"]))


def expression_cell(r):
    if r["Expression_Comparison_Status"] != "comparable":
        return na(reason(r))
    tissues = [t["tissue_name"] for t in json.loads(r["Shared_Present_Tissues"] or "[]")]
    shown = "; ".join(tissues[:TISSUES_SHOWN]) + (" ..." if len(tissues) > TISSUES_SHOWN else "")
    bgee = link("https://www.bgee.org/gene/" + quote(str(r["ENSG_IDs"]).split(",")[0].strip()),
                "Bgee (human)")
    return (score("Call agreement", r["Expression_Call_Agreement"])
            + f'<span class="small">over {r["Comparable_Tissue_Count"]} tissues</span><br>'
            + score("Present-tissue J", r["Present_Tissue_Jaccard"])
            + f'<span class="small">{len(tissues)} present in both: {e(shown)}</span><br>{bgee}')


def local_cell(r):
    status = r["Local_Comparison"]
    if status == "all_fields_agree":
        return '<span class="ok">agrees with NCBI</span>'
    if not status or status == "not_assessed":
        return na(r["Local_Note"] or "not assessed")
    return (f'<span class="warn">{e(status.replace("_", " "))}</span><br>'
            f'<span class="small">{e(r["Differences"])}</span>')


def row_html(r):
    if r["Gene_ID"]:
        human = gene_cell(r["Gene_ID"], r["Official_Symbol"], r["Official_Name"], r["Map_Location"])
    else:
        human = na("no NCBI Gene record")
    if r["Mouse_Gene_ID"]:
        mouse = gene_cell(r["Mouse_Gene_ID"], r["Mouse_Symbol"], r["Mouse_Official_Name"])
    else:
        mouse = na(str(r["Mapping_Status"] or "not assessed").replace("_", " "))
    cells = (input_cell(r), human, mouse, uniprot_cell(r), go_cell(r),
             expression_cell(r), local_cell(r))
    return "<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>"


HEADERS = ("Input", "Human gene (NCBI)", "Mouse ortholog (NCBI)", "UniProt entries",
           "Annotation overlap (Jaccard)", "Bgee expression (human vs mouse)",
           "Local API vs NCBI")

PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>{title}</title><style>
body{{font-family:Arial,sans-serif;font-size:11pt;margin:20px}}
table{{border-collapse:collapse;width:100%}}
th,td{{border:1px solid #999;padding:6px;vertical-align:top;text-align:left}}
th{{background:#2c4a6b;color:#fff}}
.small{{font-size:10pt;color:#555}} .na{{color:#888;font-style:italic}}
.ok{{background:#d6f0d6}} .warn{{background:#fff0c2}}
</style></head><body>
<h2>{title}</h2>
<table><tr>{headers}</tr>
{rows}
</table>
<p class="small">Jaccard = shared terms / all terms in human and mouse (0 = none shared,
1 = identical). MF, BP, CC = GO molecular function, biological process, cellular component;
Location = UniProt subcellular location. Call agreement = share of tissues where Bgee gives
the same present/absent call in both species. N/A always gives the reason.
Sources: NCBI Gene, NCBI orthologs, UniProtKB, Bgee, local course API.</p>
</body></html>
"""


def connect():
    settings = dict(unix_socket=SOCKET, charset="utf8mb4")
    if MY_CNF.is_file():
        settings.update(option_files=str(MY_CNF), option_groups=["client"])
    else:
        settings.update(user=getuser(), password=getpass("MySQL password: "))
    connection = mysql.connector.connect(**settings)
    if not connection.database:
        connection.database = os.environ.get("MYSQL_DATABASE") or input("Your database name: ")
    return connection


def main():
    connection = connect()
    cursor = connection.cursor(dictionary=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for dataset in DATASETS:
        cursor.execute(SUMMARY_SQL, (dataset,))
        rows = cursor.fetchall()
        title = f"Human vs mouse ortholog comparison: {dataset} gene list"
        page = PAGE.format(title=e(title), headers="".join(f"<th>{h}</th>" for h in HEADERS),
                           rows="\n".join(row_html(r) for r in rows))
        path = OUT_DIR / f"{dataset}_summary_table.html"
        path.write_text(page, encoding="utf-8")
        print(f"{dataset}: {len(rows)} rows -> {path}")
    connection.close()


if __name__ == "__main__":
    main()
