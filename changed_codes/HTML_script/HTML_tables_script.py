#!/usr/bin/env python3
"""Build one HTML summary table per gene list from the 4 MySQL tables (SQL JOIN).

Run on bioinfmsc9 after all_sql_loads:
    python3 ~/BDB_ICA/changed_codes/HTML_script/HTML_tables_script.py
It asks for your MySQL password and database name, then writes
OP3_HTML/development_summary_table.html and OP3_HTML/challenge_summary_table.html.

GenAI declaration:
  Claude (Anthropic) was used to debug my earlier version of this script.
"""
import json
import os
from getpass import getpass, getuser
from html import escape
from pathlib import Path
from urllib.parse import quote

import mysql.connector

ROOT = Path(__file__).resolve().parents[2]          # repo root (~/BDB_ICA)
OUT_DIR = ROOT / "OP3_HTML"
SOCKET = "/var/run/mysqld/mysqld.sock"
DATASETS = ("development", "challenge")
TISSUES_SHOWN = 4

# One row per human input and mouse ortholog: Part 1 (NCBI), Part 2 (orthologs),
# Part 3 (UniProt + Bgee pair summary) and Part 4 (local API check).
# The SELECT is built at run time from the tables and columns that really exist:
# a missing column is selected as NULL and a missing table is left out of the JOIN,
# so the matching cells show "N/A" and the script prints a note instead of crashing.
# "<=>" is NULL-safe equals, so inputs with no mouse ortholog keep their Part 3 row.
TABLES = {  # alias: (table, columns the HTML uses, join keys, JOIN condition)
    "a": ("p1_annotations",
          ["Input_Row", "Query_Symbol", "Resolution_Status", "Resolution_Note", "Gene_ID",
           "Official_Symbol", "Official_Name", "Map_Location", "ENSG_IDs"],
          ["Dataset", "Input_Row"], None),
    "o": ("p2_orthologs",
          ["Mouse_Gene_ID", "Mouse_Symbol", "Mouse_Official_Name", "Mapping_Status"],
          ["Dataset", "Input_Row"], "o.Dataset = a.Dataset AND o.Input_Row = a.Input_Row"),
    "s": ("p3_pair_summary",
          ["Human_UniProt_Accessions", "Mouse_UniProt_Accessions", "Protein_Comparison_Status",
           "Human_UniProt_Status", "GO_MF_Jaccard", "GO_BP_Jaccard", "GO_CC_Jaccard",
           "Location_Jaccard", "Expression_Comparison_Status", "Comparable_Tissue_Count",
           "Expression_Call_Agreement", "Present_Tissue_Jaccard", "Shared_Present_Tissues",
           "Pair_Note"],
          ["Dataset", "Input_Row"], "s.Dataset = a.Dataset AND s.Input_Row = a.Input_Row"),
    "l": ("p4_local_ncbi_comparisons",
          ["Comparison_Status", "Local_Note", "Differences"],
          ["Dataset", "Input_Row"], "l.Dataset = a.Dataset AND l.Input_Row = a.Input_Row"),
}
RENAME = {("l", "Comparison_Status"): "Local_Comparison"}   # avoid clashing names
# Used when a column is missing but can be worked out from others in the same table.
DERIVED = {("s", "Protein_Comparison_Status"):
           "CASE WHEN s.Human_UniProt_Accessions IS NOT NULL AND s.Mouse_UniProt_Accessions"
           " IS NOT NULL THEN 'annotation_overlap_only' END"}


def table_columns(cursor, table):
    """Column names of a table, or None if the table does not exist."""
    cursor.execute("SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS"
                   " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s", (table,))
    names = {list(row.values())[0] for row in cursor.fetchall()}
    return names or None


def build_sql(cursor):
    """Build SUMMARY_SQL from what is in the database. Returns None if Part 1 is missing."""
    found = {alias: table_columns(cursor, spec[0]) for alias, spec in TABLES.items()}
    if not found["a"] or not {"Dataset", "Input_Row"} <= found["a"]:
        print("Note: p1_annotations (with Dataset and Input_Row) is missing; tables will be empty")
        return None
    fields, joins = [], []
    for alias, (table, columns, keys, condition) in TABLES.items():
        existing = found[alias] or set()
        usable = existing and set(keys) <= existing
        if alias == "s" and usable:                    # also match the mouse ortholog
            if "Mouse_Gene_ID" in existing and {"Dataset", "Input_Row", "Mouse_Gene_ID"} <= (found["o"] or set()):
                condition += " AND s.Mouse_Gene_ID <=> o.Mouse_Gene_ID"
        if alias != "a":
            if usable:
                joins.append(f"LEFT JOIN `{table}` {alias} ON {condition}")
            else:
                print(f"Note: table {table} is missing or has no {'/'.join(keys)}; showing N/A")
        for name in columns:
            label = RENAME.get((alias, name), name)
            if usable and name in existing:
                fields.append(f"{alias}.`{name}` AS `{label}`")
            elif usable and (alias, name) in DERIVED:
                fields.append(f"{DERIVED[(alias, name)]} AS `{label}`")
            else:
                if usable:
                    print(f"Note: {table} has no {name} column; showing N/A")
                fields.append(f"NULL AS `{label}`")
    if "Gene_ID" in found["a"]:
        fields.append("(SELECT MIN(d.Input_Row) FROM p1_annotations d"
                      " WHERE d.Dataset = a.Dataset AND d.Gene_ID = a.Gene_ID) AS First_Input")
    else:
        fields.append("a.Input_Row AS First_Input")
    order = "a.Input_Row" + (", o.Mouse_Gene_ID" if any(" o " in j for j in joins)
                             and "Mouse_Gene_ID" in (found["o"] or ()) else "")
    return ("SELECT " + ",\n       ".join(fields) + "\nFROM p1_annotations a\n"
            + "\n".join(joins) + f"\nWHERE a.Dataset = %s\nORDER BY {order}")


def e(value):
    return escape("" if value is None else str(value))


def link(url, text):
    return f'<a href="{e(url)}" target="_blank">{e(text)}</a>'


def na(reason):
    return f'<span class="na">N/A: {e(reason)}</span>'


def score(label, value):
    if value is None:
        return f"{label} {na('no terms')}<br>"
    try:
        return f"{label} <b>{float(value):.2f}</b><br>"
    except (TypeError, ValueError):
        return f"{label} {e(value)}<br>"


def json_list(text):
    """A JSON list, or a "; "-separated text list, as a Python list."""
    text = (text or "").strip()
    if text.startswith("["):
        try:
            return json.loads(text)
        except ValueError:
            pass
    return [part.strip() for part in text.split(";") if part.strip()]


def reason(r):
    """Short reason a Part 3 comparison is missing."""
    if not r["Gene_ID"]:
        return "input not resolved"
    if not r["Mouse_Gene_ID"]:
        return "no mouse ortholog"
    if (r["Human_UniProt_Status"] == "not_applicable_gene_type"
            or "not_applicable_gene_type" in (r["Pair_Note"] or "")):
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
        accessions = json_list(r[species + "_UniProt_Accessions"])
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
    tissues = [t.get("tissue_name", "") if isinstance(t, dict) else str(t)
               for t in json_list(r["Shared_Present_Tissues"])]
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


def safe(cell, r):
    """Run one cell function; if its data is in an unexpected form, show N/A, not a crash."""
    try:
        return cell(r)
    except Exception as error:  # noqa: BLE001 - one bad value must not stop the table
        print(f"Note: input #{r.get('Input_Row')} {cell.__name__}: {error!r}; showing N/A")
        return na("value could not be read")


def row_html(r):
    if r["Gene_ID"]:
        human = gene_cell(r["Gene_ID"], r["Official_Symbol"] or r["Gene_ID"], r["Official_Name"], r["Map_Location"])
    else:
        human = na("no NCBI Gene record")
    if r["Mouse_Gene_ID"]:
        mouse = gene_cell(r["Mouse_Gene_ID"], r["Mouse_Symbol"] or r["Mouse_Gene_ID"], r["Mouse_Official_Name"])
    else:
        mouse = na(str(r["Mapping_Status"] or "not assessed").replace("_", " "))
    cells = [human, mouse] + [safe(cell, r) for cell in
                              (uniprot_cell, go_cell, expression_cell, local_cell)]
    cells.insert(0, safe(input_cell, r))
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
    """Same login as all_sql_loads: asks for the password and database name."""
    user = os.environ.get("MYSQL_USER") or getuser()
    password = getpass(f"MySQL password for {user}: ")
    database = input("Database name: ").strip()
    return mysql.connector.connect(user=user, password=password, database=database,
                                   unix_socket=SOCKET, charset="utf8mb4")


def main():
    connection = connect()
    cursor = connection.cursor(dictionary=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_sql = build_sql(cursor)
    for dataset in DATASETS:
        rows = []
        if summary_sql:
            cursor.execute(summary_sql, (dataset,))
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
