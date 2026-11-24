# Biological Databases ICA 2026: human vs mouse ortholog comparison

Annotates two lists of ten human genes (development and challenge), maps them
to mouse orthologs, compares each human/mouse pair on UniProt GO terms and
subcellular location (Jaccard) and Bgee tissue expression, cross-checks the
NCBI annotation against the local course API, loads all results into MySQL on
bioinfmsc9 and writes one HTML summary table per list from a SQL query.

## Run order

Run from the repo root (`~/BDB_ICA`) on bioinfmsc9. Each step reads the
previous step's output. Every connector processes the development list and
then the challenge list.

| Step | Script | Output |
|------|--------|--------|
| 1. NCBI basic annotation | `changed_codes/data_ret_script/part1/final_script_run_for_both_lists/final_script_for_both_end_res.py` | `OP1_data/part1_data/final_data/{development,challenge}/` |
| 2. NCBI orthologs | `changed_codes/data_ret_script/part2/final_ortho_script.py` | `OP1_data/part2_data/{development,challenge}/` |
| 3. UniProt + Bgee | `changed_codes/data_ret_script/part3/uniprot_bgee_script.py` | `OP1_data/part3_data/{development,challenge}/` |
| 4. Local API (bioinfmsc6 `/summary`) | `changed_codes/data_ret_script/part4/local_api.py` | `OP1_data/part4_data/{development,challenge}/` |
| 5. Load MySQL (Parts 1 to 4) and dump | `changed_codes/sql_table_script/all_sql_loads.py` | 9 tables (`p1_annotations`, `p2_orthologs`, `p3_gene_status`, `p3_uniprot`, `p3_bgee_calls`, `p3_pair_summary`, `p3_tissue_comparison`, `p4_local_api_records`, `p4_local_ncbi_comparisons`); dump `OP2_SQL/all_parts_sql/bdb_ica_all_parts.sql` |
| 6. Summary SQL and HTML | `changed_codes/HTML_script/make_html_tables.py` | `OP3_HTML/development_summary_table.html`, `OP3_HTML/challenge_summary_table.html` |

```bash
export NCBI_EMAIL="you@example.ac.uk"   # requested by NCBI
python3 changed_codes/data_ret_script/part1/final_script_run_for_both_lists/final_script_for_both_end_res.py
python3 changed_codes/data_ret_script/part2/final_ortho_script.py --run both
python3 changed_codes/data_ret_script/part3/uniprot_bgee_script.py --run both
python3 changed_codes/data_ret_script/part4/local_api.py --dataset both
python3 changed_codes/sql_table_script/all_sql_loads.py --check-files   # validate TSVs only
python3 changed_codes/sql_table_script/all_sql_loads.py                 # load 9 tables + write dump
python3 changed_codes/HTML_script/make_html_tables.py                   # both HTML pages
```

To restore the database from the dump: `mysql <your_database> < OP2_SQL/all_parts_sql/bdb_ica_all_parts.sql`.

## Folders

- `changed_codes/data_ret_script/`: the four REST connectors (Parts 1 to 4).
- `changed_codes/sql_table_script/`: MySQL loader (`all_sql_loads.py`).
- `changed_codes/HTML_script/`: summary query and HTML renderer.
- `OP1_data/`: TSV and JSON outputs of Parts 1 to 4.
- `OP2_SQL/all_parts_sql/`: mysqldump of all nine tables.
- `OP3_HTML/`: the two HTML summary tables.
- `set_up/README.md`: environment and MySQL connection checks.

Earlier iterations kept for reference, not part of the final run:
`part1/development/`, `part1/challenge/`, `final_script_for_both_commented.py`,
`part2/first_ortho_script.py`, the Parts 1 to 2 loaders
(`merged_part1_2.py`, `part1_sql_table.py`, `part2_sql_table.py`) and their
dumps in `OP2_SQL/part1_sql/` and `OP2_SQL/part2_sql/`.

## Requirements

Python 3, `requests`, `mysql-connector-python`. MySQL credentials are read
from `~/.my.cnf` or prompted for; none are stored in this repo.

## GenAI use

GenAI use is declared in the header of each script.

