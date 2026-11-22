# Biological Databases ICA 2026: human vs mouse ortholog comparison

Annotates two lists of ten human genes (development and challenge), maps them
to mouse orthologs, compares each human/mouse pair on UniProt function, GO and
subcellular location (Jaccard) and Bgee tissue expression, loads the results
into MySQL on bioinfmsc9 and writes an HTML summary table per list.

## Run order

Run from the repo root (`~/BDB_ICA`) on bioinfmsc9. Each step reads the
previous step's output.

| Step | Script | Output |
|------|--------|--------|
| 1. NCBI basic annotation | `changed_codes/data_ret_script/part1/final_script_run_for_both_lists/final_script_for_both_end_res.py` | `OP1_data/part1_data/final_data/{development,challenge}/` |
| 2. NCBI orthologs | `changed_codes/data_ret_script/part2/final_ortho_script.py` | `OP1_data/part2_data/{development,challenge}/` |
| 3. UniProt + Bgee | `changed_codes/data_ret_script/part3/uniprot_bgee_script.py` | `OP1_data/part3_data/{development,challenge}/` |
| 4. Local API | _to do_ | |
| 5. Load MySQL (Parts 1+2) | `changed_codes/sql_table_script/merged_part1_2.py` | tables `p1_annotations`, `p2_orthologs`; dumps in `OP2_SQL/` |
| 6. Summary SQL + HTML | _to do_ | |

```bash
export NCBI_EMAIL="you@example.ac.uk"   # requested by NCBI
python3 changed_codes/data_ret_script/part1/final_script_run_for_both_lists/final_script_for_both_end_res.py
python3 changed_codes/data_ret_script/part2/final_ortho_script.py --run both
python3 changed_codes/data_ret_script/part3/uniprot_bgee_script.py --run both
python3 changed_codes/sql_table_script/merged_part1_2.py --check-files   # validate only
python3 changed_codes/sql_table_script/merged_part1_2.py                 # load + dump
```

`part1_sql_table.py` and `part2_sql_table.py` load one part each;
`merged_part1_2.py` does both in one run.

## Folders

- `changed_codes/data_ret_script/`: API connectors (Parts 1 to 3). `part1/development` and `part1/challenge` hold earlier iterations.
- `changed_codes/sql_table_script/`: MySQL loaders.
- `OP1_data/`: TSV and JSON outputs of every part.
- `OP2_SQL/`: mysqldump files.
- `set_up/README.md`: environment and MySQL connection checks.

## Requirements

Python 3, `requests`, `mysql-connector-python`. MySQL credentials are read
from `~/.my.cnf` or prompted for; none are stored in this repo.

## GenAI use

GenAI use is declared in the header of each script.
