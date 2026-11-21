#!/usr/bin/env python3



import argparse
from collections import defaultdict
import csv
from getpass import getpass, getuser
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


PROJECT_ROOT = Path(os.environ.get("BDB_PROJECT_ROOT", "~/BDB_ICA")).expanduser()
PART1_BASE_DIR = Path(os.environ.get(
    "BDB_PART1_OUTPUT_DIR", str(PROJECT_ROOT / "OP1_data/part1_data/final_data"),
)).expanduser()
PART2_BASE_DIR = Path(os.environ.get(
    "BDB_PART2_OUTPUT_DIR", str(PROJECT_ROOT / "OP1_data/part2_data"),
)).expanduser()
PART1_SQL_OUTPUT_DIR = PROJECT_ROOT / "OP2_SQL/part1_sql"
PART2_SQL_OUTPUT_DIR = PROJECT_ROOT / "OP2_SQL/part2_sql"
DATASETS = ("development", "challenge")
PART1_TABLE_NAME = "p1_annotations"
PART2_TABLE_NAME = "p2_orthologs"
USERNAME = os.environ.get("MYSQL_USER") or getuser()
MYSQL_SOCKET = "/var/run/mysqld/mysqld.sock"
MYSQL_OPTION_FILE = Path.home() / ".my.cnf"
_FALLBACK_PASSWORD = None

PART1_INPUT_COLUMNS = (
    "Query_Symbol", "Gene_ID", "ENSG_IDs", "Official_Symbol",
    "Official_Name", "Chromosome", "Map_Location", "Other_Aliases",
    "Summary", "Same_Gene_Note", "Resolution_Status", "Resolution_Note",
)
PART1_DATABASE_COLUMNS = ("Dataset", "Input_Row") + PART1_INPUT_COLUMNS
PART2_INPUT_COLUMNS = (
    "Human_Input_ID", "Human_Gene_ID", "Human_Symbol", "Human_Ensembl_Gene_ID",
    "Mouse_Gene_ID", "Mouse_Symbol", "Mouse_Tax_ID", "Mouse_Ensembl_Gene_ID",
    "Mouse_Official_Name", "Input_Row", "Same_Gene_Note", "Resolution_Status",
    "Resolution_Note", "Mouse_Ortholog_Count", "Mapping_Status", "Mapping_Note",
    "Human_Ensembl_Status", "Mouse_Ensembl_Status",
)
PART2_DATABASE_COLUMNS = ("Dataset",) + PART2_INPUT_COLUMNS
MAPPED_STATUSES = {"ortholog_found", "multiple_mouse_orthologs"}
FAILURE_STATUSES = {
    "human_lookup_failed", "ortholog_lookup_failed", "invalid_human_gene_id",
    "input_gene_is_not_human", "human_record_not_current",
}
ALL_STATUSES = MAPPED_STATUSES | FAILURE_STATUSES | {
    "input_not_resolved", "no_mouse_ortholog_returned",
}


def read_part1_tables(base_dir):
    """Part 1 writes plain TSV: keep literal quotation marks intact."""
    tables = {}
    for dataset in DATASETS:
        path = Path(base_dir) / dataset / "ncbi_gene_summary2.tsv"
        rows = []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            if reader.fieldnames != list(PART1_INPUT_COLUMNS):
                raise ValueError(f"Unexpected Part 1 columns in {path}.")
            for row in reader:
                where = f"{path}:{reader.line_num}"
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Wrong number of fields: {where}")
                if not row["Query_Symbol"].strip() or not row["Resolution_Status"].strip():
                    raise ValueError(f"Missing input label or resolution status: {where}")
                gene_id = row["Gene_ID"]
                if gene_id and (len(gene_id) > 20 or not re.fullmatch(r"[1-9][0-9]*", gene_id)):
                    raise ValueError(f"Invalid Gene_ID: {where}")
                if row["Resolution_Status"] == "resolved" and (
                    not gene_id or not row["Official_Symbol"].strip()
                ):
                    raise ValueError(f"Incomplete resolved Part 1 row: {where}")
                # Retain the uploaded Part 1 loader's blank-cell handling.
                rows.append({key: value if value.strip() else None for key, value in row.items()})
        if not rows:
            raise ValueError(f"No Part 1 input rows in {path}.")
        tables[dataset] = rows
    return tables


def read_part2_tables(base_dir):
    """Part 2 uses csv.DictWriter: decode its normal CSV quoting."""
    tables = {}
    for dataset in DATASETS:
        path = Path(base_dir) / dataset / "ncbi_ortholog.tsv"
        rows = []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            if reader.fieldnames != list(PART2_INPUT_COLUMNS):
                raise ValueError(f"Unexpected Part 2 columns in {path}.")
            for row in reader:
                where = f"{path}:{reader.line_num}"
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Wrong number of fields: {where}")
                if not row["Human_Input_ID"].strip() or not row["Resolution_Status"].strip():
                    raise ValueError(f"Missing input label or resolution status: {where}")
                if row["Mapping_Status"] not in ALL_STATUSES:
                    raise ValueError(f"Unrecognised mapping status: {where}")
                for name in ("Human_Gene_ID", "Mouse_Gene_ID", "Mouse_Tax_ID"):
                    value = row[name]
                    if value and (len(value) > 20 or not re.fullmatch(r"[1-9][0-9]*", value)):
                        raise ValueError(f"Invalid {name}: {where}")
                    row[name] = value or None
                for name in ("Input_Row", "Mouse_Ortholog_Count"):
                    value = row[name]
                    pattern = r"[1-9][0-9]*" if name == "Input_Row" else r"[0-9]+"
                    if name == "Mouse_Ortholog_Count" and value == "":
                        row[name] = None  # Unassessed; distinct from a successful count of 0.
                    elif not re.fullmatch(pattern, value) or int(value) > 4294967295:
                        raise ValueError(f"Invalid {name}: {where}")
                    else:
                        row[name] = int(value)
                rows.append(row)
        validate_input_groups(rows, dataset)
        tables[dataset] = rows
    return tables


def validate_input_groups(rows, dataset):
    """Keep one outcome per input, or one row for each distinct mouse ortholog."""
    if not rows:
        raise ValueError(f"No Part 2 outcomes in {dataset}.")
    groups = defaultdict(list)
    for row in rows:
        groups[row["Input_Row"]].append(row)
    if set(groups) != set(range(1, len(groups) + 1)):
        raise ValueError(f"{dataset}: Input_Row must cover 1 to the number of inputs.")
    human_fields = (
        "Human_Input_ID", "Human_Gene_ID", "Human_Symbol", "Human_Ensembl_Gene_ID",
        "Same_Gene_Note", "Resolution_Status", "Resolution_Note", "Human_Ensembl_Status",
        "Mapping_Status", "Mouse_Ortholog_Count",
    )
    for number, group in groups.items():
        first = group[0]
        where = f"{dataset}, input {number} ({first['Human_Input_ID']})"
        if any(any(row[name] != first[name] for name in human_fields) for row in group):
            raise ValueError(f"Inconsistent rows for {where}.")
        status = first["Mapping_Status"]
        if status in MAPPED_STATUSES:
            mouse_ids = {row["Mouse_Gene_ID"] for row in group}
            expected = "ortholog_found" if len(group) == 1 else "multiple_mouse_orthologs"
            if (first["Resolution_Status"] != "resolved" or not first["Human_Gene_ID"]
                    or None in mouse_ids or len(mouse_ids) != len(group)
                    or any(row["Mouse_Tax_ID"] != "10090" for row in group)
                    or first["Mouse_Ortholog_Count"] != len(group) or status != expected):
                raise ValueError(f"Invalid or incomplete mouse mapping for {where}.")
        else:
            mouse_fields = ("Mouse_Gene_ID", "Mouse_Symbol", "Mouse_Tax_ID",
                            "Mouse_Ensembl_Gene_ID", "Mouse_Official_Name")
            expected_count = 0 if status == "no_mouse_ortholog_returned" else None
            if (len(group) != 1 or any(first[name] for name in mouse_fields)
                    or first["Mouse_Ortholog_Count"] != expected_count):
                raise ValueError(f"Invalid empty/failed outcome for {where}.")
            if (status == "input_not_resolved") != (first["Resolution_Status"] != "resolved"):
                raise ValueError(f"Mapping and resolution status disagree for {where}.")
            if status == "no_mouse_ortholog_returned" and not first["Human_Gene_ID"]:
                raise ValueError(f"Successful lookup lacks its human Gene ID: {where}.")


def part1_originals(rows):
    return [(number, row["Query_Symbol"], row["Gene_ID"], row["Resolution_Status"],
             row["Same_Gene_Note"], row["Resolution_Note"])
            for number, row in enumerate(rows, start=1)]


def compare_with_part1(dataset, originals, rows):
    """Use the same identity comparison for the files and stored Part 1 rows."""
    original = {record[0]: record[1:] for record in originals}
    if None in original or len(original) != len(originals):
        raise ValueError("Part 1 needs distinct, populated Input_Row values.")
    if not original or set(original) != {row["Input_Row"] for row in rows}:
        raise ValueError(f"Part 2 {dataset} does not cover the Part 1 inputs.")
    for row in rows:
        source = original[row["Input_Row"]]
        expected = (source[0], str(source[1]) if source[1] else None,
                    source[2], source[3] or "", source[4] or "")
        observed = tuple(row[name] for name in (
            "Human_Input_ID", "Human_Gene_ID", "Resolution_Status",
            "Same_Gene_Note", "Resolution_Note",
        ))
        if expected != observed:
            raise ValueError(f"Part 2 {dataset}, input {row['Input_Row']} differs from Part 1. "
                             "Regenerate the matching results first.")


def source_counts(part, rows):
    if part == 1:
        resolved = [row for row in rows if row["Resolution_Status"] == "resolved"]
        return len(rows), len(resolved), len({row["Gene_ID"] for row in resolved})
    mapped = [row for row in rows if row["Mapping_Status"] in MAPPED_STATUSES]
    return (len({row["Input_Row"] for row in rows}), len(rows),
            len({row["Input_Row"] for row in mapped}),
            len({(row["Human_Gene_ID"], row["Mouse_Gene_ID"]) for row in mapped}))


def print_counts(part, dataset, rows):
    counts = source_counts(part, rows)
    if part == 1:
        print(f"Part 1 {dataset}: {counts[0]} input rows; {counts[1]} resolved inputs; "
              f"{counts[2]} unique resolved genes.")
    else:
        print(f"Part 2 {dataset}: {counts[0]} inputs; {counts[1]} outcome rows; "
              f"{counts[2]} mapped inputs; {counts[3]} unique human-mouse pairs.")


def table_details(part):
    if part == 1:
        return PART1_TABLE_NAME, PART1_DATABASE_COLUMNS
    if part == 2:
        return PART2_TABLE_NAME, PART2_DATABASE_COLUMNS
    raise ValueError("Part must be 1 or 2.")


def build_table_sql(part, collation=None):
    table, _ = table_details(part)
    suffix = " ENGINE=InnoDB DEFAULT CHARACTER SET utf8mb4"
    if collation is not None:
        if not re.fullmatch(r"utf8mb4_[A-Za-z0-9_]+", collation):
            raise ValueError("Dataset must use an utf8mb4 collation.")
        suffix += f" COLLATE {collation}"
    fields = ["`Row_ID` INT UNSIGNED AUTO_INCREMENT PRIMARY KEY",
              "`Dataset` VARCHAR(20) NOT NULL"]
    if part == 1:
        types = {"Query_Symbol": "VARCHAR(255) NOT NULL", "Gene_ID": "VARCHAR(20)",
                 "Official_Symbol": "VARCHAR(255)", "Chromosome": "VARCHAR(20)",
                 "Map_Location": "VARCHAR(100)", "Resolution_Status": "VARCHAR(30)"}
        fields.append("`Input_Row` INT UNSIGNED NOT NULL")
        fields.extend(f"`{name}` {types.get(name, 'TEXT')}" for name in PART1_INPUT_COLUMNS)
        fields.extend(["INDEX `idx_dataset` (`Dataset`)", "INDEX `idx_gene_id` (`Gene_ID`)",
                       "UNIQUE INDEX `idx_dataset_input` (`Dataset`, `Input_Row`)"])
    else:
        types = {"Input_Row": "INT UNSIGNED NOT NULL",
                 "Mouse_Ortholog_Count": "INT UNSIGNED NULL"}
        for name in PART2_INPUT_COLUMNS:
            datatype = types.get(name, "VARCHAR(20)" if name in {
                "Human_Gene_ID", "Mouse_Gene_ID", "Mouse_Tax_ID"} else "TEXT")
            fields.append(f"`{name}` {datatype}")
        # Non-unique: one input can produce several distinct ortholog rows.
        fields.extend(["INDEX `idx_dataset_input` (`Dataset`, `Input_Row`)",
                       "INDEX `idx_human_gene` (`Human_Gene_ID`)",
                       "INDEX `idx_mouse_gene` (`Mouse_Gene_ID`)"])
    return f"CREATE TABLE IF NOT EXISTS `{table}` (\n  " + ",\n  ".join(fields) + "\n)" + suffix


def check_table(cursor, part):
    table, columns = table_details(part)
    cursor.execute(f"SHOW COLUMNS FROM `{table}`")
    found = {row[0] for row in cursor.fetchall()}
    if part == 1 and "Input_Row" not in found:
        # Migrate before checking required fields; the import populates the values.
        cursor.execute(f"ALTER TABLE `{table}` ADD COLUMN `Input_Row` "
                       "INT UNSIGNED NULL AFTER `Dataset`")
        found.add("Input_Row")
        print("[INFO] Added Part 1 Input_Row column.")
    missing = {"Row_ID", *columns} - found
    if missing:
        raise ValueError(f"Existing {table} is missing columns: {', '.join(sorted(missing))}")
    cursor.execute("SELECT ENGINE FROM information_schema.TABLES "
                   "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s", (table,))
    engine = cursor.fetchone()
    if not engine or str(engine[0]).upper() != "INNODB":
        raise ValueError(f"{table} must use InnoDB for rollback.")


def ensure_indexes(cursor, part):
    table, _ = table_details(part)
    cursor.execute(f"SHOW INDEX FROM `{table}`")
    indexes = defaultdict(list)
    for row in cursor.fetchall():
        indexes[row[2]].append(row)
    observed = {}
    for name, rows in indexes.items():
        if all(row[7] is None for row in rows):
            columns = tuple(row[4] for row in sorted(rows, key=lambda item: item[3]))
            observed[name] = (columns, rows[0][1] == 0)
    if part == 1:
        required = [("idx_dataset", ("Dataset",), False),
                    ("idx_gene_id", ("Gene_ID",), False),
                    ("idx_dataset_input", ("Dataset", "Input_Row"), True)]
    else:
        required = [("idx_dataset_input", ("Dataset", "Input_Row"), False),
                    ("idx_human_gene", ("Human_Gene_ID",), False),
                    ("idx_mouse_gene", ("Mouse_Gene_ID",), False)]
        if any(unique and columns == ("Dataset", "Input_Row")
               for columns, unique in observed.values()):
            raise ValueError("Part 2 has a UNIQUE Dataset/Input_Row index, which prevents "
                             "multiple orthologs per input. That index must be non-unique.")
    for name, columns, unique in required:
        if any(existing == columns and (not unique or is_unique)
               for existing, is_unique in observed.values()):
            continue
        if name in indexes:
            raise ValueError(f"Existing index {table}.{name} has an incompatible definition.")
        fields = ", ".join(f"`{column}`" for column in columns)
        kind = "UNIQUE INDEX" if unique else "INDEX"
        cursor.execute(f"ALTER TABLE `{table}` ADD {kind} `{name}` ({fields})")


def check_stored_part1(cursor, tables):
    for dataset, rows in tables.items():
        cursor.execute(f"SELECT `Input_Row`, `Query_Symbol`, `Gene_ID`, `Resolution_Status`, "
                       f"`Same_Gene_Note`, `Resolution_Note` FROM `{PART1_TABLE_NAME}` "
                       "WHERE `Dataset` = %s ORDER BY `Input_Row`", (dataset,))
        compare_with_part1(dataset, cursor.fetchall(), rows)


def input_values(part, dataset, rows):
    if part == 1:
        return [(dataset, number) + tuple(row[name] for name in PART1_INPUT_COLUMNS)
                for number, row in enumerate(rows, start=1)]
    return [(dataset,) + tuple(row[name] for name in PART2_INPUT_COLUMNS) for row in rows]


def replace_and_verify(connection, cursor, part, tables):
    """Replace both datasets for one part in one verified data transaction."""
    table, columns = table_details(part)
    names = ", ".join(f"`{name}`" for name in columns)
    placeholders = ", ".join(["%s"] * len(columns))
    insert = f"INSERT INTO `{table}` ({names}) VALUES ({placeholders})"
    connection.start_transaction()
    try:
        if part == 2:
            check_stored_part1(cursor, tables)
        for dataset, rows in tables.items():
            cursor.execute(f"DELETE FROM `{table}` WHERE `Dataset` = %s", (dataset,))
            cursor.executemany(insert, input_values(part, dataset, rows))
        for dataset, rows in tables.items():
            cursor.execute(f"SELECT {names} FROM `{table}` "
                           "WHERE `Dataset` = %s ORDER BY `Row_ID`", (dataset,))
            if cursor.fetchall() != input_values(part, dataset, rows):
                raise ValueError(f"Imported Part {part} {dataset} values differ from the TSV.")
            if part == 1:
                cursor.execute(
                    "SELECT COUNT(*), COALESCE(SUM(`Resolution_Status` = 'resolved'), 0), "
                    "COUNT(DISTINCT CASE WHEN `Resolution_Status` = 'resolved' "
                    f"THEN `Gene_ID` END) FROM `{table}` WHERE `Dataset` = %s", (dataset,),
                )
                if tuple(int(value) for value in cursor.fetchone()) != source_counts(1, rows):
                    raise ValueError(f"Part 1 {dataset} database counts differ from the TSV.")
        connection.commit()
    except BaseException:
        connection.rollback()
        raise


def load_part(connection, cursor, part, tables):
    collation = None
    if part == 2:
        check_stored_part1(cursor, tables)
        cursor.execute("SELECT COLLATION_NAME FROM information_schema.COLUMNS "
                       "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s "
                       "AND COLUMN_NAME = 'Dataset'", (PART1_TABLE_NAME,))
        record = cursor.fetchone()
        if not record or not record[0]:
            raise ValueError("Cannot determine Part 1's Dataset collation.")
        collation = record[0]  # Match text JOIN keys when creating Part 2.
    cursor.execute(build_table_sql(part, collation))
    check_table(cursor, part)
    ensure_indexes(cursor, part)
    connection.commit()  # Finish DDL before starting a data-only transaction.
    replace_and_verify(connection, cursor, part, tables)
    print(f"[INFO] Part {part} committed; both datasets verified.")


def connect_mysql():
    try:
        from mysql.connector import connect, Error
    except ImportError:
        raise RuntimeError("Install in this Python environment: "
                           "python -m pip install mysql-connector-python") from None
    global _FALLBACK_PASSWORD
    settings = dict(unix_socket=MYSQL_SOCKET, charset="utf8mb4",
                    collation="utf8mb4_unicode_ci", autocommit=False, connection_timeout=15)
    if MYSQL_OPTION_FILE.is_file():
        settings.update(option_files=str(MYSQL_OPTION_FILE), option_groups=["client"])
        if os.environ.get("MYSQL_USER"):
            settings["user"] = USERNAME
    else:
        settings["user"] = USERNAME
    if _FALLBACK_PASSWORD is not None:
        settings["password"] = _FALLBACK_PASSWORD
    try:
        return connect(**settings)
    except Error as error:
        if error.errno != 1045:
            raise
        if os.environ.get("BDB_NON_INTERACTIVE") == "1":
            raise ValueError("MySQL authentication failed. Configure [client] in ~/.my.cnf.") from None
        _FALLBACK_PASSWORD = getpass(f"MySQL password for {settings.get('user', USERNAME)}: ")
        settings["password"] = _FALLBACK_PASSWORD
        return connect(**settings)


def choose_database(connection, cursor):
    cursor.execute("SHOW DATABASES")
    available = [row[0] for row in cursor.fetchall()]
    print("Accessible databases: " + ", ".join(available))
    database = os.environ.get("MYSQL_DATABASE", "").strip() or connection.database
    if not database:
        if os.environ.get("BDB_NON_INTERACTIVE") == "1":
            raise ValueError("Set MYSQL_DATABASE to your personal database name.")
        database = input("Your personal MySQL database name: ").strip()
    if database not in available or database.casefold() in {
        "genbank", "mysql", "sys", "information_schema", "performance_schema",
    }:
        raise ValueError("Choose your personal database from the names shown above.")
    connection.database = database
    print(f"[INFO] Using personal database: {database}")


def quote_option_value(value):
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"').replace(
        "\n", "\\n").replace("\r", "\\r").replace("\t", "\\t") + '"'


def export_sql_dump(database, username, password, part):
    """Save one table's schema/data atomically; credentials stay off the CLI."""
    table, _ = table_details(part)
    output_dir = PART1_SQL_OUTPUT_DIR if part == 1 else PART2_SQL_OUTPUT_DIR
    filename = "part1_annotations.sql" if part == 1 else "part2_orthologs.sql"
    program = shutil.which("mysqldump") or shutil.which("mariadb-dump")
    temporary_output = None
    try:
        if not program:
            raise RuntimeError("mysqldump/mariadb-dump was not found on this server.")
        output_dir = Path(output_dir).expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / filename
        with tempfile.TemporaryDirectory(prefix="bdb_mysql_dump_") as private_dir:
            if password is None and MYSQL_OPTION_FILE.is_file():
                options_path = MYSQL_OPTION_FILE
            else:
                options_path = Path(private_dir) / "client.cnf"
                descriptor = os.open(options_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
                with os.fdopen(descriptor, "w", encoding="utf-8") as options:
                    options.write("[client]\n")
                    for name, value in (("user", username), ("password", password or ""),
                                        ("socket", MYSQL_SOCKET)):
                        options.write(f"{name}={quote_option_value(value)}\n")
            help_result = subprocess.run([program, "--no-defaults", "--help"],
                                         capture_output=True, text=True, timeout=15, check=False)
            command = [program, f"--defaults-file={options_path}", "--single-transaction",
                       "--skip-lock-tables", "--no-tablespaces", "--skip-comments",
                       "--skip-extended-insert", "--default-character-set=utf8mb4",
                       f"--user={username}", f"--socket={MYSQL_SOCKET}"]
            for option, setting in (("--column-statistics", "0"), ("--set-gtid-purged", "OFF")):
                if option in help_result.stdout:
                    command.append(f"{option}={setting}")
            command.extend([database, table])
            with tempfile.NamedTemporaryFile(mode="wb", dir=output_dir, prefix=f"part{part}_",
                                             suffix=".tmp", delete=False) as output:
                temporary_output = Path(output.name)
                result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE,
                                        timeout=60, check=False)
            if result.returncode or not temporary_output.stat().st_size:
                raise RuntimeError(f"Dump failed or was empty (exit code {result.returncode}).")
            temporary_output.replace(destination)
    except Exception as error:
        raise RuntimeError(f"Part {part} is committed in MySQL, but its SQL file "
                           f"could not be saved: {error}") from error
    finally:
        if temporary_output is not None:
            temporary_output.unlink(missing_ok=True)
    print(f"[INFO] Saved SQL dump: {destination}")
    return destination


# Example SELECT for the later HTML step; all original input labels remain.
JOIN_SQL = f"""
SELECT p1.Dataset, p1.Input_Row, p1.Query_Symbol, p1.Gene_ID,
       p1.Official_Symbol, p1.Resolution_Status,
       p2.Mouse_Gene_ID, p2.Mouse_Symbol, p2.Mouse_Ensembl_Gene_ID,
       p2.Mouse_Ortholog_Count, p2.Mapping_Status, p2.Mapping_Note
FROM `{PART1_TABLE_NAME}` AS p1
LEFT JOIN `{PART2_TABLE_NAME}` AS p2
  ON p1.Dataset = p2.Dataset AND p1.Input_Row = p2.Input_Row
WHERE p1.Dataset = %s
ORDER BY p1.Input_Row, p2.Row_ID
"""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part1-dir", type=Path, default=PART1_BASE_DIR,
                        help="Folder containing Part 1 development/ and challenge/ TSVs.")
    parser.add_argument("--part2-dir", type=Path, default=PART2_BASE_DIR,
                        help="Folder containing Part 2 development/ and challenge/ TSVs.")
    parser.add_argument("--check-files", action="store_true",
                        help="Validate all four TSVs and their alignment without MySQL.")
    args = parser.parse_args(argv)
    part1_tables = read_part1_tables(args.part1_dir.expanduser())
    part2_tables = read_part2_tables(args.part2_dir.expanduser())
    for dataset, rows in part2_tables.items():
        compare_with_part1(dataset, part1_originals(part1_tables[dataset]), rows)
    print("[INFO] All four TSVs checked; Part 2 inputs match Part 1.")
    for part, tables in ((1, part1_tables), (2, part2_tables)):
        for dataset, rows in tables.items():
            print_counts(part, dataset, rows)
    for dataset, rows in part2_tables.items():
        failures = sorted({row["Input_Row"] for row in rows
                           if row["Mapping_Status"] in FAILURE_STATUSES})
        if failures:
            print(f"[WARN] {dataset}: retaining lookup failures for inputs {failures}; "
                  "these do not establish that no ortholog exists.")
    if args.check_files:
        print("[INFO] File checks passed. MySQL was not contacted.")
        return

    committed = []
    problems = []
    stage = "MySQL connection/database selection"
    database = username = None
    try:
        with connect_mysql() as connection:
            with connection.cursor() as cursor:
                choose_database(connection, cursor)
                database, username = connection.database, connection.user
                for part, tables in ((1, part1_tables), (2, part2_tables)):
                    stage = f"Part {part}"
                    load_part(connection, cursor, part, tables)
                    committed.append(part)
    except Exception as error:
        problems.append(f"{stage} failed: {error}")
        if committed:
            problems.append("Already committed: " + ", ".join(f"Part {p}" for p in committed))
    for part in committed:
        try:
            export_sql_dump(database, username, _FALLBACK_PASSWORD, part)
        except Exception as error:
            problems.append(str(error))
    if problems:
        raise RuntimeError(" | ".join(problems))
    print("[INFO] Both MySQL tables verified and both SQL dumps saved.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[ERROR] Cancelled.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
