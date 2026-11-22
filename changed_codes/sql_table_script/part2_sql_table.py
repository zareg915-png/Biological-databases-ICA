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
BASE_DIR = Path(os.environ.get(
    "BDB_PART2_OUTPUT_DIR", str(PROJECT_ROOT / "OP1_data/part2_data"),
)).expanduser()
SQL_OUTPUT_DIR = PROJECT_ROOT / "OP2_SQL/part2_sql"
DATASETS = ("development", "challenge")
TABLE_NAME = "p2_orthologs"
PART1_TABLE_NAME = "p1_annotations"
USERNAME = os.environ.get("MYSQL_USER") or getuser()
MYSQL_SOCKET = "/var/run/mysqld/mysqld.sock"
MYSQL_OPTION_FILE = Path.home() / ".my.cnf"
_FALLBACK_PASSWORD = None

# Exactly the OUTPUT_COLUMNS written by the attached Part 2 script.
INPUT_COLUMNS = (
    "Human_Input_ID", "Human_Gene_ID", "Human_Symbol", "Human_Ensembl_Gene_ID",
    "Mouse_Gene_ID", "Mouse_Symbol", "Mouse_Tax_ID", "Mouse_Ensembl_Gene_ID",
    "Mouse_Official_Name", "Input_Row", "Same_Gene_Note", "Resolution_Status",
    "Resolution_Note", "Mouse_Ortholog_Count", "Mapping_Status", "Mapping_Note",
    "Human_Ensembl_Status", "Mouse_Ensembl_Status",
)
DATABASE_COLUMNS = ("Dataset",) + INPUT_COLUMNS
MAPPED_STATUSES = {"ortholog_found", "multiple_mouse_orthologs"}
FAILURE_STATUSES = {
    "human_lookup_failed", "ortholog_lookup_failed", "invalid_human_gene_id",
    "input_gene_is_not_human", "human_record_not_current",
}
ALL_STATUSES = MAPPED_STATUSES | FAILURE_STATUSES | {
    "input_not_resolved", "no_mouse_ortholog_returned",
}


def read_input_tables(base_dir=BASE_DIR):
    """Validate both files completely before making any database changes."""
    tables = {}
    for dataset in DATASETS:
        path = Path(base_dir) / dataset / "ncbi_ortholog.tsv"
        rows = []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            # Part 2 uses csv.DictWriter: decode its normal CSV quoting.
            # Part 1's QUOTE_NONE reader would corrupt quoted Part 2 notes.
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            if reader.fieldnames != list(INPUT_COLUMNS):
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
                        row[name] = None  # Not assessed; distinct from a successful count of 0.
                    elif not re.fullmatch(pattern, value) or int(value) > 4294967295:
                        raise ValueError(f"Invalid {name}: {where}")
                    else:
                        row[name] = int(value)
                rows.append(row)
        validate_input_groups(rows, dataset)
        tables[dataset] = rows
    return tables


def validate_input_groups(rows, dataset):
    """One outcome per input, or one row for each distinct returned mouse gene."""
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
            if (
                first["Resolution_Status"] != "resolved" or not first["Human_Gene_ID"]
                or None in mouse_ids or len(mouse_ids) != len(group)
                or any(row["Mouse_Tax_ID"] != "10090" for row in group)
                or first["Mouse_Ortholog_Count"] != len(group) or status != expected
            ):
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


def source_counts(rows):
    inputs = {row["Input_Row"] for row in rows}
    mapped = [row for row in rows if row["Mapping_Status"] in MAPPED_STATUSES]
    return (len(inputs), len(rows), len({row["Input_Row"] for row in mapped}),
            len({(row["Human_Gene_ID"], row["Mouse_Gene_ID"]) for row in mapped}))


def print_counts(dataset, counts):
    print(f"{dataset}: {counts[0]} input rows; {counts[1]} stored outcome rows; "
          f"{counts[2]} mapped inputs; {counts[3]} unique human-mouse pairs.")


def build_table_sql(collation=None):
    """Keep the lecturer's table builder, with suitable keys and field types."""
    suffix = " ENGINE=InnoDB DEFAULT CHARACTER SET utf8mb4"
    if collation is not None:
        # Dataset is a text JOIN key. Match the existing Part 1 collation;
        # otherwise MySQL can reject the JOIN with 'Illegal mix of collations'.
        if not re.fullmatch(r"utf8mb4_[A-Za-z0-9_]+", collation):
            raise ValueError("Part 1 Dataset must use an utf8mb4 collation.")
        suffix += f" COLLATE {collation}"
    numeric = {"Input_Row": "INT UNSIGNED NOT NULL",
               "Mouse_Ortholog_Count": "INT UNSIGNED NULL"}
    identifiers = {"Human_Gene_ID", "Mouse_Gene_ID", "Mouse_Tax_ID"}
    fields = ["`Row_ID` INT UNSIGNED AUTO_INCREMENT PRIMARY KEY",
              "`Dataset` VARCHAR(20) NOT NULL"]
    for name in INPUT_COLUMNS:
        datatype = numeric.get(name, "VARCHAR(20) NULL" if name in identifiers else "TEXT")
        fields.append(f"`{name}` {datatype}")
    # This index is deliberately non-unique: an input can have several orthologs.
    fields.extend([
        "INDEX `idx_dataset_input` (`Dataset`, `Input_Row`)",
        "INDEX `idx_human_gene` (`Human_Gene_ID`)",
        "INDEX `idx_mouse_gene` (`Mouse_Gene_ID`)",
    ])
    return (f"CREATE TABLE IF NOT EXISTS `{TABLE_NAME}` (\n  "
            + ",\n  ".join(fields)
            + "\n)" + suffix)


def build_table_data_sql():
    names = ", ".join(f"`{name}`" for name in DATABASE_COLUMNS)
    placeholders = ", ".join(["%s"] * len(DATABASE_COLUMNS))
    # Values are passed separately, preserving quotes, tabs and biological text.
    return f"INSERT INTO `{TABLE_NAME}` ({names}) VALUES ({placeholders})"


def input_values(dataset, rows):
    return [(dataset,) + tuple(row[name] for name in INPUT_COLUMNS) for row in rows]


def check_table(cursor):
    cursor.execute(f"SHOW COLUMNS FROM `{TABLE_NAME}`")
    missing = {"Row_ID", *DATABASE_COLUMNS} - {row[0] for row in cursor.fetchall()}
    if missing:
        raise ValueError(f"Existing {TABLE_NAME} is missing columns: {', '.join(sorted(missing))}")
    cursor.execute(
        "SELECT ENGINE FROM information_schema.TABLES "
        "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s", (TABLE_NAME,),
    )
    engine = cursor.fetchone()
    if not engine or str(engine[0]).upper() != "INNODB":
        raise ValueError(f"{TABLE_NAME} must use InnoDB for rollback.")
    cursor.execute(f"SHOW INDEX FROM `{TABLE_NAME}`")
    indexes = defaultdict(list)
    for row in cursor.fetchall():
        indexes[row[2]].append(row)
    existing = {tuple(row[4] for row in sorted(parts, key=lambda item: item[3]))
                for parts in indexes.values() if all(row[7] is None for row in parts)}
    for name, columns in (
        ("idx_dataset_input", ("Dataset", "Input_Row")),
        ("idx_human_gene", ("Human_Gene_ID",)), ("idx_mouse_gene", ("Mouse_Gene_ID",)),
    ):
        if columns not in existing:
            fields = ", ".join(f"`{column}`" for column in columns)
            cursor.execute(f"ALTER TABLE `{TABLE_NAME}` ADD INDEX `{name}` ({fields})")


def check_input_alignment(cursor, dataset, rows):
    """Match each original input, keeping MYC and its aliases separate."""
    cursor.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS "
        "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s", (PART1_TABLE_NAME,),
    )
    required = {"Dataset", "Input_Row", "Query_Symbol", "Gene_ID", "Resolution_Status",
                "Same_Gene_Note", "Resolution_Note"}
    if required - {row[0] for row in cursor.fetchall()}:
        raise ValueError("Run the updated Part 1 SQL loader first: p1_annotations "
                         "must contain Dataset, Input_Row and its annotation columns.")
    cursor.execute(
        f"SELECT `Input_Row`, `Query_Symbol`, `Gene_ID`, `Resolution_Status`, "
        f"`Same_Gene_Note`, `Resolution_Note` FROM `{PART1_TABLE_NAME}` "
        "WHERE `Dataset` = %s ORDER BY `Input_Row`", (dataset,),
    )
    originals = cursor.fetchall()
    original = {row[0]: row[1:] for row in originals}
    if None in original or len(original) != len(originals):
        raise ValueError("Re-run the updated Part 1 loader to populate distinct Input_Row "
                         "values; adding the column alone does not populate it.")
    if not original or set(original) != {row["Input_Row"] for row in rows}:
        raise ValueError(f"Part 2 {dataset} does not cover the stored Part 1 inputs.")
    for row in rows:
        expected = original[row["Input_Row"]]
        # Older Part 1 loaders stored blank notes as NULL; both mean empty text.
        expected = (expected[0], str(expected[1]) if expected[1] else None,
                    expected[2], expected[3] or "", expected[4] or "")
        observed = tuple(row[name] for name in (
            "Human_Input_ID", "Human_Gene_ID", "Resolution_Status",
            "Same_Gene_Note", "Resolution_Note",
        ))
        if observed != expected:
            raise ValueError(f"Part 2 {dataset}, input {row['Input_Row']} differs from "
                             "stored Part 1. Reload the matching Part 1/Part 2 results.")


def replace_and_verify(connection, cursor, tables):
    """Commit both imports together only after reading every field back."""
    connection.start_transaction()
    try:
        # Check before deleting; a Part 1 update between runs must be detected.
        for dataset, rows in tables.items():
            check_input_alignment(cursor, dataset, rows)
        for dataset, rows in tables.items():
            cursor.execute(f"DELETE FROM `{TABLE_NAME}` WHERE `Dataset` = %s", (dataset,))
            cursor.executemany(build_table_data_sql(), input_values(dataset, rows))
        names = ", ".join(f"`{name}`" for name in DATABASE_COLUMNS)
        for dataset, rows in tables.items():
            cursor.execute(f"SELECT {names} FROM `{TABLE_NAME}` "
                           "WHERE `Dataset` = %s ORDER BY `Row_ID`", (dataset,))
            if cursor.fetchall() != input_values(dataset, rows):
                raise ValueError(f"Imported {dataset} values differ from the Part 2 TSV.")
        connection.commit()
    except BaseException:
        connection.rollback()
        raise


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


def export_sql_dump(database, username, password, output_dir=SQL_OUTPUT_DIR):
    """Export schema and rows; keep the previous dump if the export fails."""
    program = shutil.which("mysqldump") or shutil.which("mariadb-dump")
    temporary_output = None
    try:
        if not program:
            raise RuntimeError("mysqldump/mariadb-dump was not found on this server.")
        output_dir = Path(output_dir).expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / "part2_orthologs.sql"
        with tempfile.TemporaryDirectory(prefix="bdb_part2_dump_") as private_dir:
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
                       "--default-character-set=utf8mb4", f"--user={username}",
                       f"--socket={MYSQL_SOCKET}"]
            for option, setting in (("--column-statistics", "0"), ("--set-gtid-purged", "OFF")):
                if option in help_result.stdout:
                    command.append(f"{option}={setting}")
            command.extend([database, TABLE_NAME])
            with tempfile.NamedTemporaryFile(mode="wb", dir=output_dir, prefix="part2_",
                                             suffix=".tmp", delete=False) as output:
                temporary_output = Path(output.name)
                result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE,
                                        timeout=60, check=False)
            if result.returncode or not temporary_output.stat().st_size:
                raise RuntimeError(f"Dump failed or was empty (exit code {result.returncode}).")
            temporary_output.replace(destination)
    except Exception as error:
        raise RuntimeError("The Part 2 data were committed in MySQL, but the SQL file "
                           f"could not be saved: {error}") from error
    finally:
        if temporary_output is not None:
            temporary_output.unlink(missing_ok=True)
    print(f"[INFO] Saved SQL dump: {destination}")
    return destination


# Use this SELECT for the later HTML integration. LEFT JOIN keeps every Part 1
# input, including unresolved inputs and those with no returned mouse ortholog.
JOIN_SQL = f"""
SELECT p1.Dataset, p1.Input_Row, p1.Query_Symbol, p1.Gene_ID,
       p1.Official_Symbol, p1.Resolution_Status,
       p2.Mouse_Gene_ID, p2.Mouse_Symbol, p2.Mouse_Ensembl_Gene_ID,
       p2.Mouse_Ortholog_Count, p2.Mapping_Status, p2.Mapping_Note
FROM `{PART1_TABLE_NAME}` AS p1
LEFT JOIN `{TABLE_NAME}` AS p2
  ON p1.Dataset = p2.Dataset AND p1.Input_Row = p2.Input_Row
WHERE p1.Dataset = %s
ORDER BY p1.Input_Row, p2.Row_ID
"""


def build_query_database(tables):
    with connect_mysql() as connection:
        with connection.cursor() as cursor:
            choose_database(connection, cursor)
            # Reject stale/mismatched sources before creating or changing Part 2.
            for dataset, rows in tables.items():
                check_input_alignment(cursor, dataset, rows)
            cursor.execute(
                "SELECT COLLATION_NAME FROM information_schema.COLUMNS "
                "WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = %s "
                "AND COLUMN_NAME = 'Dataset'", (PART1_TABLE_NAME,),
            )
            collation = cursor.fetchone()[0]
            cursor.execute(build_table_sql(collation))
            check_table(cursor)
            # DDL commits separately; the following transaction is data-only.
            connection.commit()
            replace_and_verify(connection, cursor, tables)
            print("[INFO] Both Part 2 datasets committed; every imported field verified.")
            database, username = connection.database, connection.user
    export_sql_dump(database, username, _FALLBACK_PASSWORD)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=BASE_DIR,
                        help="Folder containing development/ and challenge/ ortholog TSVs.")
    parser.add_argument("--check-files", action="store_true",
                        help="Validate both TSVs without connecting to MySQL.")
    args = parser.parse_args(argv)
    tables = read_input_tables(args.base_dir.expanduser())
    print("[INFO] Part 2 source files checked:")
    for dataset, rows in tables.items():
        print_counts(dataset, source_counts(rows))
        failures = sorted({row["Input_Row"] for row in rows
                           if row["Mapping_Status"] in FAILURE_STATUSES})
        if failures:
            print(f"[WARN] {dataset}: lookup failures retained for input rows {failures}; "
                  "these are not evidence that no ortholog exists.")
    if args.check_files:
        print("[INFO] File checks passed; Part 1 database alignment is checked during import.")
        return
    build_query_database(tables)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[ERROR] Cancelled.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
