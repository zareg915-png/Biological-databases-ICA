#setup and connection checks
## 1. Check Python and the required libraries
```bash
python3 --version

python3 -c "import sys; print(sys.executable)"

python3 -c "import requests; import xml.etree.ElementTree; print('NCBI libraries OK')"

python3 -c "import mysql.connector; print('MySQL connector OK')"

python3 -c "import pandas, keyring; print('Other libraries OK')"
```
## 2. Open the MySQL console- will prompt for user name and password 

```bash
read -r -p "MySQL username: " BDB_DB_USER
read -r -p "MySQL database name: " BDB_DB_NAME

mysql --no-defaults --socket=/var/run/mysqld/mysqld.sock --user="$BDB_DB_USER" --password --database="$BDB_DB_NAME" --connect-timeout=10
```

## 3. Check the database and permissions in MySQL

```sql
SHOW DATABASES;

SELECT DATABASE();

SHOW TABLES;

SHOW GRANTS;

exit;
```

## 4. Check Python connection to MySQL

```bash
python3 -c '
import getpass
import mysql.connector

username = input("MySQL username: ").strip()
database = input("MySQL database name: ").strip()
db_password = getpass.getpass("MySQL password: ")

try:
    with mysql.connector.connect(
        unix_socket="/var/run/mysqld/mysqld.sock",
        user=username,
        password=db_password,
        database=database,
        connection_timeout=10,
    ) as connection:
        with connection.cursor() as cursor:
            cursor.execute("SELECT DATABASE(), 1")
            print("Connection successful:", cursor.fetchone())

except mysql.connector.Error as error:
    print("Connection failed:", error)
    raise SystemExit(1)'
```
