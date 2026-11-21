# #########################################################################
#  Custom script to handle db connectivity                               #
#  Notice how the key is kept private                                    #
#  Simon Tomlinson Nov 2026                                              #
# #########################################################################


import keyring
import pandas as pd

#
#This code assumes we are running on a local machine with a local tunnel
#port 9999 to 3306 on the server
#make this tunnel first before running  the code.

USERNAME="guest101"
DATABASE="guest101"
HOST="127.0.0.1"
PORT="9999"

#Set the password on the system from the Python Console using
#keyring.set_password("DATABASE", "USERNAME", "password")

#Then it can be retrieved later using
#keyring.get_password("DATABASE", "USERNAME")

#Build a simple table to store the data
def build_table_sql(tablename, fieldnames,fieldtypes= None):


    #check the fieldtypes
    print("[INFO] Setting to default field type")
    if fieldtypes is None:
        fieldtypes = ["TEXT"]* len(tablename)

    elif len(fieldtypes) < len(tablename):
        print("[WARNING] not enough fields supplied- using defaults")
        fieldtypes = ["TEXT"] * len(tablename)

    tab = "CREATE TABLE " + tablename + "(\n"

    for i in range(0,len(fieldnames)-1):
        tab = tab + fieldnames[i] + " " + fieldtypes[i] + ",\n"

    tab = tab + fieldnames[len(fieldnames)-1] + " " + fieldtypes[len(fieldnames)-1] + "\n"
    tab =tab+ ");\n"

    return tab
#this build the SQL to insert data
def build_table_data_sql(table_name,table):

    #INSERT INTO table_name (c1, c2, c3, ...)
    #VALUES (v1, v2, v3, ...);

    #build for first part
    sql ="INSERT INTO " + table_name +" ("

    for i in range(0,len(table.columns)-1):
        sql=sql + table.columns[i] +","
    sql = sql + table.columns[len(table.columns)-1] +")\n"

    sql_header =sql

    sql2=""

    for row in range(0,len(table)):
        sql2 = sql2+ sql_header
        sql2 = sql2 + "VALUES("

        crow = table.iloc[row]

        #use the headers to iterate the rows

        for value in range(0,len(table.columns)-1):

            ccol =table.columns[value]
            vval = crow[ccol]
            vval =str(vval)

            #replace internal ' with
            vval = vval.replace("'"," ")
            vval =vval.replace('"'," ")
            vval =vval.replace(";",".")

            vval = "\'" + vval + "\'"
            sql2=sql2 + vval +", "

        ccol = table.columns[len(table.columns)-1]
        vval = crow[ccol]
        vval = str(vval)
        # replace internal ' with
        vval = vval.replace("'", " ")
        vval = vval.replace('"', " ")
        vval = vval.replace(";", ".")

        vval = "\'" + vval + "\'"
        sql2=sql2 + vval + ");\n"

    return sql2




#this build the database
def build_query_database(gene_table,orth_table,string_table):

    #first read in the three tables
    gt=pd.read_csv(gene_table+".tsv", sep='\t',dtype=str)
    ot=pd.read_csv(orth_table+".tsv", sep='\t',dtype=str)
    st=pd.read_csv(string_table+".tsv", sep='\t',dtype=str)

    header_gt = gt.columns
    header_ot = ot.columns
    header_st = st.columns


    from mysql.connector import connect, Error

    try:
        with connect(
            host=HOST,
            user=USERNAME,
            port=PORT,
            database=DATABASE,
            password=keyring.get_password(DATABASE, USERNAME),
        ) as connection:
            print(connection)

            #add rest of connection code here...

            with connection.cursor() as cursor:
                print("[INFO] DROP existing tables")
                cursor.execute("DROP TABLE IF EXISTS "+gene_table)
                cursor.execute("DROP TABLE IF EXISTS "+orth_table)
                cursor.execute("DROP TABLE IF EXISTS "+string_table)


                #upload new table descriptions
                print("[INFO] Uploading table descriptions")
                sq =build_table_sql(gene_table,header_gt)


                cursor.execute(sq)
                sq = build_table_sql(orth_table, header_ot)
                cursor.execute(sq)


                sq = build_table_sql(string_table, header_st)
                cursor.execute(sq)


                print("[INFO] Adding table data")

                sq = build_table_data_sql(gene_table, gt)
                sq = sq +" \n "+ build_table_data_sql(orth_table, ot)
                sq = sq + " \n " + build_table_data_sql(string_table, st)

                print(sq)

                sq2= sq.split(";")

                for val in sq2:
                    if val != "" and val != "\n":
                        cursor.execute(val+";")
                        connection.commit()

                #perform join command
                qury ="SELECT Gene_ID,Official_Symbol,Chromosome,Map_Location,Interactions FROM ncbi_gene_summary2,string_interactions WHERE ncbi_gene_summary2.Official_Symbol=string_interactions.Gene;"
                cursor.execute(qury)

                #download final results
                results = cursor.fetchall()

                df = pd.DataFrame(results, columns=['Gene_ID','Official_Symbol','Chromosome','Map_Location','Interactions'])

                print("[INFO] Finished query, final result, returning....")

                print(df.to_string())

                return(df)

    except Error as e:
        print(e)

build_query_database("ncbi_gene_summary2","ncbi_ortholog","string_interactions")