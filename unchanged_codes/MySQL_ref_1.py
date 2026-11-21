#Import MySQL Connector Functionality
import mysql.connector

#Specify the Connection & Connect
db = mysql.connector.connect (
    host="localhost",
    user="stomlins",
    password="########",
    database="stomlins_db"
)

#Issue Select commands on the table & print the results
csr =db.cursor()
csr.execute("SELECT * FROM Mart1 LIMIT 5")
for x in csr:
    print(x)
	


#Inserting data

toinsert= (
        "INSERT INTO Mart1(GeneID,Genename,Start,End,Strand)"
        "VALUES (%s,%s,%s,%s,%s)"
)

data = ('ENSTEST','Test','1','1','1')

csr.execute (toinsert,data)

#if you update the database then call
db.commit()


#check for the data
csr.execute("SELECT * from Mart1 WHERE GeneID='ENSTEST'")

res = csr.fetchall()
for i in res:
    print(i)

#Tidy up at the end
db.close()

