CREATE TABLE Mart1(
    GeneID varchar(20),  
    Genename varchar(50),  
    Start int,  
    End	int,
    Strand char(1),	
    PRIMARY KEY (GeneID));
    
LOAD DATA LOCAL INFILE './test_mart.txt' INTO TABLE Mart1;
