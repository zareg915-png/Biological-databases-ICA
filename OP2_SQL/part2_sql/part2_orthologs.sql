
/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;
DROP TABLE IF EXISTS `p2_orthologs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `p2_orthologs` (
  `Row_ID` int unsigned NOT NULL AUTO_INCREMENT,
  `Dataset` varchar(20) NOT NULL,
  `Human_Input_ID` text,
  `Human_Gene_ID` varchar(20) DEFAULT NULL,
  `Human_Symbol` text,
  `Human_Ensembl_Gene_ID` text,
  `Mouse_Gene_ID` varchar(20) DEFAULT NULL,
  `Mouse_Symbol` text,
  `Mouse_Tax_ID` varchar(20) DEFAULT NULL,
  `Mouse_Ensembl_Gene_ID` text,
  `Mouse_Official_Name` text,
  `Input_Row` int unsigned NOT NULL,
  `Same_Gene_Note` text,
  `Resolution_Status` text,
  `Resolution_Note` text,
  `Mouse_Ortholog_Count` int unsigned DEFAULT NULL,
  `Mapping_Status` text,
  `Mapping_Note` text,
  `Human_Ensembl_Status` text,
  `Mouse_Ensembl_Status` text,
  PRIMARY KEY (`Row_ID`),
  KEY `idx_dataset_input` (`Dataset`,`Input_Row`),
  KEY `idx_human_gene` (`Human_Gene_ID`),
  KEY `idx_mouse_gene` (`Mouse_Gene_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=41 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

LOCK TABLES `p2_orthologs` WRITE;
/*!40000 ALTER TABLE `p2_orthologs` DISABLE KEYS */;
INSERT INTO `p2_orthologs` VALUES (21,'development','BRCA1','672','BRCA1','ENSG00000012048','12189','Brca1','10090','ENSMUSG00000017146','breast cancer 1, early onset',1,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (22,'development','BRCA2','675','BRCA2','ENSG00000139618','12190','Brca2','10090','ENSMUSG00000041147','breast cancer 2, early onset',2,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (23,'development','TP53','7157','TP53','ENSG00000141510','22059','Trp53','10090','ENSMUSG00000059552','transformation related protein 53',3,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (24,'development','EGFR','1956','EGFR','ENSG00000146648','13649','Egfr','10090','ENSMUSG00000020122','epidermal growth factor receptor',4,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (25,'development','MYC','4609','MYC','ENSG00000136997','17869','Myc','10090','ENSMUSG00000022346','Myc proto-oncogene, bHLH transcription factor',5,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (26,'development','APOE','348','APOE','ENSG00000130203','11816','Apoe','10090','ENSMUSG00000002985','apolipoprotein E',6,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (27,'development','CFTR','1080','CFTR','ENSG00000001626','12638','Cftr','10090','ENSMUSG00000041301','CF transmembrane conductance regulator',7,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (28,'development','ACE2','59272','ACE2','ENSG00000130234','70008','Ace2','10090','ENSMUSG00000015405','angiotensin converting enzyme 2',8,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (29,'development','TNF','7124','TNF','ENSG00000232810','21926','Tnf','10090','ENSMUSG00000024401','tumor necrosis factor',9,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (30,'development','VEGFA','7422','VEGFA','ENSG00000112715','22339','Vegfa','10090','ENSMUSG00000023951','vascular endothelial growth factor A',10,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (31,'challenge','MYC','4609','MYC','ENSG00000136997','17869','Myc','10090','ENSMUSG00000022346','Myc proto-oncogene, bHLH transcription factor',1,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (32,'challenge','c-Myc','4609','MYC','ENSG00000136997','17869','Myc','10090','ENSMUSG00000022346','Myc proto-oncogene, bHLH transcription factor',2,'Resolves to the same gene as \'MYC\' (Gene ID 4609).','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (33,'challenge','Blah',NULL,'','',NULL,'',NULL,'','',3,'','not found','No current NCBI Gene record matched this input in the selected organism.',NULL,'input_not_resolved','Not queried: Part 1 did not mark this input as resolved; see Resolution_Status and Resolution_Note.','not_requested','not_requested');
INSERT INTO `p2_orthologs` VALUES (34,'challenge','ENSG00000136997','4609','MYC','ENSG00000136997','17869','Myc','10090','ENSMUSG00000022346','Myc proto-oncogene, bHLH transcription factor',4,'Resolves to the same gene as \'MYC\' (Gene ID 4609).','resolved','Verified the exact Ensembl gene cross-reference.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (35,'challenge','AADACL2','344752','AADACL2','ENSG00000197953','639634','Aadacl2','10090','ENSMUSG00000091376','arylacetamide deacetylase like 2',5,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (36,'challenge','Plet1','349633','PLET1','ENSG00000188771','76509','Plet1','10090','ENSMUSG00000032068','placenta expressed transcript 1',6,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (37,'challenge','HYDIN','54768','HYDIN','ENSG00000157423','244653','Hydin','10090','ENSMUSG00000059854','HYDIN, axonemal central pair apparatus protein',7,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
INSERT INTO `p2_orthologs` VALUES (38,'challenge','ENSG00000261846',NULL,'','',NULL,'',NULL,'','',8,'','not found','No current NCBI Gene record matched this input in the selected organism.',NULL,'input_not_resolved','Not queried: Part 1 did not mark this input as resolved; see Resolution_Status and Resolution_Note.','not_requested','not_requested');
INSERT INTO `p2_orthologs` VALUES (39,'challenge','XIST','7503','XIST','ENSG00000229807',NULL,'',NULL,'','',9,'','resolved','Selected the unique exact official-symbol match.',0,'no_mouse_ortholog_returned','NCBI successfully returned no mouse ortholog for this human Gene ID.','from_part1','not_requested');
INSERT INTO `p2_orthologs` VALUES (40,'challenge','MALAT1','378938','MALAT1','ENSG00000251562','72289','Malat1','10090','ENSMUSG00000092341','metastasis associated lung adenocarcinoma transcript 1 (non-coding RNA)',10,'','resolved','Resolved from one current symbol/alias search match.',1,'ortholog_found','Retained every mouse ortholog returned by NCBI; one row per human-mouse pair.','from_part1','from_ncbi_datasets');
/*!40000 ALTER TABLE `p2_orthologs` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

