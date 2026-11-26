/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `p1_annotations` (
  `Dataset` varchar(20) NOT NULL,
  `Input_Row` bigint DEFAULT NULL,
  `Query_Symbol` longtext,
  `Official_Symbol` longtext,
  `Official_Name` longtext,
  `Gene_ID` bigint DEFAULT NULL,
  `ENSG_IDs` longtext,
  `Map_Location` longtext,
  `Other_Aliases` longtext,
  `Summary` longtext,
  `Resolution_Status` longtext,
  `Resolution_Note` longtext,
  `Same_Gene_Note` longtext
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `p2_orthologs` (
  `Dataset` varchar(20) NOT NULL,
  `Human_Input_ID` longtext,
  `Human_Gene_ID` bigint DEFAULT NULL,
  `Human_Symbol` longtext,
  `Human_Ensembl_Gene_ID` longtext,
  `Mouse_Gene_ID` bigint DEFAULT NULL,
  `Mouse_Symbol` longtext,
  `Mouse_Tax_ID` bigint DEFAULT NULL,
  `Mouse_Ensembl_Gene_ID` longtext,
  `Mouse_Official_Name` longtext,
  `Input_Row` bigint DEFAULT NULL,
  `Same_Gene_Note` longtext,
  `Resolution_Status` longtext,
  `Resolution_Note` longtext,
  `Mouse_Ortholog_Count` bigint DEFAULT NULL,
  `Mapping_Status` longtext,
  `Mapping_Note` longtext,
  `Human_Ensembl_Status` longtext,
  `Mouse_Ensembl_Status` longtext
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `p3_pair_summary` (
  `Dataset` varchar(20) NOT NULL,
  `Input_Row` bigint DEFAULT NULL,
  `Query_Symbol` longtext,
  `Human_Gene_ID` bigint DEFAULT NULL,
  `Human_Symbol` longtext,
  `Mouse_Gene_ID` bigint DEFAULT NULL,
  `Mouse_Symbol` longtext,
  `Resolution_Status` longtext,
  `Mapping_Status` longtext,
  `Human_UniProt_Status` longtext,
  `Human_UniProt_Accessions` longtext,
  `Human_Protein_Names` longtext,
  `Human_GO_Terms` longtext,
  `Human_Function` longtext,
  `Human_Subcellular_Locations` longtext,
  `Human_Location_Evidence` longtext,
  `Human_Bgee_Status` longtext,
  `Mouse_UniProt_Status` longtext,
  `Mouse_UniProt_Accessions` longtext,
  `Mouse_Protein_Names` longtext,
  `Mouse_GO_Terms` longtext,
  `Mouse_Function` longtext,
  `Mouse_Subcellular_Locations` longtext,
  `Mouse_Location_Evidence` longtext,
  `Mouse_Bgee_Status` longtext,
  `GO_MF_Jaccard` double DEFAULT NULL,
  `GO_BP_Jaccard` double DEFAULT NULL,
  `GO_CC_Jaccard` double DEFAULT NULL,
  `Location_Jaccard` double DEFAULT NULL,
  `Expression_Comparison_Status` longtext,
  `Comparable_Tissue_Count` longtext,
  `Both_Present_Count` longtext,
  `Both_Absent_Count` longtext,
  `Shared_Present_Tissues` longtext,
  `Expression_Call_Agreement` longtext,
  `Present_Tissue_Jaccard` longtext,
  `Pair_Note` longtext
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `p4_local_ncbi_comparisons` (
  `Dataset` varchar(20) NOT NULL,
  `Input_Row` bigint DEFAULT NULL,
  `Query_Symbol` longtext,
  `Gene_ID` bigint DEFAULT NULL,
  `Official_Symbol` longtext,
  `NCBI_Map_Location` longtext,
  `Local_Gene_ID` bigint DEFAULT NULL,
  `Local_Official_Symbol` longtext,
  `Local_Map_Location` longtext,
  `Local_Status` longtext,
  `Symbol_Comparison` longtext,
  `Map_Location_Comparison` longtext,
  `Comparison_Status` longtext,
  `Differences` longtext,
  `Local_Note` longtext
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
