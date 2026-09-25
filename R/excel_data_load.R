library(readxl)
library(tidyverse)
library(purrr)

# Purpose:
# To load the Excel metrics data in a standardised format into [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data_Excel]
# The data are then inserted into [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data]
# The data are finally deduplicated
# Input file path: //Mlcsu-bi-fs/bsolccg/Reports/02_Routine/BBCS Oversight Framework/SQL scripts/Excel Input


#1.  List all Excel files ------------------------------------------------------

cli::cli_h1("Listing Excel files")

excel_files <- list.files(
  path = "//Mlcsu-bi-fs/bsolccg/Reports/02_Routine/BBCS Oversight Framework/SQL scripts/Excel Input",
  pattern = "\\.(xlsx|xlsm|xls)$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Don't read the Data Input Template
excel_files <- excel_files[basename(excel_files) != "Data Input Template.xlsx"]

#2. Establish SQL connection ---------------------------------------------------

sql_connection <- dbConnect(
  odbc::odbc(),
  Driver   = "SQL Server",
  Server   = "MLCSU-BI-SQL",
  Database = "Cluster_BBCS",
  Trusted_Connection = "Yes"
)

#3. Read all Excel files -------------------------------------------------------

cli::cli_h1("Reading Excel files")

all_data <- purrr::map_dfr(
  excel_files,
  metricengineR::read_excel_file,
  sheet_name = "Data Input"
)

cli::cli_alert_success("Process completed.")

head(all_data)

#4. Loading data into SQL ------------------------------------------------------

# Append data to Oversight_Framework_Fact_SQL_Staging_Data_Excel

cli::cli_h1("Loading Excel data into SQL")

DBI::dbWriteTable(
  conn = sql_connection,
  name = DBI::Id(
    schema = "BBCS",
    table = "Oversight_Framework_Fact_SQL_Staging_Data_Excel"
  ),
  value = all_data,
  append = TRUE
)

# Insert data [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data]
dbExecute(sql_connection,
          "  INSERT INTO [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data] (
       [Reference_ID]
      ,[Start_Date]
      ,[End_Date]
      ,[Provider_Code]
      ,[Provider_Site_Code]
      ,[ICB]
      ,[Age]
      ,[Numerator]
      ,[Denominator]
      ,[Ethnicity_Code]
      ,[IMD_Quintile]
      )
	  (
  SELECT [Reference_ID]
      ,[Start_Date]
      ,[End_Date]
      ,[Provider_Code]
      ,[Provider_Site_Code]
      ,[ICB]
      ,[Age]
      ,[Numerator]
      ,[Denominator]
      ,[Ethnicity_Code]
      ,[IMD_Quintile]
	  FROM [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data_Excel]
	  ) "
)

# Remove the Excel staging table from the database
dbExecute(sql_connection,
          "DROP TABLE IF EXISTS [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data_Excel]")


cli::cli_alert_success("Process completed.")

#5. Deduplicate data -----------------------------------------------------------

cli::cli_h1("Data Deduplication")

dbExecute(
  sql_connection,
  "DROP TABLE IF EXISTS #Duplicates

 SELECT PK_ID
           ,ROW_NUMBER() OVER
           (
               PARTITION BY
                   Reference_ID,
                   Start_Date,
                   End_Date,
                   Provider_Code,
                   Provider_Site_Code,
                   ICB,
                   Age,
                   Ethnicity_Code,
                   IMD_Quintile
               ORDER BY PK_ID DESC
           ) AS rn
INTO #Duplicates
    FROM [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data]

  DELETE T1
    FROM [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data] T1
   INNER JOIN #Duplicates T2
      ON T1.PK_ID = T2.PK_ID
   WHERE T2.rn > 1"
)

cli::cli_alert_success("Process completed.")