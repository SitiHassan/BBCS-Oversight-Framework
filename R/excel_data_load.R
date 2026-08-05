library(readxl)
library(tidyverse)
library(purrr)

excel_files <- list.files(
  path = "//Mlcsu-bi-fs/bsolccg/Reports/02_Routine/BBCS Oversight Framework/SQL scripts/Excel Input",
  pattern = "\\.(xlsx|xlsm|xls)$",
  full.names = TRUE,
  ignore.case = TRUE
)

excel_files <- excel_files[basename(excel_files) != "Data Input Template.xlsx"]

sql_connection <- dbConnect(
  odbc::odbc(),
  Driver   = "SQL Server",
  Server   = "MLCSU-BI-SQL",
  Database = "Cluster_BBCS",
  Trusted_Connection = "Yes"
)

read_input_file <- function(file_path) {
  
  file_name <- basename(file_path)
  message("Processing file: ", file_name)
  
  tryCatch(
    {
      read_excel(
        path = file_path,
        sheet = "Data Input"
      ) |>
        mutate(Source_File = file_name)
    },
    error = function(e) {
      warning(
        "Could not process ", file_name,
        ": ", conditionMessage(e)
      )
      
      NULL
    }
  )
}

all_data <- map_dfr(
  excel_files,
  read_input_file
)

head(all_data)

DBI::dbWriteTable(
  conn = sql_connection,
  name = DBI::Id(
    schema = "BBCS",
    table = "Oversight_Framework_Fact_SQL_Staging_Data_Excel"
  ),
  value = all_data,
  append = TRUE
)

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

dbExecute(sql_connection,
          "DROP TABLE IF EXISTS [Cluster_BBCS].[BBCS].[Oversight_Framework_Fact_SQL_Staging_Data_Excel]")

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