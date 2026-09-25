# Establish sql connection
sql_connection <-
  dbConnect(
    odbc(),
    Driver = "SQL Server",
    Server = "MLCSU-BI-SQL",
    Database = "Cluster_BBCS",
    Trusted_Connection = "True"
  )

# Fingertips IDs
ids <- c(93725, 93726,92600, 94063, 30311)

# Map Fingertips ids to NOF reference ids
ref_table <- tibble(
  fingertips_id = c(93725, 93726, 92600, 94063, 30311),
  reference_id = c(8.03, 8.03, 8.04, 8.05, 8.06)
)

# Fetch data
fingertips_data <- metricengineR::get_fingertips_indicators(ids)


# Get metadata to populate key columns
metadata <- dbGetQuery(
  sql_connection, "
SELECT e.indicator_id
, a.reference_id
, b.age_code
, c.sex_code
, 999 as imd_code
, 999 as ethnicity_code
, a.precalculated
, d.value_type_code
, a.source_code
FROM [Cluster_BBCS].[BBCS].[Oversight_Framework_Reference_Metadata] a
LEFT JOIN [Cluster_BBCS].[BBCS].[Metric_Engine_Reference_Age_Group] b
ON a.age = b.age_group_label
LEFT JOIN [Cluster_BBCS].[BBCS].[Metric_Engine_Reference_Sex] c
ON a.sex = c.sex
LEFT JOIN [Cluster_BBCS].[BBCS].[Metric_Engine_Reference_Value_Type] d
ON a.value_type = d.value_type
LEFT JOIN [Cluster_BBCS].[BBCS].[Oversight_Framework_Reference_Indicator_List] e
ON a.reference_id = e.Reference_ID
"
)

# Update metadata to add Fingertips ids
metadata <- metadata %>%
  left_join(
    ref_table, by = "reference_id"
  )

final_dt <- fingertips_data %>%
  filter(
    Area.Type == "ICBs",
    Area.Name %in% c(
      "NHS Birmingham and Solihull Integrated Care Board - QHL",
      "NHS Black Country Integrated Care Board - QUA"
    )
  ) %>%
  left_join(
    metadata,
    by = c("Indicator.ID" = "fingertips_id")
  ) %>%
  group_by(
    indicator_id,
    reference_id,
    precalculated,
    Area.Name,
    Time.period,
    age_code,
    sex_code,
    imd_code,
    ethnicity_code,
    value_type_code,
    source_code
  ) %>%
  summarise(
    numerator = if_else(
      all(is.na(Count)), NA_real_, sum(Count, na.rm = TRUE)), # To prevent NULL numerators from being assigned as zeros
    denominator = sum(Denominator, na.rm = TRUE),
    indicator_value = if_else(
      first(precalculated) == "Yes",
      first(Value),
      NA_real_
    ),
    lower_ci95 = if_else(
      first(precalculated) == "Yes",
      first(Lower.CI.95.0.limit),
      NA_real_
    ),
    upper_ci95 = if_else(
      first(precalculated) == "Yes",
      first(Upper.CI.95.0.limit),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    start_date = dmy(paste0("01-04-", substr(Time.period, 1, 4))),
    end_date = dmy(paste0("31-03-20", substr(Time.period, 6, 7))),
    aggregation_id = case_when(
      Area.Name == "NHS Birmingham and Solihull Integrated Care Board - QHL" ~ 151L,
      Area.Name == "NHS Black Country Integrated Care Board - QUA" ~ 163L,
      TRUE ~ NA_integer_
    ),
    age_group_code = age_code,
    creation_date = Sys.time()
  ) %>%
  select(
    indicator_id,
    start_date,
    end_date,
    numerator,
    denominator,
    indicator_value,
    lower_ci95,
    upper_ci95,
    imd_code,
    aggregation_id,
    age_group_code,
    sex_code,
    ethnicity_code,
    creation_date,
    value_type_code,
    source_code
  ) %>%
  mutate(
    indicator_id     = as.integer(indicator_id),
    numerator        = as.numeric(numerator),
    denominator      = as.numeric(denominator),
    indicator_value  = as.numeric(indicator_value),
    lower_ci95       = as.numeric(lower_ci95),
    upper_ci95       = as.numeric(upper_ci95),
    age_group_code   = as.integer(age_group_code),
    sex_code         = as.integer(sex_code),
    value_type_code  = as.integer(value_type_code),
    source_code      = as.integer(source_code)
  )

# Insert into API staging table
dbWriteTable(
  sql_connection,
  name = DBI::Id(
    schema = "BBCS",
    table = "Oversight_Framework_Fact_API_Data"
  ),
  value = final_dt,
  append = TRUE,
  row.names = FALSE
)

DBI::dbDisconnect(conn)

