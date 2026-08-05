get_age_lookup <- function(metadata, age_lookup){
  
  # Inputs:
  #   metadata   - a dataframe containing at least column `age` and `indicator_id`
  #   age_lookup - a dataframe containing `age_group_label` and `age_code`
  # Output:
  #   A dataframe with `indicator_id` and `age_code` columns
  
  lookup <- metadata |>
    mutate(age_group_label = age) |>
    left_join(age_lookup, by = "age_group_label") |>
    select(indicator_id, age_code)
  
  return(lookup)
}

# Function to create combination id column -------------------------------------
create_combination_id <- function(data) {
  
  # Inputs:
  #   data   - a dataframe containing columns `imd_code` and `ethnicity_code`
  # Output:
  #   The same dataframe with an added `combination_id` column
  
  data |>
    mutate(
      combination_id = case_when(
        imd_code != 999 & ethnicity_code != 999 ~ 1L,  # both splits present
        imd_code != 999 & ethnicity_code == 999 ~ 2L,  # IMD only
        imd_code == 999 & ethnicity_code != 999 ~ 3L,  # Ethnicity only
        imd_code == 999 & ethnicity_code == 999 ~ 4L,  # neither (overall)
        TRUE ~ NA_integer_                           # anything unexpected/NA
      )
    )
}


# Function to clean data types -------------------------------------------------
clean_data_types <- function (data){
  
  # Inputs:
  #   data   - a dataframe containing columns to be cleansed
  # Output:
  #   The same dataframe with corrected datatypes
  
  data |>
    mutate(indicator_id    = as.integer(indicator_id),
           start_date      = as.Date(start_date),
           end_date        = as.Date(end_date),
           numerator       = as.numeric(numerator),
           denominator     = as.numeric(denominator),
           indicator_value = as.numeric(indicator_value),
           lower_ci95      = as.numeric(lower_ci95),
           upper_ci95      = as.numeric(upper_ci95),
           imd_code        = as.integer(imd_code),
           aggregation_id  = as.integer(aggregation_id),
           age_group_code  = as.integer(age_group_code),
           sex_code        = as.integer(sex_code),
           ethnicity_code  = as.integer(ethnicity_code),
           creation_date   = as.POSIXct(creation_date),
           value_type_code = as.integer(value_type_code),
           source_code     = as.integer(source_code),
           time_period_type= as.character(time_period_type),
           combination_id  = as.integer(combination_id))
}


# Function to  tidy output -----------------------------------------------------
TIDY_COLS = c("indicator_id","start_date","end_date","numerator","denominator",
              "indicator_value","lower_ci95","upper_ci95","imd_code","aggregation_id",
              "age_group_code","sex_code","ethnicity_code","creation_date",
              "value_type_code","source_code", "time_period_type", "combination_id")

tidy_output <- function(data) {
  
  # Inputs:
  #   data - a dataframe containing calculated results
  #          (e.g. DASR, crude ratio, percentage)
  # Output:
  #   A dataframe with corrected column names and standardised schema
  #   to enable binding results from different calculations
  
  out <- data
  
  # create indicator_value
  if(!"indicator_value" %in% names(out)){
    out$indicator_value <- if("value" %in% names(out)){
      out$value
    } else{
      NA_real_
    }
  }
  
  # create lower confidence interval
  if(!"lower_ci95" %in% names(out)){
    out$lower_ci95 <- if("lowercl" %in% names(out)){
      out$lowercl
    } else{
      NA_real_
    }
  }
  
  # create upper confidence interval
  if(!"upper_ci95" %in% names(out)){
    out$upper_ci95 <- if("uppercl" %in% names(out)){
      out$uppercl
    } else{
      NA_real_
    }
  }
  
  # replace numerator and denominator for DASR outputs
  if("total_count" %in% names(out)){
    out$numerator <- out$total_count
  }
  
  if("total_pop" %in% names(out)){
    out$denominator <- out$total_pop
  }
  
  # check all required output columns exist
  missing_cols <- setdiff(TIDY_COLS, names(out))
  
  if(length(missing_cols) > 0){
    stop(
      "The following required columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out |> 
    dplyr::select(
      dplyr::all_of(
        TIDY_COLS
      )
    )
}

# Function to normalise indicator ids ------------------------------------------
# To normalise indicator ids so that function accepts either
# "All" / "all" / "*" : process all indicators, or
# a vector of IDs (numeric or character), or
# a single comma-separated string like "10, 11, 12"
# especially useful when specifying which ids to extract from SQL using
# get_indicators_from_sql function and which ids to process using run_all function

normalize_indicator_ids <- function(indicator_ids) {
  # Treat NULL or "all"/"*"(any case) as no filter
  if (is.null(indicator_ids)) return(NULL)
  if (is.character(indicator_ids) && length(indicator_ids) == 1) {
    if (tolower(indicator_ids) %in% c("all", "*")) return(NULL)
    # If a single comma-separated string, split it
    if (grepl(",", indicator_ids)) {
      indicator_ids <- trimws(unlist(strsplit(indicator_ids, ",")))
    }
  }
  # Coerce numerics when possible; keep characters that aren't numeric
  suppressWarnings({
    as_num <- suppressWarnings(as.numeric(indicator_ids))
  })
  out <- ifelse(!is.na(as_num), as_num, indicator_ids)
  out <- unique(out)
  out[!is.na(out) & out != ""]
}
