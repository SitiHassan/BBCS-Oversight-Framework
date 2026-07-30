library(tidyverse)
library(lubridate)
library(PHEindicatormethods)

# calculate_dsr2 ---------------------------------------------------------------
# This function ONLY allows the multiplier argument to be passed as a scalar value

calculate_dsr2 <-
  function (data, x, n, stdpop = NULL, type = "full", confidence = 0.95,
            multiplier = 1e+05, independent_events = TRUE, eventfreq = NULL,
            ageband = NULL) {
    if (missing(data) | missing(x) | missing(n) | missing(stdpop)) {
      stop("function calculate_dsr requires at least 4 arguments: data, x, n, stdpop")
    }
    if (!is.data.frame(data)) {
      stop("data must be a data frame object")
    }
    if (!deparse(substitute(x)) %in% colnames(data)) {
      stop("x is not a field name from data")
    }
    if (!deparse(substitute(n)) %in% colnames(data)) {
      stop("n is not a field name from data")
    }
    if (!deparse(substitute(stdpop)) %in% colnames(data)) {
      stop("stdpop is not a field name from data")
    }
    data <- data %>% rename(x = {
      {
        x
      }
    }, n = {
      {
        n
      }
    }, stdpop = {
      {
        stdpop
      }
    })
    if (!is.numeric(data$x)) {
      stop("field x must be numeric")
    }
    else if (!is.numeric(data$n)) {
      stop("field n must be numeric")
    }
    else if (!is.numeric(data$stdpop)) {
      stop("field stdpop must be numeric")
    }
    else if (anyNA(data$n)) {
      stop("field n cannot have missing values")
    }
    else if (anyNA(data$stdpop)) {
      stop("field stdpop cannot have missing values")
    }
    else if (any(pull(data, x) < 0, na.rm = TRUE)) {
      stop("numerators must all be greater than or equal to zero")
    }
    else if (any(pull(data, n) <= 0)) {
      stop("denominators must all be greater than zero")
    }
    else if (any(pull(data, stdpop) < 0)) {
      stop("stdpop must all be greater than or equal to zero")
    }
    else if (!(type %in% c("value", "lower", "upper",
                           "standard", "full"))) {
      stop("type must be one of value, lower, upper, standard or full")
    }
    else if (!is.numeric(confidence)) {
      stop("confidence must be numeric")
    }
    else if (length(confidence) > 2) {
      stop("a maximum of two confidence levels can be provided")
    }
    else if (length(confidence) == 2) {
      if (!(confidence[1] == 0.95 & confidence[2] == 0.998)) {
        stop("two confidence levels can only be produced if they are specified as 0.95 and 0.998")
      }
    }
    else if ((confidence < 0.9) | (confidence > 1 & confidence <
                                   90) | (confidence > 100)) {
      stop("confidence level must be between 90 and 100 or between 0.9 and 1")
    }
    else if (!is.numeric(multiplier)) {
      stop("multiplier must be numeric")
    }
    else if (multiplier <= 0) {
      stop("multiplier must be greater than 0")
    }
    else if (!rlang::is_bool(independent_events)) {
      stop("independent_events must be TRUE or FALSE")
    }
    if (!independent_events) {
      if (missing(eventfreq)) {
        stop(paste0("function calculate_dsr requires an eventfreq column ",
                    "to be specified when independent_events is FALSE"))
      }
      else if (!deparse(substitute(eventfreq)) %in% colnames(data)) {
        stop("eventfreq is not a field name from data")
      }
      else if (!is.numeric(data[[deparse(substitute(eventfreq))]])) {
        stop("eventfreq field must be numeric")
      }
      else if (anyNA(data[[deparse(substitute(eventfreq))]])) {
        stop("eventfreq field must not have any missing values")
      }
      if (missing(ageband)) {
        stop(paste0("function calculate_dsr requires an ageband column ",
                    "to be specified when independent_events is FALSE"))
      }
      else if (!deparse(substitute(ageband)) %in% colnames(data)) {
        stop("ageband is not a field name from data")
      }
      else if (anyNA(data[[deparse(substitute(ageband))]])) {
        stop("ageband field must not have any missing values")
      }
    }
    if (independent_events) {
      dsrs <- dsr_inner2(data = data, x = x, n = n, stdpop = stdpop,
                         type = type, confidence = confidence, multiplier = multiplier)
    }
    else {
      data <- data %>% rename(eventfreq = {
        {
          eventfreq
        }
      }, ageband = {
        {
          ageband
        }
      }) %>% group_by(eventfreq, .add = TRUE)
      grps <- group_vars(data)[!group_vars(data) %in% "eventfreq"]
      check_groups <- filter(summarise(group_by(data, pick(all_of(c(grps,
                                                                    "ageband")))), num_n = n_distinct(.data$n),
                                       num_stdpop = n_distinct(.data$stdpop), .groups = "drop"),
                             .data$num_n > 1 | .data$num_stdpop > 1)
      if (nrow(check_groups) > 0) {
        stop(paste0("There are rows with the same grouping variables and ageband",
                    " but with different populations (n) or standard populations",
                    "(stdpop)"))
      }
      freq_var <- data %>% dsr_inner2(x = x, n = n, stdpop = stdpop,
                                      type = type, confidence = confidence, multiplier = multiplier,
                                      rtn_nonindependent_vardsr = TRUE) %>% mutate(freqvars = .data$vardsr *
                                                                                     .data$eventfreq^2) %>% group_by(pick(all_of(grps))) %>%
        summarise(custom_vardsr = sum(.data$freqvars), .groups = "drop")
      event_data <- data %>% mutate(events = .data$eventfreq *
                                      .data$x) %>% group_by(pick(all_of(c(grps, "ageband",
                                                                          "n", "stdpop")))) %>% summarise(x = sum(.data$events,
                                                                                                                  na.rm = TRUE), .groups = "drop")
      dsrs <- event_data %>% left_join(freq_var, by = grps) %>%
        group_by(pick(all_of(grps))) %>% dsr_inner2(x = x,
                                                    n = n, stdpop = stdpop, type = type, confidence = confidence,
                                                    multiplier = multiplier, use_nonindependent_vardsr = TRUE)
    }
    return(dsrs)
  }



# calculate_dsr3 ---------------------------------------------------------------
# This function allows the multiplier argument to be passed as a column name

calculate_dsr3 <-
  function (data, x, n, stdpop = NULL, type = "full", confidence = 0.95,
            multiplier = 1e+05, independent_events = TRUE, eventfreq = NULL,
            ageband = NULL) {
    
    # ---- basic argument checks ----
    if (missing(data) | missing(x) | missing(n) | missing(stdpop)) {
      stop("function calculate_dsr requires at least 4 arguments: data, x, n, stdpop")
    }
    if (!is.data.frame(data)) stop("data must be a data frame object")
    if (!deparse(substitute(x)) %in% colnames(data)) stop("x is not a field name from data")
    if (!deparse(substitute(n)) %in% colnames(data)) stop("n is not a field name from data")
    if (!deparse(substitute(stdpop)) %in% colnames(data)) stop("stdpop is not a field name from data")
    
    # ---- standardise core columns ----
    data <- dplyr::rename(
      data,
      x      = {{ x }},
      n      = {{ n }},
      stdpop = {{ stdpop }}
    )
    
    if (!is.numeric(data$x))      stop("field x must be numeric")
    else if (!is.numeric(data$n)) stop("field n must be numeric")
    else if (!is.numeric(data$stdpop)) stop("field stdpop must be numeric")
    else if (anyNA(data$n))       stop("field n cannot have missing values")
    else if (anyNA(data$stdpop))  stop("field stdpop cannot have missing values")
    else if (any(dplyr::pull(data, x) < 0, na.rm = TRUE)) stop("numerators must all be >= 0")
    else if (any(dplyr::pull(data, n) <= 0))              stop("denominators must all be > 0")
    else if (any(dplyr::pull(data, stdpop) < 0))          stop("stdpop must all be >= 0")
    else if (!(type %in% c("value", "lower", "upper", "standard", "full")))
      stop("type must be one of value, lower, upper, standard or full")
    else if (!is.numeric(confidence))
      stop("confidence must be numeric")
    else if (length(confidence) > 2)
      stop("a maximum of two confidence levels can be provided")
    else if (length(confidence) == 2) {
      if (!(confidence[1] == 0.95 & confidence[2] == 0.998)) {
        stop("two confidence levels can only be produced if they are specified as 0.95 and 0.998")
      }
    } else if ((confidence < 0.9) | (confidence > 1 & confidence < 90) | (confidence > 100)) {
      stop("confidence level must be between 90 and 100 or between 0.9 and 1")
    }
    
    # ---- allow multiplier to be a column OR a scalar ----
    mult_is_col <- deparse(substitute(multiplier)) %in% colnames(data)
    
    if (mult_is_col) {
      data <- dplyr::rename(data, multiplier = {{ multiplier }})
      if (!is.numeric(data$multiplier)) stop("multiplier column must be numeric")
      if (anyNA(data$multiplier))       stop("multiplier column cannot have missing values")
      if (any(data$multiplier <= 0))    stop("multiplier column values must be > 0")
    } else {
      # scalar path -> inject as a column so downstream is uniform
      if (!is.numeric(multiplier)) stop("multiplier must be numeric")
      if (length(multiplier) != 1) stop("multiplier must be length 1 when not a column")
      if (multiplier <= 0)         stop("multiplier must be greater than 0")
      data <- dplyr::mutate(data, multiplier = multiplier)
    }
    
    # ---- handle independent / non-independent events ----
    if (!rlang::is_bool(independent_events)) {
      stop("independent_events must be TRUE or FALSE")
    }
    
    if (independent_events) {
      dsrs <- dsr_inner2(
        data = data, x = x, n = n, stdpop = stdpop,
        type = type, confidence = confidence,
        rtn_nonindependent_vardsr = FALSE,
        use_nonindependent_vardsr = FALSE
      )
    } else {
      # validate eventfreq and ageband
      if (missing(eventfreq)) {
        stop("function calculate_dsr requires an eventfreq column when independent_events is FALSE")
      } else if (!deparse(substitute(eventfreq)) %in% colnames(data)) {
        stop("eventfreq is not a field name from data")
      } else if (!is.numeric(data[[deparse(substitute(eventfreq))]])) {
        stop("eventfreq field must be numeric")
      } else if (anyNA(data[[deparse(substitute(eventfreq))]])) {
        stop("eventfreq field must not have any missing values")
      }
      
      if (missing(ageband)) {
        stop("function calculate_dsr requires an ageband column when independent_events is FALSE")
      } else if (!deparse(substitute(ageband)) %in% colnames(data)) {
        stop("ageband is not a field name from data")
      } else if (anyNA(data[[deparse(substitute(ageband))]])) {
        stop("ageband field must not have any missing values")
      }
      
      data <- data %>%
        dplyr::rename(
          eventfreq = {{ eventfreq }},
          ageband   = {{ ageband }}
        ) %>%
        dplyr::group_by(eventfreq, .add = TRUE)
      
      grps <- dplyr::group_vars(data)[!dplyr::group_vars(data) %in% "eventfreq"]
      
      check_groups <- dplyr::group_by(data, dplyr::pick(dplyr::all_of(c(grps, "ageband")))) %>%
        dplyr::summarise(
          num_n      = dplyr::n_distinct(.data$n),
          num_stdpop = dplyr::n_distinct(.data$stdpop),
          .groups = "drop"
        ) %>%
        dplyr::filter(.data$num_n > 1 | .data$num_stdpop > 1)
      
      if (nrow(check_groups) > 0) {
        stop(paste0("There are rows with the same grouping variables and ageband",
                    " but with different populations (n) or standard populations (stdpop)"))
      }
      
      # vardsr under non-independent events (multiplier not required here)
      freq_var <- data %>%
        dsr_inner2(
          x = x, n = n, stdpop = stdpop, type = type,
          confidence = confidence,
          rtn_nonindependent_vardsr = TRUE,
          use_nonindependent_vardsr = FALSE
        ) %>%
        dplyr::mutate(freqvars = .data$vardsr * .data$eventfreq^2) %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(grps))) %>%
        dplyr::summarise(custom_vardsr = sum(.data$freqvars), .groups = "drop")
      
      # collapse events and retain multiplier per group
      event_data <- data %>%
        dplyr::mutate(events = .data$eventfreq * .data$x) %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(c(grps, "ageband", "n", "stdpop")))) %>%
        dplyr::summarise(
          x = sum(.data$events, na.rm = TRUE),
          multiplier = dplyr::first(.data$multiplier),
          .groups = "drop"
        )
      
      dsrs <- event_data %>%
        dplyr::left_join(freq_var, by = grps) %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(grps))) %>%
        dsr_inner2(
          x = x, n = n, stdpop = stdpop, type = type,
          confidence = confidence,
          rtn_nonindependent_vardsr = FALSE,
          use_nonindependent_vardsr = TRUE
        )
    }
    
    return(dsrs)
  }

# dsr_inner2 --------------------------------------------------------------------
dsr_inner2 <-
  function (data, x, n, stdpop, type, confidence,
            rtn_nonindependent_vardsr = FALSE,
            use_nonindependent_vardsr = FALSE) {
    
    if (isTRUE(rtn_nonindependent_vardsr) &&
        ("custom_vardsr" %in% names(data) || isTRUE(use_nonindependent_vardsr))) {
      stop("cannot get nonindependent vardsr and use nonindependent vardsr in the same execution")
    }
    
    # normalise confidence inputs (allow 95 and 0.95 forms)
    confidence[confidence >= 90] <- confidence[confidence >= 90] / 100
    conf1 <- confidence[1]
    conf2 <- confidence[2]
    
    # ensure multiplier column present & valid
    if (!"multiplier" %in% names(data)) {
      stop("internal error: multiplier column not found. Ensure calculate_dsr2() injected/renamed it.")
    }
    if (!is.numeric(data$multiplier)) stop("multiplier column must be numeric")
    if (anyNA(data$multiplier))       stop("multiplier column cannot have missing values")
    if (any(data$multiplier <= 0))    stop("multiplier column values must be > 0")
    
    # multiplier must be constant within groups
    grps <- dplyr::group_vars(data)
    if (length(grps) > 0) {
      mult_check <- dplyr::summarise(data, nuniq_mult = dplyr::n_distinct(.data$multiplier), .groups = "drop")
      if (any(mult_check$nuniq_mult > 1)) stop("Within-group multiplier must be constant.")
    }
    
    # choose method
    if (!isTRUE(use_nonindependent_vardsr)) {
      method <- "Dobson"
      data <- dplyr::mutate(data, custom_vardsr = NA_real_)
    } else {
      method <- "Dobson, with confidence adjusted for non-independent events"
    }
    
    dsrs <- data %>%
      dplyr::mutate(
        wt_rate = PHEindicatormethods:::na.zero(.data$x) * .data$stdpop / .data$n,
        sq_rate = PHEindicatormethods:::na.zero(.data$x) * (.data$stdpop / (.data$n))^2
      ) %>%
      dplyr::summarise(
        total_count = sum(.data$x, na.rm = TRUE),
        total_pop   = sum(.data$n),
        base_value  = sum(.data$wt_rate) / sum(.data$stdpop),
        vardsr = if (isTRUE(use_nonindependent_vardsr)) {
          unique(.data$custom_vardsr)
        } else {
          1 / sum(.data$stdpop)^2 * sum(.data$sq_rate)
        },
        # vardsr = dplyr::case_when(
        #   isTRUE(use_nonindependent_vardsr) ~ unique(.data$custom_vardsr),
        #   TRUE ~ 1 / sum(.data$stdpop)^2 * sum(.data$sq_rate)
        # ),
        multiplier = dplyr::first(.data$multiplier),
        .groups = "keep"
      )
    
    if (!isTRUE(rtn_nonindependent_vardsr)) {
      dsrs <- dsrs %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          value    = .data$base_value * .data$multiplier,
          lowercl  = .data$value + sqrt(.data$vardsr/.data$total_count) *
            (PHEindicatormethods:::byars_lower(.data$total_count, conf1) - .data$total_count) *
            .data$multiplier,
          uppercl  = .data$value + sqrt(.data$vardsr/.data$total_count) *
            (PHEindicatormethods:::byars_upper(.data$total_count, conf1) - .data$total_count) *
            .data$multiplier,
          lower99_8cl = .data$value + sqrt(.data$vardsr/.data$total_count) *
            (PHEindicatormethods:::byars_lower(.data$total_count, 0.998) - .data$total_count) *
            .data$multiplier,
          upper99_8cl = .data$value + sqrt(.data$vardsr/.data$total_count) *
            (PHEindicatormethods:::byars_upper(.data$total_count, 0.998) - .data$total_count) *
            .data$multiplier,
          confidence = paste0(confidence * 100, "%", collapse = ", "),
          statistic  = paste("dsr per", format(.data$multiplier, scientific = FALSE)),
          method     = method
        )
      
      # rename/drops to mirror PHE functions' dual-CI behavior
      if (!is.na(conf2)) {
        names(dsrs)[names(dsrs) == "lowercl"] <- "lower95_0cl"
        names(dsrs)[names(dsrs) == "uppercl"] <- "upper95_0cl"
      } else {
        dsrs <- dplyr::select(dsrs, !c("lower99_8cl", "upper99_8cl"))
      }
    }
    
    if (isTRUE(rtn_nonindependent_vardsr)) {
      dsrs <- dplyr::select(dsrs, dplyr::group_cols(), "vardsr")
    } else if (type == "lower") {
      dsrs <- dplyr::select(dsrs, !c("total_count","total_pop","value",
                                     dplyr::starts_with("upper"),
                                     "vardsr","confidence","statistic","method"))
    } else if (type == "upper") {
      dsrs <- dplyr::select(dsrs, !c("total_count","total_pop","value",
                                     dplyr::starts_with("lower"),
                                     "vardsr","confidence","statistic","method"))
    } else if (type == "value") {
      dsrs <- dplyr::select(dsrs, !c("total_count","total_pop",
                                     dplyr::starts_with("lower"), dplyr::starts_with("upper"),
                                     "vardsr","confidence","statistic","method"))
    } else if (type == "standard") {
      dsrs <- dplyr::select(dsrs, !c("vardsr","confidence","statistic","method"))
    } else if (type == "full") {
      dsrs <- dplyr::select(dsrs, !c("vardsr"))
    }
    
    return(dsrs)
  }

stat_mode <- function(x) {
  tx <- table(x, useNA = "no")              # count frequency of each unique value
  if (!length(tx)) return(NA_character_)    # handle empty input gracefully
  names(tx)[which.max(tx)]                  # return the name (value) with max frequency
}


# Function to check which indicators requiring pooled data ---------------------
# Only data requiring Census population will need to be aggregated to create pooled-years

require_pooled_data <- function(metadata) {
  # Inputs:
  #   metadata   - a dataframe containing at least column `denominator_source` and `indicator_id`
  # Output:
  #   A vector of indicator ids requiring pooled data
  ids <- metadata |>
    filter(str_detect(tolower(denominator_source), "census")) |>
    distinct(indicator_id) |>
    pull(indicator_id)

  return(ids)
}

# Global keys to pool data by
POOL_KEYS <- c(
  "indicator_id","imd_code","aggregation_id",
  "age_group_code","sex_code","ethnicity_code",
  "creation_date","value_type_code","source_code","combination_id"
)

# Function to generate year series ---------------------------------------------
generate_year_series <- function(min_year, max_year, span_years = 3L) {
  stopifnot(is.numeric(min_year), is.numeric(max_year), is.numeric(span_years))
  min_year     <- as.integer(min_year)
  max_year     <- as.integer(max_year)
  window_years <- as.integer(span_years)
  if (window_years < 1L) stop("span_years must be >= 1")

  gap <- window_years - 1L
  if ((max_year - min_year) < gap) {
    return(data.frame(from = integer(0), to = integer(0), k = integer(0)))
  }

  starts <- seq.int(min_year, max_year - gap)
  data.frame(
    from = starts,
    to   = starts + gap,
    k    = window_years
  )
}

# Example
# generate_year_series(2014, 2024, 3)
#
# generate_year_series(2014, 2024, 5)

# Function to create pooled data -----------------------------------------------
# Inputs:
#   df  - data frame containing at least POOL_KEYS, `start_date`, `numerator`,
#         `denominator`, `period_type` ("Calendar"/"Financial"), and `population_type`
#   ks  - integer vector of window sizes to pool over (e.g., c(3, 5))
#   time_period_3yrs (optional) - data frame with columns `from`, `to` defining
#         explicit 3-year windows; if NULL, windows are auto-generated from data
#   time_period_5yrs (optional) - data frame with columns `from`, `to` defining
#         explicit 5-year windows; if NULL, windows are auto-generated from data
#
# Output:
#   A data frame containing only the 3- and/or 5-year pooled rows (no yearly rows),
#   where `numerator` and `denominator` are summed over each window, dates are set
#   from `from`/`to` according to `period_type`, and `time_period_type` is
#   "3 year pooled" or "5 year pooled".
#
# Notes/Assumptions:
#   - Only rows with `population_type == "Census"` are considered.
#   - One yearly record per POOL_KEYS ? year is constructed internally before pooling.
#   - If custom window tables are provided, they must have integer `from`/`to` years.


create_pooled_with_ranges <- function(df, ks = c(3,5),
                                      POOL_KEYS,
                                      time_period_3yrs = NULL,
                                      time_period_5yrs = NULL) {

  # Select only indicators which use Census as denominator source
  df <- df |> filter(population_type == "Census")

  # Group data by POOL KEYS and year
  yearly <- df |>
    mutate(period_year = lubridate::year(as.Date(start_date))) |>
    filter(period_type %in% c("Calendar","Financial"),
           !is.na(period_year)) |>
    group_by(across(all_of(POOL_KEYS)), period_type, period_year, ) |>
    summarise(
      numerator   = sum(numerator,   na.rm = TRUE),
      denominator = sum(denominator, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(yearly) == 0) return(slice(df, 0)) # empty in, empty out

  # Build range windows (use data span if not provided)
  miny <- min(yearly$period_year, na.rm = TRUE)
  maxy <- max(yearly$period_year, na.rm = TRUE)

  ranges_list <- list()
  if (3 %in% ks) {
    rng3 <- if (is.null(time_period_3yrs)) {
      # auto-generate from data
      tibble::tibble(from = seq.int(miny, maxy - 2), to = from + 2, k = 3L)
    } else {
      time_period_3yrs |>  transmute(from = as.integer(from), to = as.integer(to), k = 3L)
    }
    ranges_list <- c(ranges_list, list(rng3))
  }
  if (5 %in% ks) {
    rng5 <- if (is.null(time_period_5yrs)) {
      tibble(from = seq.int(miny, maxy - 4), to = from + 4, k = 5L)
    } else {
      time_period_5yrs |>  transmute(from = as.integer(from), to = as.integer(to), k = 5L)
    }
    ranges_list <- c(ranges_list, list(rng5))
  }
  ranges <- bind_rows(ranges_list)

  # If the ranges sit outside the data years, drop those that cant match
  ranges <- ranges |>  filter(from >= miny, to <= maxy)

  if (nrow(ranges) == 0) {
    return(df |>  slice(0))
  }

  # Non-equi join each yearly row to all windows that cover its year
  pooled <- yearly |>
    inner_join(ranges, by = join_by(period_year >= from, period_year <= to)) |>
    group_by(across(all_of(POOL_KEYS)), period_type, from, to, k) |>
    summarise(
      numerator   = sum(numerator,   na.rm = TRUE),
      denominator = sum(denominator, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      # Map window start/end to actual dates by period_type
      start_date = if_else(
        period_type == "Calendar",
        lubridate::make_date(from, 1, 1),
        lubridate::make_date(from, 4, 1)
      ),
      end_date = if_else(
        period_type == "Calendar",
        lubridate::make_date(to, 12, 31),
        lubridate::make_date(to + 1, 3, 31)
      ),
      time_period_type = paste0(k, " year pooled"),
      indicator_value = NA_real_, lower_ci95 = NA_real_, upper_ci95 = NA_real_
    ) |>
    transmute(
      indicator_id, start_date, end_date,
      numerator, denominator,
      indicator_value, lower_ci95, upper_ci95,
      imd_code, aggregation_id, age_group_code, sex_code, ethnicity_code,
      creation_date, value_type_code, source_code, time_period_type, combination_id
    )

  pooled
}


# Function to get age lookup ---------------------------------------------------

age_lookup <- dbGetQuery(
  conn,
  "SELECT * FROM [Cluster_BBCS].[BBCS].[Metric_Engine_Reference_Age_Group]"
)


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

  out <- data |>
    mutate(
      indicator_value = value,
      lower_ci95      = lowercl,
      upper_ci95      = uppercl
    )

  # If DASR output: replace numerator/denominator
  if ("total_count" %in% names(out)) {
    out <- out |> mutate(numerator = .data$total_count)
  }
  if ("total_pop" %in% names(out)) {
    out <- out |> mutate(denominator = .data$total_pop)
  }

  out <- out |> select(all_of(TIDY_COLS))

  return(out)
}

# Function to identify time period type ----------------------------------------

get_duration_label <- function(start_date, end_date) {

  # Inputs:
  #   start_date - a Date vector
  #   end_date   - a Date vector
  # Output:
  #   A character vector of duration labels, aligned with input length
  
  # helper to check if period matches monthly pattern
  is_n_months <- function(n){
    anniv <- start_date %m+% months(n)
    (end_date == anniv) | (end_date == anniv - days(1))
  }

  # helper to check if period matches exactly n years
  is_n_years <- function(n) {
    anniv <- start_date %m+% years(n)
    (end_date == anniv) | (end_date == anniv - days(1))
  }

  case_when(
    is_n_months(1) ~ "Monthly",
    is_n_years(1) ~ "1 year",
    is_n_years(3) ~ "3 year pooled",
    is_n_years(5) ~ "5 year pooled",
    TRUE ~ NA_character_
  )
}

# Function to exclude incomplete data ------------------------------------------
# Automatically exclude incomplete data cycles (e.g., financial, calendar, or other)
# based on the typical end date (mode month-day) and typical duration between start and end dates.
exclude_incomplete_data <- function(
    df,
    period_type_col, # column name: year type (e.g., "Calendar", "Financial", "Other")
    start_date_col, # column name: start date
    end_date_col, # column name: end date
    time_period_col = NULL, # optional column name for time period type (e.g., "1 year", "3 year pooled")
    indicator_col = NULL,     # optional column name for indicator (e.g., "indicator_id")
    current_date = Sys.Date(), # current date reference for detecting incomplete cycles
    enforce_duration = TRUE, # whether to filter by typical duration
    tolerance_days = 5L  # allowed deviation (days) from typical duration
) {
  # Capture column names for tidy evaluation
  period_type <- rlang::ensym(period_type_col)
  start_date <- rlang::ensym(start_date_col)
  end_date <- rlang::ensym(end_date_col)
  time_period <- if (!is.null(time_period_col)) rlang::ensym(time_period_col) else NULL
  indicator   <- if (!is.null(indicator_col)) rlang::ensym(indicator_col) else NULL

  # Standardize year type and time period text for grouping
  df2 <- df |>
    mutate(
      .row_id  = dplyr::row_number(), # row id for tracking
      yt_lower = stringr::str_to_lower(as.character(!!period_type)),
      tp_lower = if (!is.null(time_period)) stringr::str_to_lower(as.character(!!time_period)) else NA_character_
    )

  # Determine which columns to group by (year type + optional time period)
  group_keys <- c("yt_lower", if (!is.null(time_period)) "tp_lower" else NULL)

  meta <- df2 |>
    filter(!is.na(!!end_date), !!end_date <= current_date) |> # keep only past data
    mutate(
      mmdd = format(!!end_date, "%m-%d"),  # extract month-day part of end date
      duration_days = if_else(!is.na(!!start_date), as.numeric(!!end_date - !!start_date) + 1, NA_real_)
    ) |>
    group_by(across(all_of(group_keys))) |>
    summarise(
      end_mmdd = stat_mode(mmdd), # most common end month-day (e.g., 03-31)
      typical_duration = stats::median(duration_days, na.rm = TRUE), # median number of days per cycle
      .groups = "drop"
    )

  # Compute filtered result
  result <- df2 |>
    left_join(meta, by = group_keys) |>  # Join metadata back to the original data and filter incomplete periods
    mutate(
      end_mmdd = coalesce(end_mmdd, "12-31"), # Default to Dec 31 if no mode end date found
      latest_month_end = floor_date(current_date, "month") - days(1),
      # Construct this year's expected cycle end date
      cycle_end_this_year = as.Date(sprintf("%d-%s", lubridate::year(current_date), end_mmdd)),
      # If this years cycle hasnt ended yet, use last years end date
      latest_cycle_end = case_when(
        tp_lower == "monthly" ~ latest_month_end,
        cycle_end_this_year <= current_date ~ cycle_end_this_year,
        TRUE ~ as.Date(sprintf("%d-%s", lubridate::year(current_date) - 1L, end_mmdd))
      ),
      # latest_cycle_end = if_else(
      #   cycle_end_this_year <= current_date, cycle_end_this_year,
      #   as.Date(sprintf("%d-%s", lubridate::year(current_date) - 1L, end_mmdd))
      # ),
      # Calculate actual duration for each record
      duration_days = if_else(!is.na(!!start_date) & !is.na(!!end_date),
                                     as.numeric(!!end_date - !!start_date) + 1, NA_real_),
      # Check if duration is within the expected tolerance range
      duration_ok = if_else(
        enforce_duration & !is.na(typical_duration) & !is.na(duration_days),
        abs(duration_days - typical_duration) <= tolerance_days,
        TRUE
      )
    ) |>
    # Keep only rows that end on or before the latest completed cycle and have valid duration
    filter(!is.na(!!end_date), !!end_date <= latest_cycle_end, duration_ok)

  # Determine excluded rows (anti-join by row id)
  excluded <- df2 |> dplyr::anti_join(result |> dplyr::select(.row_id), by = ".row_id")

  # Print summary of exclusions
  total_excluded <- nrow(excluded)
  total_rows     <- nrow(df)
  if (total_excluded == 0) {
    message("No rows were excluded (0/", total_rows, ").")
  } else {
    message("Excluded ", total_excluded, " row(s) out of ", total_rows, ".")
    
    # Show excluded rows and unique start and end dates
    print(head(excluded))
    print(paste0("unique start dates: (", unique(excluded$start_date)))
    print(paste0("unique end dates: (", unique(excluded$end_date)))
    print(unique(excluded$indicator_id))
    
    if (!is.null(indicator)) {
      by_indicator <- excluded |>
        dplyr::mutate(`Indicator` = as.character(!!indicator)) |>
        dplyr::count(`Indicator`, name = "Excluded_Rows") |>
        dplyr::arrange(dplyr::desc(Excluded_Rows))
      print(by_indicator)
    }
  }

  # Return filtered data without helper columns
  result |>
    dplyr::select(-yt_lower, -tp_lower, -cycle_end_this_year) |>
    dplyr::select(-.row_id)
}



# Function to calculate difference ---------------------------------------

calc_difference <- function(df_in){
  
  df_in |> 
    filter(value_type_code %in% c(10, 11)) |> # 10 = Percentage point difference, 11 = Crude difference
    mutate(
      value = numerator - denominator,
      lowercl = NA_real_,
      uppercl = NA_real_
    )
}

# Function to calculate percentage change ---------------------------------------

calc_percentage_change <- function(df_in){
  
  df_in |> 
    filter(value_type_code == 9) |> 
    mutate(
      value = if_else(
        denominator == 0 | is.na(denominator),
        NA_real_, # ignore 0/NA denominator as sometimes there are no plan/previous/baseline values
        (numerator - denominator) / denominator * value_multiplier 
      ),
      lowercl = NA_real_,
      uppercl = NA_real_
    )
}

# Function to calculate percentage ---------------------------------------------

calc_percentage <- function(df_in){

  # Inputs:
  #   df_in - dataframe eligible for processing, containing `numerator`, `denominator`,
  #           `value_multiplier`, and `value_type_code` (2 = percentage)
  # Output:
  #   The same rows (where value_type_code == 2) with percentage and 95% CI
  #   columns added by PHEindicatormethods::phe_proportion:
  #   `value`, `lowercl`, `uppercl`, `method`

  df_perc <- df_in %>%
    filter(value_type_code == 2)
  
  df_proportion <- df_perc |> 
    filter(!is.na(denominator),
           denominator != 0,
           numerator <= denominator)
  
  multiplier_value <- df_proportion$value_multiplier[1]
  
  df_proportion <- df_proportion |> 
    phe_proportion(
      x = numerator,
      n = denominator,
      confidence = 0.95,
      type = "standard",
      multiplier = multiplier_value
    )
  
  # Achievement against target, growth or increase, or repeat events
  df_other_percentage <- df_perc |> 
    filter(is.na(denominator) |
             denominator == 0 |
             numerator > denominator) |> 
    mutate(
      value = if_else(
        !is.na(denominator) & denominator != 0,
        (numerator / denominator) * value_multiplier,
        NA_real_
      ),
      lowercl = NA_real_,
      uppercl = NA_real_,
    )
  
  output <- bind_rows(
    df_proportion,
    df_other_percentage
  )
  
  return(output)

}

# Function to calculate count --------------------------------------------------
calc_count <- function(df_in){
  
  output <- df_in |>
    filter(value_type_code == 1) |>  
    mutate(
      value = numerator,
      lowercl = NA_real_,
      uppercl = NA_real_,
      method = NA_character_
    )
  
  return(output)
  
}

calc_count_poisson <- function(df_in){
  
  output <- df_in |>
    filter(value_type_code == 1) |>
    mutate(
      indicator_value = numerator,
      lowercl = qchisq(0.025, 2 * numerator) / 2,
      uppercl = qchisq(0.975, 2 * (numerator + 1)) / 2,
      method = "Poisson exact"
    )
  
  return(output)
  
}

# Function to calculate crude or ratio -----------------------------------------

calc_crude_or_ratio <- function(df_in){

  # Inputs:
  #   df_in - dataframe eligible for processing, containing `numerator`, `denominator`,
  #           `value_multiplier`, and `value_type_code` (3 = crude rate, 7 = ratio)
  # Output:
  #   The same rows (where value_type_code is either 3 or 7) with crude rates or ratios calculated
  #   and 95% CI columns added by PHEindicatormethods::phe_rate:
  #   `value`, `lowercl`, `uppercl`, `method`

  df_rate <- df_in |>
    filter(value_type_code %in% c(3, 7)) |>
    mutate(
      num_sign = sign(coalesce(numerator, 0)),
      num_abs  = abs(coalesce(numerator, 0))
    )
  
  output <-  df_rate |>
    phe_rate(
      x = num_abs,
      n = denominator,
      confidence = 0.95,
      type = "standard",
      multiplier = df_rate$value_multiplier
    ) |>
    mutate(numerator = num_abs * num_sign,
           value = value * num_sign,
           lower_tmp = pmin(lowercl, uppercl),
           upper_tmp = pmax(lowercl, uppercl),
           lowercl   = if_else(num_sign < 0, -upper_tmp, lower_tmp),
           uppercl   = if_else(num_sign < 0, -lower_tmp, upper_tmp)) |>
    select(-lower_tmp, -upper_tmp)

  
  return(output)
}

# Function to calculate dasr ---------------------------------------------------

calc_dasr <- function(df_in, metadata, age_lookup){

  # Inputs:
  #   df_in      - dataframe containing `numerator`, `denominator`, `value_multiplier`,
  #                and `value_type_code` (4 = DASR)
  #   metadata   - dataframe containing `indicator_id` and `age` columns,
  #                used to derive the final age code
  #   age_lookup - dataframe containing age lookups (e.g., `age_code`,
  #                `age_group_label`)
  #
  # Output:
  #   A dataframe filtered to rows where `value_type_code == 4`, with DASR
  #   calculated and 95% confidence interval columns added using
  #   `calculate_dsr3` adapted from pheindicatormethods (utils.R). The output also
  #   includes the correct `age_group_code` derived from `metadata` and
  #   `age_lookup`.

  esp2013_lookup <- PHEindicatormethods::esp2013 |>
    as_tibble() |>
    rename(stdpop = value) |>
    mutate(age_group_code = as.integer(c(1:18, 18))) |>
    group_by(age_group_code) |>
    summarise(stdpop = sum(stdpop), .groups = "drop")

  dasr_keys <- c(
    "indicator_id","start_date","end_date","imd_code","aggregation_id", "sex_code",
    "ethnicity_code","creation_date", "value_type_code","source_code","time_period_type", "combination_id"
  )

  df_calc <- df_in |>
    filter(value_type_code == 4) |>
    left_join(esp2013_lookup, by = "age_group_code")

  output <- df_calc |>
    group_by(across(all_of(dasr_keys))) |>
    calculate_dsr3(
      x = numerator,
      n = denominator,
      stdpop = stdpop,
      type = "standard",
      multiplier = value_multiplier # this function accept either value_multiplier or a scalar value e.g., 100000
    ) |>
    ungroup() |>
    left_join(get_age_lookup(metadata = metadata, age_lookup = age_lookup), by = "indicator_id") |>
    mutate(age_group_code = age_code)


  return(output)

}

# Function to calculate SII ----------------------------------------------------
calculate_sii <- function(data, group_cols = c("indicator_id", "start_date", "end_date",
                                               "aggregation_id", "age_group_code", "sex_code",
                                               "ethnicity_code", "source_code"),
                          quintile_col, numerator_col, denominator_col){
  
  #1. Check inputs
  
  # Missing function arguments
  if(missing(data) || missing(quintile_col) || missing(numerator_col) || missing(denominator_col)){
    stop("function calculate_sii() requires the arguments: `data`,  `quintile_col`, `numerator_col`, `denominator_col`")
  }
  
  # Capture the supplied column names
  group_names <- group_cols # already a character vector
  quintile_name <- rlang::as_name(rlang::ensym(quintile_col))
  numerator_name <- rlang::as_name(rlang::ensym(numerator_col))
  denominator_name <- rlang::as_name(rlang::ensym(denominator_col))
  
  required_columns <- c(group_names, quintile_name, numerator_name, denominator_name)
  
  missing_columns <- setdiff(required_columns, names(data))
  
  if(length(missing_columns) > 0){
    stop("The following supplied columns do not exist in `data`: ", 
         paste(missing_columns, collapse = ", "))
  }
  
  # Extract the supplied columns needed to calculate SII
  quintiles <- suppressWarnings(
    as.integer(as.character(data[[quintile_name]])) # Convert IMD quintiles to integer if in character
    ) 
  numerators <- data[[numerator_name]]
  denominators <- data[[denominator_name]]
  
  # Check for missing numerators or denominators
  if(anyNA(numerators) || anyNA(denominators)){
    stop("Numerator and denominator cannot contain missing values.")
  }
  
  # Check for non-numeric numerators and denominators
  if(!is.numeric(numerators) || !is.numeric(denominators)){
    stop("Numerators and denominators must be numeric columns.")
  }
  
  # Check for negative denominators
  if(any(denominators <= 0)){
    stop("All denominators must be greater than zero.")
  }
  
  # Check for negative numerators
  if(any(numerators < 0)){
    stop("Numerators cannot be negative.")
  }
  
  # Check for invalid numerators being greater than denominators
  if(any(numerators > denominators)){
    stop("The numerator cannot be greater than the denominator.")
  }
  
  # Check quintiles within every grouping level
  quintile_check <- data |> 
    dplyr::mutate(
      .quintile = quintiles
    ) |> 
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_names))
    ) |> 
    dplyr::summarise(
      row_count = dplyr::n(),
      distinct_quintiles = dplyr::n_distinct(.quintile),
      has_missing_quintile = anyNA(.quintile),
      duplicated_quintile = anyDuplicated(.quintile) > 0,
      correct_quintile_set = setequal(.quintile, 1:5),
      quintiles_found = paste(sort(unique(.quintile)), collapse = ", "),
      .groups = "drop"
    ) |> 
    dplyr::mutate(
      valid_quintiles = row_count == 5 & 
        distinct_quintiles == 5 &
        !has_missing_quintile &
        !duplicated_quintile &
        correct_quintile_set
        
    )
  
  invalid_groups <- quintile_check |> 
    dplyr::filter(!valid_quintiles)
  
  if(nrow(invalid_groups) > 0){
    print(invalid_groups)
    stop(
      nrow(invalid_groups),
      " grouping level(s) do not contain exactly one row ",
      "for each IMD quintile 1, 2, 3, 4 and 5. ",
      "See the printed table for details."
    )
  }
  
  #2. Create variables needed to fit a linear model
  model_data <- data |> 
    dplyr::mutate(
      deprivation_position = dplyr::recode(
        quintiles,
        `1` = 0.9, # Most deprived
        `2` = 0.7,
        `3` = 0.5,
        `4` = 0.3,
        `5` = 0.1 # Least deprived
      ),
      proportion = numerators / denominators,
      variance = proportion * (1 - proportion) / denominators
    )
  
  if (any(model_data$variance == 0)){
    stop(
      "At least one quintile has variance equal to zero. ",
      "This occurs when the proportion is exactly 0 or 1."
    )
  }
  
  if(any(!is.finite(model_data$variance)) || any(model_data$variance <0)){
    stop("One or more variance values are invalid.")
  }
  
  #3. Calculate inverse variance weight
  model_data <- model_data |> 
    dplyr::mutate(
      inverse_variance_weight = 1 / variance
    )
  
  #4. Create one nested dataset per grouping level
  model_data <- model_data |> 
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_names))
    ) |>
    tidyr::nest()
    
  #5. Fit one linear model per grouping level
  model_data <- model_data |> 
    dplyr::mutate(
      model = purrr::map(
        .data$data,
        ~ stats::lm(
          proportion ~ deprivation_position,
          data = .x,
          weights = .x$inverse_variance_weight
          )
        )
      )
  
  #6. Predict the deprivation extremes for each model
  model_data <- model_data |> 
    dplyr::mutate(
      predictions = purrr::map(
        model, 
        ~ stats::predict(
          .x,
          newdata = tibble::tibble(
            deprivation_position = c(0,1)
            )
          )
        )
      )
  
  
  #7. Extract slope and intercept
  model_data <- model_data |> 
    dplyr::mutate(
      slope = purrr:map_dbl(
        model,
        ~ unname(
          stats::coef(.x)[["deprivation_position"]]
        )
      ),
      intercept =  purrr:map_dbl(
        model,
        ~ unname(
          stats::coef(.x)[["(Intercept)"]]
        )
      )
    )
  
  #8. Create final output
  model_data <- model_data |> 
    dplyr::mutate(
    sii_signed_percentage_points = slope * 100,
    sii_absolute_percentage_points = abs(slope * 100),
    predicted_least_deprived_percent = purrr::map_dbl(
      predictions,
      ~ unname(.x[1]) * 100
    ),
    predicted_most_deprived_percent = purrr::map_dbl(
      predictions,
      ~ unname(.x[2]) * 100
    )
  )
  
  #9. Return final table
  model_data |> 
    dplyr::select(
      dplyr::all_of(group_names),
      intercept,
      slope,
      sii_signed_percentage_points,
      sii_absolute_percentage_points,
      predicted_least_deprived_percent,
      predicted_most_deprived_percent
    ) |> 
    dplyr::ungroup()

}


# Function to calculate different value types ----------------------------------

# Purpose:
#   Calculates and standardises indicator values (Percentage, Crude/Ratio, DASR)
#   and returns a unified schema ready for row-binding across methods.
#
# Inputs:
#   data         - dataframe with (at minimum) the following columns:
#                  indicator_id, start_date, end_date, numerator, denominator,
#                  indicator_value (may be NA), value_type_code, imd_code,
#                  aggregation_id, sex_code, ethnicity_code, creation_date,
#                  source_code
#   metadata     - OPTIONAL dataframe used to enrich `data`; must contain:
#                  <metadata_key>, age, period_type, value_multiplier, status_code
#   metadata_key - character scalar giving the join key name in `metadata`
#                  (default: "indicator_id")
#
# Output:
#   A dataframe in the standardised schema `TIDY_COLS`, containing:
#     - Rows not requiring calculation (kept as-is, standardised)
#     - Rows that required calculation but were ineligible (standardised)
#     - Newly calculated results for:
#         * Percentage (via calc_percentage)
#         * Crude rate / Ratio (via calc_crude_or_ratio)
#         * DASR (via calc_dasr)
#   All outputs have consistent columns (e.g., indicator_value, lower_ci95,
#   upper_ci95) and cleaned datatypes.
#
# Eligibility rules for calculation:
#   - indicator_value is NA
#   - denominator is present and > 0
#   - value_multiplier is present (non-NA)
#   - status_code == 1
#   - value_type_code is present (non-NA)

calculate_values <- function(data, metadata, metadata_key = "indicator_id"){

  # -------- Validate inputs / bring in value_multiplier ---------------------------
  message("\u25B6 Cleaning data types and validating inputs...")
  df <- data |>
    mutate(time_period_type = get_duration_label(as.Date(start_date), as.Date(end_date))) |>
    create_combination_id() |>
    clean_data_types()

  # metadata is REQUIRED
  if (is.null(metadata)) {
    stop("`metadata` is required and must not be NULL")
  }

  # Basic checks
  if (!metadata_key %in% names(metadata)) {
    stop(sprintf("\u274C metadata must contain the key column '%s'", metadata_key))
  }
  if (!"value_multiplier" %in% names(metadata)) {
    stop("\u274C metadata must contain a 'value_multiplier' column")
  }
  if (!"status_code" %in% names(metadata)) {
    stop("\u274C metadata must contain a 'status_code' column")
  }

  # Select required columns from metadata
  metadata_cols = c(metadata_key, "age", "period_type", "value_multiplier", "status_code", "population_type", "precalculated")
  metadata2 <- metadata |> select(all_of(metadata_cols))

  # Bring metadata onto staging
  message("\u25B6 Bringing metadata onto staging table...")
  df <- df |>
      left_join(metadata2, by = setNames(metadata_key, metadata_key))

  # Exclude incomplete data
  # Comment this out as we don't want to exclude incomplete data
  # message("\u25B6 Excluding incomplete data..")
  # df <- exclude_incomplete_data(df = df, period_type_col = "period_type",
  #                                   start_date_col = "start_date", end_date_col = "end_date",
  #                                   time_period_col = "time_period_type",
  #                                   current_date = Sys.Date(), enforce_duration = TRUE, tolerance_days = 5L)


  # Create 3 and 5 year pooled data
  message("\u25B6 Creating 3 and 5 year pooled data...")
  # pooled_df <- create_pooled_with_ranges(df, ks = c(3,5), POOL_KEYS = POOL_KEYS) |>
  #   left_join(metadata2, by = setNames(metadata_key, metadata_key))
  pooled_df <- create_pooled_with_ranges(df, ks = c(3,5), POOL_KEYS = POOL_KEYS)

  # Combine pooled data with yearly data
  message("\u25B6 Combining yearly and pooled data...")
  all_df <- bind_rows(pooled_df, df) |>
    select(
      -ends_with(".x"),
      -ends_with(".y")
    )
    

  # --- Partition rows -------------------------------------------------------
  
  message("\u25B6 Determining rows requiring calculations...")
  # needs_calc <- is.na(all_df$indicator_value) & (!is.na(all_df$denominator) & all_df$denominator != 0) # this will bypass count data as denominator will always be zero/null
  needs_calc <- !is.na(all_df$precalculated) &
    all_df$precalculated == "No"
  
  # Only process where status_code == 1 or 2 AND value_multiplier is present
  eligible_for_processing <- needs_calc &
    !is.na(all_df$value_multiplier) &
    all_df$status_code %in% c(1,2) &
    !is.na(all_df$value_type_code)

  # To be calculated
  df_calc <- all_df[eligible_for_processing, ] |>
    mutate(numerator = if_else(is.na(numerator), 0, numerator))

  # Rows that we keep as-is (no calc needed)
  df_keep <- all_df[!needs_calc, ] |>
    # mutate(denominator = if_else(is.na(denominator), 0, denominator)) |>
    left_join(get_age_lookup(metadata = metadata2, age_lookup = age_lookup), by = "indicator_id") |>
    mutate(age_group_code = age_code) |>
    select(all_of(TIDY_COLS))

  
  # Print unique indicator ids for which rows are kept as-is 
  print(paste0("Unique indicator ids for which rows are kept as-is (no calc required): ",
               unique(df_keep$indicator_id)))

  # Rows that needed calc but were NOT eligible (will be appended to the output later)
  df_skip <- all_df[needs_calc & !eligible_for_processing, ] |>
    # mutate(denominator = if_else(is.na(denominator), 0, denominator)) |>
    left_join(get_age_lookup(metadata = metadata2, age_lookup = age_lookup), by = "indicator_id") |>
    mutate(age_group_code = age_code) |>
    select(all_of(TIDY_COLS))


  # -------- Percentage ------------------------------------------------------
  message("\u25B6 Calculating percentage...")
  temp1 <- df_calc |> calc_percentage() |>
    tidy_output()

  # -------- Crude rate & Ratio ----------------------------------------------
  message("\u25B6 Calculating crude rate and ratio...")
  temp2 <- df_calc |> calc_crude_or_ratio() |>
    tidy_output()

  # -------- Directly age standardised rate (DASR) ---------------------------
  message("\u25B6 Calculating DASR...")
  temp3 <- df_calc |> calc_dasr(metadata = metadata2, age_lookup = age_lookup) |>
    tidy_output()

  # -------- Count -----------------------------------------------------------
  message("\u25B6 Calculating count...")
  temp4 <- df_calc |> calc_count() |>
    tidy_output()
  
  # -------- Percentage change -------------------------------------------------
  message("\u25B6 Calculating count...")
  temp5 <- df_calc |> calc_percentage_change() |>
    tidy_output()
  
  # -------- Difference --------------------------------------------------------
  message("\u25B6 Calculating difference...")
  temp6 <- df_calc |> calc_difference() |>
    tidy_output()
  
  
  # -------- Combine all outputs ---------------------------------------------
  message("\u25B6 Combining all outputs...")
  output <- bind_rows(
    df_keep,          # Rows not calculated
    df_skip,          # Rows that needed calc but were ineligible (status/value_multiplier)
    temp1,            # Percentage
    temp2,            # Crude rate & Ratio
    temp3,            # DASR
    temp4,            # Count
    temp5,            # Percentage change
    temp6,            # Difference
  ) |> clean_data_types()

  # -------- Reporting of output ---------------------------------------------
  message("\u25B6 Reporting output total rows...")
  print(paste("Total rows for df_keep:", nrow(df_keep)))
  print(paste("Total rows for df_skip:", nrow(df_skip)))
  print(paste("Total rows for temp1 (perc):", nrow(temp1)))
  print(paste("Total rows for temp2 (crude & ratio):", nrow(temp2)))
  print(paste("Total rows for temp3 (DASR):", nrow(temp3)))
  print(paste("Total rows for temp4 (Count):", nrow(temp4)))
  print(paste("Total rows for temp5 (Percentage change):", nrow(temp5)))
  print(paste("Total rows for temp6 (Difference):", nrow(temp6)))

  # -------- Reporting of skipped items --------------------------------------
  message("\u25B6 Reporting skipped items...")
  if (any(needs_calc & !eligible_for_processing)) {
    reasons_df <- all_df[needs_calc & !eligible_for_processing, ] |>
      transmute(
        indicator_id = indicator_id,
        status_code = status_code,
        reason_missing_value_multiplier = is.na(value_multiplier),
        reason_status_not_1 = is.na(status_code) | status_code != 1,
        reason_missing_value_type = is.na(value_type_code)
      ) |>
      group_by(indicator_id, status_code) |>
      summarise(
        reasons = paste0(
          c(
            if (any(reason_missing_value_multiplier)) "missing value_multiplier" else NULL,
            if (any(reason_status_not_1)) paste0("status_code = ", unique(status_code)) else NULL,
            if (any(reason_missing_value_type)) "missing value_type_code" else NULL
          ),
          collapse = "; "
        ),
        .groups = "drop"
      )

    message("\u26A0\uFE0F  Some indicators were not processed:")
    apply(
      reasons_df,
      1,
      function(r) message("  ID ", r[["indicator_id"]], ": ", r[["reasons"]])
    )
  }

  message("\u2705 Process completed!")
  # Return a list of all dfs for tracking
  return(list(combined_calc_dfs = output,
              df_keep = df_keep,
              df_skip = df_skip,
              perc = temp1,
              crude_ratio = temp2,
              dasr = temp3,
              count = temp4,
              perc_change = temp5,
              difference = temp6))
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


# Function to run SQL script ---------------------------------------------------
run_sql_file <- function(conn, path) {
  # Inputs:
  #   conn - SQL connection
  #   path - path to the SQL file
  #
  # Output:
  #   Prints a message indicating the SQL script has run successfully.
  sql_text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  dbExecute(conn, sql_text)
  message("\u2705 SQL script has succesfully run!")
}

# To extract indicators from sql table
get_indicators_from_sql <- function(conn, table_name, indicator_ids = NULL) {
  # Base query
  sql_query <- paste0("SELECT * FROM ", table_name)

  # WHERE clause if IDs provided
  if (!is.null(indicator_ids) && length(indicator_ids) > 0) {
    # Quote each literal safely for this connection (handles numbers/strings)
    quoted <- vapply(indicator_ids, DBI::dbQuoteLiteral, character(1), conn = conn)
    sql_query <- paste0(sql_query, " WHERE Indicator_ID IN (", paste(quoted, collapse = ", "), ")")
  }

  # Execute
  result <- DBI::dbGetQuery(conn, sql_query)
  message("\u2705 Indicators extracted from SQL!")
  print(paste("Total rows for data extracted:", nrow(result)))
  result
}


insert_data_into_sql_table <- function(conn, database, schema, table, data, indicator_ids = NULL,
                                               id_column = "indicator_id") {
  # Delete existing data and append new data
  # Inputs:
  #   conn          - Active SQL connection (e.g., from DBI::dbConnect)
  #   database      - Name of database (as a string)
  #   schema        - Name of schema (as a string)
  #   table         - Name of the SQL table (as a string)
  #   data          - Data frame to append to SQL
  #   indicator_ids - Optional vector of indicator IDs to delete before appending

  # Normalize the ids (NULL or "all"/"*" to delete all rows or certain ids)
  ids <- normalize_indicator_ids(indicator_ids)

  tbl_id  <- DBI::Id(catalog = database, schema = schema, table = table)
  tbl_sql <- DBI::dbQuoteIdentifier(conn, tbl_id)
  col_sql <- DBI::dbQuoteIdentifier(conn, id_column)

  # Build DELETE
  sql_query <- paste0("DELETE FROM ", tbl_sql)
  if (!is.null(ids) && length(ids) > 0) {
    # Quote each literal (works for numbers and strings)
    quoted_vals <- vapply(ids, DBI::dbQuoteLiteral, character(1), conn = conn)
    sql_query <- paste0(sql_query, " WHERE ", col_sql, " IN (", paste(quoted_vals, collapse = ", "), ")")
  } # else: NULL/empty then delete ALL rows in the table

  # Execute DELETE and append
  DBI::dbExecute(conn, sql_query)
  DBI::dbWriteTable(conn, name = tbl_id, value = data, append = TRUE)

  # Messages
  if (is.null(ids) || length(ids) == 0) {
    message("\u2705 All rows deleted and new data appended to ", DBI::dbQuoteIdentifier(conn, DBI::Id(schema=schema, table=table)), ".")
  } else {
    message("\u2705 Deleted and replaced data for selected indicator_ids in ", DBI::dbQuoteIdentifier(conn, DBI::Id(schema=schema, table=table)), ".")
  }
}

# DQ functions------------------------------------------------------------------
## Function to check row counts ------------------------------------------------

check_row_counts <- function(df, reference_data) {
  # Inputs:
  #   input_data     - data frame whose row count will be checked
  #   reference_data - data frame to compare against
  #
  # Output:
  #   Prints a message indicating whether the row counts match (with counts).

  input_rows <- nrow(df)
  reference_rows <- nrow(reference_data)

  if (input_rows == reference_rows) {
    message("\u2705 PASS: Row counts match: ", input_rows, " rows")
  } else {
    message("\u274C FAIL: Row counts do NOT match. ",
            "Input: ", input_rows, " rows | Reference: ", reference_rows, " rows")
  }
}

## Function to check non-populated columns -------------------------------------

# ignoring numerator and denominator as these can be empty columns for pre-calculated indicators

check_rows_with_missing <- function(df, metadata, cols = NULL,
                              ignore = c("numerator", "denominator",
                                         "lower_ci95", "upper_ci95"),
                              show_n = 10
                              ) {
  # Inputs:
  #   df     - data frame whose rows will be checked
  #   metadata - metadata to filter indicators to be checked
  #
  # Output:
  #   A dataframe with columns containing missing rows

  current_ids <- metadata |>
    filter(status_code == 1 & precalculated == "No") |> # only want to check non pre-calculated indicators
    distinct(indicator_id) |>
    pull(indicator_id)

  # columns to check
  cols_to_check <-
    if (is.null(cols)) setdiff(names(df), ignore) else intersect(cols, names(df))

  miss_df <- df |>
    filter(denominator != 0 & !is.na(denominator)) |> # rows required processing
    filter(indicator_id %in% current_ids) |>
    filter(if_any(all_of(cols_to_check), ~ is.na(.x)))

  if (nrow(miss_df) == 0) {
    message("\u2705 PASS: No rows with missing values in the checked columns.")
  } else {
    message("\u274C FAIL: Found rows with missing values in the checked columns: ", nrow(miss_df))
    print(utils::head(miss_df, show_n))
  }

  return(miss_df)

}

## Function to check unique age group code for dasr indicator ------------------

check_age_group_code <- function(df) {
  #  Inputs:
  #   df - data frame containing DASR records; must include
  #        `indicator_id`, `age_group_code`, and `value_type_code`
  #
  # Output:
  #   Prints a message indicating whether each DASR (`value_type_code == 4`)
  #   indicator_id has exactly one unique `age_group_code` (OK) or if any
  #   indicator_id has more than one (problem).

  unique_age_count <- df |>
    filter(value_type_code == 4) |>
    distinct(indicator_id, age_group_code) |>
    group_by(indicator_id) |>
    summarise(count = n()) |>
    arrange(indicator_id) |>
    filter(count > 1)

  if(nrow(unique_age_count) == 0){
    message("\u2705 PASS: One unique age group code per dasr indicator ID")
  } else{
    message("\u274C FAIL: More than 1 age group code per dasr indicator ID")
  }
}

## Function to check all time period types are populated -----------------------
check_time_period_type <- function(df) {
  #  Inputs:
  #   df - data frame whose rows will be checked
  #
  # Output:
  #   Prints a message indicating whether time period column is populated or not
  missing_rows <- df |>
    filter(is.na(time_period_type)) |>
    distinct(indicator_id, start_date, end_date)

  if (nrow(missing_rows) == 0) {
    message("\u2705 PASS: time_period_type is populated for all rows.")
  } else {
    message("\u274C FAIL: Some indicators are missing time_period_type. Details below:")
    print(missing_rows)
  }
}


## Function to check value columns are populated for active indicators ---------

check_active_indicator_values <- function(df, metadata) {
  active_ids <- metadata |> 
    filter(status_code == 1) |> 
    distinct(indicator_id) |> 
    pull(indicator_id)

  # Return rows with missing values in the key indicator_value column
  failures <- df |> 
    filter(indicator_id %in% active_ids,
           is.na(indicator_value))

  if (nrow(failures) == 0) {
    message("\u2705 PASS: All active indicators have populated indicator_value.")
    return(invisible(NULL))  
  } else {
    failed_ids <- failures |> 
      distinct(indicator_id) |> 
      pull(indicator_id)
    
    message("\u274C FAIL: Found active indicators with missing indicator_value. However, this may be acceptable for percentage change metrics, particularly when there are no baseline/previous/plan values for the actuals to be compared against.
             Indicator IDs:", paste(failed_ids, collapse = ", "))
    print(head(failures))
    return(failures)
  }
}

# Function to check missing confidence intervals -------------------------------
check_missing_confidence_intervals <- function(df){
  warnings <- df |> 
    filter(!is.na(indicator_value),
           is.na(lower_ci95) | is.na(upper_ci95)
           )
  
  if(nrow(warnings) == 0){
    message("\u2705 PASS: No missing confidence intervals." )
  } else{
    failed_ids <- warnings |> 
      distinct(indicator_id) |> 
      pull(indicator_id)
    
    message("\u274C FAIL: Some rows have missing confidence intervals, but this may be acceptable.
            Indicator IDs:", paste(sort(failed_ids), collapse = ", "))
    print(head(warnings))
  }
  
  return(warnings)
}

## Function to check combination id is populated -------------------------------
check_combination_splits <- function(df) {
  #  Inputs:
  #   df - data frame whose rows will be checked
  #
  # Output:
  #   Prints a message indicating whether `combination_id` column is populated or not
  missing_rows <- df |>
    filter(is.na(combination_id))

  if (nrow(missing_rows) == 0) {
    message("\u2705 PASS: combination_id is populated for all rows.")
  } else {
    message("\u274C FAIL: Some indicators are missing combination_id. Details below:")
    print(missing_rows)
  }
}

## Function to check duplicates ------------------------------------------------

check_duplicates <- function(df, key_cols=c(
  "indicator_id",
  "time_period_type",
  "aggregation_id",
  "start_date",
  "end_date",
  "imd_code",
  "ethnicity_code"
), indicator_filter = NULL) {

  result <- df %>%

    # Optional indicator filter
    {
      if (!is.null(indicator_filter)) {
        filter(., indicator_id %in% indicator_filter)
      } else {
        .
      }
    } |>

    # Group by key columns
    group_by(across(all_of(key_cols))) |>

    # Count rows
    summarise(row_count = n(), .groups = "drop") |>

    # Keep only duplicates
    filter(row_count > 1) |>

    arrange(indicator_id)

  if (nrow(result) > 0) {
    message(
      "\u274C FAIL: Duplicate records found: ",
      nrow(result),
      " duplicated key combination(s)."
    )
  } else {
    message("\u2705 PASS: No Duplicate records found for the specified key columns")
  }

  return(result)
}

## Function to check source codes ----------------------------------------------

check_source_code <- function(
    df,
    indicator_col = "indicator_id",
    source_code_col = "source_code"){

  # Count distinct source codes per indicator
  results <- df |>
    group_by(.data[[indicator_col]]) |>
    summarise(
      n_source_codes = n_distinct(.data[[source_code_col]]),
      source_codes = paste(unique(.data[[source_code_col]]), collapse = ", "),
      .groups = "drop"
    ) |>
    filter(n_source_codes > 1)

  if(nrow(results) > 0){
    message("\u274C FAIL: Some indicators have more than one source code.")
    return(results)
  }else{
    message("\u2705 PASS: Every indicator has exactly one source code.")
    return(NULL)
  }
}

# df <- data.frame(
#   indicator_id = c(1,1,2,2,3,3),
#   source_code = c(1,2, 3,3, 1, 2 )
# )
#
# check_source_code(df)

# Function to check valid percentages ------------------------------------------
check_percentages <- function(df){
  invalid_rows <- df |>
    filter(value_type_code == 2,
           !is.na(indicator_value),
           indicator_value > 100)

  if(nrow(invalid_rows) > 0){

    failed_indicators <- unique(invalid_rows$indicator_id)
    message("\u274C FAIL: Found percentages greather than 100 for indicator(s):",
            paste(failed_indicators, collapse = ", ")
            )

  } else{
    message("\u2705 PASS: All percentage values are valid.")
  }

  return(invalid_rows)
}

# df <- data.frame(
#   indicator_id = c(1,1,2,2),
#   value_type_code = c(2,2,2,1),
#   indicator_value = c(95,120,105,300)
# )
#
# check_percentages(df)

## Run all DQ checks ------------------------------------------------------------

run_all_dq_checks <- function(df, reference_data, metadata, cols_to_check = NULL, show_n = 10) {
  message("\n==== Running Data Quality Checks ====\n\n")

  message("1) Row counts\n")
  check_row_counts(df, reference_data)
  cat("\n")

  message("2) Rows with missing required columns\n")
  check_rows_with_missing(df, metadata = metadata)
  cat("\n")

  message("3) Unique age_group_code for DASR indicators\n")
  check_age_group_code(df)
  cat("\n")

  message("4) time_period_type populated\n")
  check_time_period_type(df)
  cat("\n")

  message("5) Active indicator values populated\n")
  check_active_indicator_values(df, metadata)
  cat("\n")

  message("6) combination_id populated\n")
  check_combination_splits(df)
  cat("\n")

  message("7) Identify duplicates\n")
  check_duplicates(df)
  cat("\n")

  message("8) Check number of unique source codes \n")
  check_source_code(df)
  cat("\n")

  message("9) Check invalid percentage values \n")
  check_percentages(df)
  
  message("10) Check missing confidence intervals \n")
  check_missing_confidence_intervals(df)

  message("\n==== DQ checks completed ====\n")
}


































