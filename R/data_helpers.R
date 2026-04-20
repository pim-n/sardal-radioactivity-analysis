### From the file path, obtain species and isotope strings

get_species_name <- function(path) {
  str_remove_all(path, glue("^{sardal_path}/.+sardal_activity_")) |> str_remove("_Mattsson_2025.csv$")
}

get_isotope_name <- function(path) {
  str_remove_all(path, glue("^{discharges_path}/discharges_")) |> str_remove(".csv$")
}

### From the list of time series dfs, get a full date sequence from
### the earliest to latest observation across all dfs

get_full_time_range_by_day <- function(list_of_df) {

  min_date <- as_date(min(unlist(map(list_of_df, ~min(.x[[1]])))))
  max_date <- as_date(max(unlist(map(list_of_df, ~max(.x[[1]])))))
  seq(min_date, max_date, by = "day")
}

get_full_time_range_by_year <- function(list_of_df) {
  min_year <- min(unlist(map(list_of_df, ~min(.x[[1]]))))
  max_year <- max(unlist(map(list_of_df, ~max(.x[[1]]))))
  seq(from = min_year, to = max_year)
}

### Interpolation function for discharges. mvgam does not allow NA values in
### predictors (as of December 2025)

interpolate_data <- function(isotope, df_long_discharges) {
  # subset only the isotope specified
  df_subset <- df_long_discharges %>%
    filter(str_detect(variable, isotope))
  
  min_max <- df_subset %>%
    filter(!is.na(value)) %>%
    group_by(variable) %>%
    summarise(min = min(year),
              max = max(year))
  
  # Apply interpolation between min and max dates for each variable
  df_interpolated <- df_subset %>%
    group_by(variable) %>%
    mutate(value = ifelse(
      year >= min_max$min[min_max$variable == first(variable)] &
        year <= min_max$max[min_max$variable == first(variable)],
      na.approx(value, rule = 2),
      value
    )) %>%
    ungroup()
  
  return(df_interpolated)
}

time_weighted_mean <- function(date, value) {
  o <- order(date)
  date <- date[o]
  value <- value[o]
  
  t <- as.numeric(date)
  
  if(length(value) == 1) return(value)  # fallback
  
  dt <- diff(t)
  
  # trapezoidal integration
  integral <- sum((value[-1] + value[-length(value)]) / 2 * dt)
  
  year_length <- max(t) - min(t)
  
  integral / year_length
}

### Combination of discharge and measurement.
### This code summarizes Särdal to a yearly scale and returens unified df.

combine_data <- function(isotope, df_long_sardal, df_long_discharges, scale_data = FALSE) {
  
  # filter out the isotope in measurements & summarize to yearly values.
  tmp_sardal <- df_long_sardal %>%
    filter(str_detect(variable, isotope)) %>%
    mutate(year = lubridate::year(obs_date),
           variable = factor(variable)) %>%
    distinct() %>%
    group_by(variable, year) %>%
    summarize(sardal = mean(value, na.rm=TRUE), .groups = "drop")
  
  # select isotope from discharges, and pivot to make discharge columns
  tmp_discharges <- df_long_discharges %>%
    filter(str_detect(variable, isotope)) %>%
    interpolate_data(isotope, .) %>%
    filter(!is.na(value)) %>%
    mutate(variable = factor(variable)) %>%
    distinct() %>%
    pivot_wider(id_cols = "year", names_from = "variable", values_from = "value")
  
  # join the two temporary data structures.
  # rename_with is flexible and won't throw an error if columns are not extisting.
  # this is why we use it for bärseback and winfrith.
  return_data <- tmp_discharges %>%
    left_join(tmp_sardal, by=c('year')) %>%
    mutate(time = year - min(year) + 1) %>%
    rename(series = variable,
           sellafield = contains("sellafield"),
           la_hague = contains("la_hague"))
  
  if (any(grepl("winfrith", names(return_data)))) {
    return_data <- return_data %>%
      rename_with(~ "winfrith", contains("winfrith")) %>%
      rename_with(~ "barseback", contains("barseback")) 
  }
  
  if (scale_data) {
    return_data <- return_data %>% 
      group_by(series) %>%
      mutate(sardal = as.numeric(sardal),
             sellafield = as.numeric(scale(sellafield)),
             la_hague = as.numeric(scale(la_hague)))
  }
  
  return_data
}

combine_data_2 <- function(isotope, df_long_sardal, df_long_discharges, scale_data = FALSE) {
  
  # filter out the isotope in measurements & summarize to yearly values.
  tmp_sardal <- df_long_sardal %>%
    filter(str_detect(variable, isotope)) %>%
    mutate(year = lubridate::year(obs_date),
           variable = factor(variable)) %>%
    distinct() %>%
    #group_by(variable, year) %>%
    rename(sardal = value)
    #summarize(sardal = mean(value, na.rm=TRUE), .groups = "drop")
  
  # select isotope from discharges, and pivot to make discharge columns
  tmp_discharges <- df_long_discharges %>%
    filter(str_detect(variable, isotope)) %>%
    interpolate_data(isotope, .) %>%
    filter(!is.na(value)) %>%
    mutate(variable = factor(variable)) %>%
    distinct() %>%
    pivot_wider(id_cols = "year", names_from = "variable", values_from = "value")
  
  # join the two temporary data structures.
  # rename_with is flexible and won't throw an error if columns are not extisting.
  # this is why we use it for bärseback and winfrith.
  return_data <- tmp_sardal %>%
    left_join(tmp_discharges, by=c('year')) %>%
    mutate(time = year - min(year) + 1) %>%
    rename(series = variable,
           sellafield = contains("sellafield"),
           la_hague = contains("la_hague"))
  
  if (any(grepl("winfrith", names(return_data)))) {
    return_data <- return_data %>%
      rename_with(~ "winfrith", contains("winfrith")) %>%
      rename_with(~ "barseback", contains("barseback")) 
  }
  
  if (scale_data) {
    return_data <- return_data %>% 
      group_by(series) %>%
      mutate(sardal = as.numeric(sardal),
             sellafield = as.numeric(scale(sellafield)),
             la_hague = as.numeric(scale(la_hague)))
  }
  
  return_data
}

### Creating the lagged matrices
add_lagged_matrices <- function(data, n_max = 6, n_min = 0) {
  
  stopifnot(n_min >= 0, n_min < n_max)
  
  n_lag <- n_max - n_min
  
  # identify predictor columns to lag
  ignore_set <- c("time", "year", "sardal", "series", "bomb_pulse")
  preds_to_lag <- setdiff(names(data), ignore_set)
  
  # function to create lag number matrices
  lagard <- function(x, n_lag) {
    n <- length(x)
    X <- matrix(NA, n, n_lag)
    for (i in seq_len(n_lag)) {
      X[i:n, i] <- x[1:(n - i + 1)]
    }
    X
  }
  
  # compute lag matrices for each predictor in data
  lagged_preds <- lapply(preds_to_lag, function(pred) {
    
    do.call(
      rbind,
      lapply(seq_along(levels(data$series)), function(x) {
        
        tempdat <- data %>%
          dplyr::filter(series == levels(data$series)[x]) %>%
          dplyr::arrange(time) %>%
          dplyr::pull(!!rlang::sym(pred))
        
        lag_mat <- lagard(tempdat, n_lag)
        tail(lag_mat, NROW(lag_mat) - n_max + 1)
      })
    )
  })
  
  names(lagged_preds) <- preds_to_lag
  
  # align base data
  data_trim <- data %>%
    dplyr::arrange(series, time) %>%
    dplyr::filter(time > n_max - 1)
  
  # assemble final output
  data_all <- c(
    list(
      lag = matrix(
        seq(n_min, n_max - 1),
        nrow = nrow(data_trim),
        ncol = n_max - n_min,
        byrow = TRUE
      ),
      sardal = data_trim$sardal,
      time = data_trim$time,
      series = data_trim$series
    ),
    lagged_preds
  )
  
  if ("bomb_pulse" %in% names(data_trim)) {
    data_all$bomb_pulse <- data_trim$bomb_pulse
  }
  
  return(data_all)
}

### Apply weights for hierarchy

apply_weights <- function(data_all) {
  weights_s <- weights_v <- weights_var <- weights_ss <- weights_sw <-
    matrix(1, ncol = ncol(data_all$lag), nrow = nrow(data_all$lag))
  
  weights_s[!grepl("Serratus", data_all$series, ignore.case = TRUE), ] <- 0
  weights_v[!grepl("Vesiculosus", data_all$series, ignore.case = TRUE), ] <- 0
  weights_var[!grepl("Various", data_all$series, ignore.case = TRUE), ] <- 0
  weights_ss[!grepl("Summer", data_all$series, ignore.case = TRUE), ] <- 0
  weights_sw[!grepl("Winter", data_all$series, ignore.case = TRUE), ] <- 0
  
  data_all$weights_s <- weights_s
  data_all$weights_v <- weights_v
  data_all$weights_var <- weights_var
  data_all$weights_ss <- weights_ss
  data_all$weights_sw <- weights_sw
  data_all
}