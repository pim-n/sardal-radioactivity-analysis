### takes the lowercased and underscored names for species in the dataframes,
### and formats them with proper casing and spacing.
### (ex. Tc99_F_Serratus would become Tc99 (F. Serratus) .

format_species_label <- function(label) {
  parts <- strsplit(label, "_")[[1]]
  if (length(parts) == 3) {
    formatted <- paste(parts[1], " (", parts[2], ". ", parts[3], ")", sep = "")
    return(str_to_title(formatted))
  } else {
    return(label)
  }
}

### Same but for discharge sites

format_site_label <- function(label) {
  parts <- strsplit(label, "_")[[1]]
  if (length(parts) == 3) {
    formatted <- paste0(parts[1], " (", parts[2], ")")
    return(str_to_title(formatted))
  } else if (length(parts) == 4) {
    formatted <- paste0(parts[1], " (", parts[2], " ", parts[3], ")")
    return(str_to_title(formatted))
  } else {
    return(label)
  }
}

format_series_name <- function(x) {
  x %>%
    sub("^[^_]+_", "", .) %>%
    gsub("_", ". ", .)
}

### Plotting all data for one isotope.
plot_isotope <- function(df_iso) {
  
  title <- format_series_name(substitute(df_iso)) # extract name of isotope for title
  
  # we un-factor by giving each fucus species its own column
  df_iso <- df_iso %>%
    pivot_wider(
      names_from  = series,
      values_from = sardal,
      names_glue = "{format_series_name(series)}"
    )
  
  # just some fancy labels
  row_labels <- c(
    `F. Serratus` = "F. Serratus\n(Bq / kg d.w.)",
    `F. Serratus. Summer` = "F. Serratus (summer)\n(10^14 atoms / kg d.w.)",
    `F. Serratus. Winter` = "F. Serratus (winter)\n(10^14 atoms / kg d.w.)",
    `F. Vesiculosus` = "F. Vesiculosus\n(Bq / kg d.w.)",
    `F. Various` = "Fucus (various)\n(F14C)",
    sellafield = "Sellafield\n(TBq/y)",
    la_hague = "La Hague\n(TBq/y)",
    winfrith = "Winfrith\n(TBq/y)",
    barseback = "Bärseback\n(TBq/y)",
    baltic = "Baltic sea\n(Bq/m^3)"
  )
  
  # check which cols are present
  source_cols <- intersect(names(row_labels), names(df_iso))
  
  # first we make it long format again, and ensure correct order of levels by mutate()
  df_long <- df_iso %>%
    pivot_longer(
      cols = any_of(source_cols), # pivot on all available columns
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = factor(source, levels = source_cols)#,
      #series = factor(series, levels = levels(df_iso$series))
    )
  
  format_labels <- function(x) {
    gsub("_", " ", x)
  }
  
  p <- ggplot(df_long, aes(x = year, y = value)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    facet_wrap(
      vars(source),
      ncol = 1,
      scales = "free_y",
      labeller = labeller(source = row_labels),
      strip.position = "left"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9),
      strip.placement = "outside",
    ) +
    labs(x = NULL, y = NULL) +
    ggtitle(title)
  return(p)
}

### Plotting code for the NA values

plot_na_values <- function(df_long) {
  tmp <- df_long %>%
    mutate(tf = !is.na(value))
  
  # if obs_date, it is särdal. otherwise it is discharges
  if("obs_date" %in% colnames(tmp)) {
    tmp_yearly <- tmp %>%
      mutate(year = format(obs_date, "%Y")) %>%
      group_by(variable, year) %>%
      summarize(data_available = any(tf), .groups = "drop")
    
    tmp_yearly <- tmp_yearly %>%
      mutate(formatted_variable = sapply(variable, format_species_label))
  }
  else {
    tmp_yearly <- tmp %>%
      group_by(variable, year) %>%
      summarize(data_available = any(tf), .groups = "drop") %>%
      mutate(formatted_variable = sapply(variable, format_site_label))
  }
  
  min_year <- min(as.numeric(tmp_yearly$year))
  max_year <- max(as.numeric(tmp_yearly$year))
  
  # Create breaks for every 5 years (adjust the 'by' parameter as needed)
  breaks <- seq(min_year, max_year, by = 5)
  
  ggplot(tmp_yearly, aes(x = year, y = formatted_variable, fill = data_available)) +
    geom_tile() +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#555555")) +
    labs(x = "Year", y = "Series", fill = "Data Available") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),  # Adjust the size of y-axis labels
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better 
    scale_x_discrete(breaks = breaks)
}