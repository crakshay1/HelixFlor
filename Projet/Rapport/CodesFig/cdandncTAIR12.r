library(gt)
library(tidyverse)

# 1. Create the dataset
tair12_data <- tibble(
  Category = c(
    "Coding genes", 
    "Non-coding genes", "nc_RNA", "tRNA", "miscRNA", "precursor", "rRNA",
    "TOTAL genes"
  ),
  Count = c(
    26867, 
    13580, 5602, 526, 26, 881, 6545,
    40447
  ),
  Type = c(
    "Primary", 
    "Primary", "Subtype", "Subtype", "Subtype", "Subtype", "Subtype",
    "Total"
  )
)

# 2. Build the gt table
tair12_table <- tair12_data %>%
  gt() %>%
  tab_header(
    title = md("**TAIR12 Genome Annotation Summary**"),
    subtitle = "Gene count distribution extracted from GBOT"
  ) %>%
  tab_style(
    style = cell_text(indent = px(20), style = "italic", color = "gray30"),
    locations = cells_body(
      columns = Category,
      rows = Type == "Subtype"
    )
  ) %>%
  # Highlight the Total and Primary rows
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = everything(),
      rows = Type %in% c("Primary", "Total")
    )
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#f9f9f9"),
      cell_borders(sides = "top", color = "black", weight = px(2))
    ),
    locations = cells_body(
      columns = everything(),
      rows = Category == "TOTAL genes"
    )
  ) %>%
  # Format numbers 
  fmt_number(
    columns = Count,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  # Clean up columns
  cols_hide(columns = Type) %>%
  cols_label(
    Category = md("**Feature Type**"),
    Count = md("**Gene Count**")
  )

# Display the table
tair12_table
