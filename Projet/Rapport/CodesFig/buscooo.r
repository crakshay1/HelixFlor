library(gt)

busco_table <- data.frame(
  Subset = c("Common HELIXER", "Common TAIR12", "TAIR12-specific (Missed)", "HELIXER-specific (Added)"),
  Complete = c(98.9, 99.5, 0.2, 0.1),
  Single_copy = c(98.0, 79.6, 0.2, 0.1),
  Duplicated = c(0.9, 19.9, 0.0, 0.0),
  Fragmented = c(0.4, 0.1, 0.0, 0.2),
  Missing = c(0.7, 0.4, 99.8, 99.7)
)

busco_table |>
  gt() |>
  cols_label(
    Subset = "Subset",
    Complete = "Complete (C)",
    Single_copy = "Single-copy (S)",
    Duplicated = "Duplicated (D)",
    Fragmented = "Fragmented (F)",
    Missing = "Missing (M)"
  ) |>
  fmt_number(
    columns = -Subset,
    decimals = 1,
    pattern = "{x}%"
  ) |>
  data_color(
    columns = Complete,
    palette = c("white", "#2a78d6")
  ) |>
  data_color(
    columns = Missing,
    palette = c("white", "#e34948")
  ) |>
  tab_header(
    title = "BUSCO Assessment Results",
    subtitle = "Lineage : embryophyta_odb10 | n = 1614"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  opt_stylize(style = 1)