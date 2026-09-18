BLASTnSummary <- function(directory) {
  suppressMessages(
    suppressWarnings(
      library(readr))
  )
  suppressMessages(
    suppressWarnings(
      library(dplyr))
  )
  suppressMessages(
    suppressWarnings(
      library(ggplot2))
  )
    suppressMessages(
    suppressWarnings(
      library(gt))
  )

  input_file  <- paste0(directory,"/found_TEs.tsv")  # Sacré BlastN
  input_file1  <- paste0(directory,"/Protein_Analysis_Added.tsv")  # Sacré BlastN
  final_tes <- read_delim(input_file, delim = "\t", col_types = cols())
  blast_results <- read_delim(input_file1, delim = "\t", col_types = cols())

  df_true4 <- final_tes[
    final_tes$QCover == 100 &
    final_tes$PIdent > 90 &
    final_tes$TEChr != "Unknown",
  ]

  write.table(
    df_true4,
    file = paste0(directory, "/TRUEfiltered_TEs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary_table <- data.frame(
    Metric = c(
      "Number of genes",
      "E-value range",
      "QCover range",
      "Identity range"
    ),
    
    Value = c(
      nrow(final_tes),
      paste0(min(final_tes$Evalue, na.rm = TRUE), " - ", round(max(final_tes$Evalue, na.rm = TRUE), 46)),
      paste0(min(final_tes$QCover, na.rm = TRUE), " - ", max(final_tes$QCover, na.rm = TRUE), "%"),
      paste0(round(min(final_tes$PIdent, na.rm = TRUE), 1), " - ", max(final_tes$PIdent, na.rm = TRUE), "%")
    )
  )

  gt(summary_table) %>%
    tab_header(
      title = "BLASTn Summary Statistics",
      subtitle = "Global ranges across all Helixer-specific predictions having the same genomic positions as TEs"
    ) %>%
    cols_align(
      align = "center",
      columns = everything()
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    )


  blast_results$Uniprot <- ifelse(is.na(blast_results$Evalue), "No hit", "Uniprot Hit")
  blast_results$TEs <- ifelse(blast_results$Gene %in% final_tes$Gene,
                                "Present", "Absent")
  blast_filtered <- blast_results %>%
    filter(
      (is.na(QCover) | QCover > 60) &
        (is.na(PIdent) | PIdent > 70) &
        (is.na(Evalue) |Evalue < 1e-5)
    )
  dfB <- blast_filtered %>%
    count(Uniprot, TEs) %>%
    group_by(Uniprot) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  pB <- ggplot(dfB, aes(x = Uniprot, y = n, fill = TEs)) +
    geom_col() +
    geom_text(
      aes(
        label = paste0(round(pct, 1), "%","\n(", n, "/", nrow(blast_results), ")"),
        color = TEs
      ),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(values = c(
      "Present" = "#63003C",
      "Absent" = "#D6DBDF"
    )) +
    scale_color_manual(values = c(
      "Present" = "#D6DBDF",
      "Absent" = "#63003C"
    )) +
    guides(
      fill = guide_legend(title = "Found TEs"),
      color = "none"
    ) +
    theme_bw() +
    labs(
      title = "Association between TEs overlap and absence of UniProt protein hits",
      subtitle = paste0(
        "QCover filtered > ", 60, "%",
        " | Evalue <", 1e-5,
        " | Identity filtered >", 70, "%",
        " | Total proteins: ", length(unique(blast_results$Gene)),
        " | Hits: ", sum(blast_filtered$Uniprot == "Uniprot Hit"),
        " | TEs: ",sum(blast_filtered$TEs == "Present"), "/" , nrow(final_tes[unique(final_tes$Gene),])
      ),
      x = "",
      y = ""
    )

  ggsave(file.path(directory, "/figs/TEsAnalysisEveryLevel.png"),
        pB,
        width = 10,
        height = 7,
        dpi = 600)
}

args <- commandArgs(trailingOnly = TRUE)

chemin <- args[1]
result <- BLASTnSummary(chemin)