# Mini script to check the number of PFAM genes on GBOT (so a gene with at least 1 PFAM)

PFAMNumber <- function(name, save, number1, number2, number3) {
  
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(ggplot2)))
  
  summary_df <- data.frame(
    PFAM = c("Helixer", "Official", "Both"),
    n = as.numeric(c(number2, number1, number3))
  )
  
  summary_df <- summary_df %>%
    mutate(
      freq = n / sum(n),
      label = paste0(round(freq * 100, 1), "%\n(", n, ")")
    )
  
  # We hate pie charts
  p <- ggplot(summary_df, aes(x = PFAM, y = n, fill = PFAM)) +
    geom_col(width = 0.6) +
    
    geom_text(
      aes(label = label, color = PFAM),
      vjust = -0.4,
      size = 3,
      family = "Courier New",
      show.legend = FALSE
    ) +
    
    scale_fill_manual(values = c(
      "Helixer" = "#63003C",
      "Official" = "#D6DBDF",
      "Both" = "#00a3a6"
    )) +
    
    scale_colour_manual(values = c(
      "Helixer" = "black",
      "Official" = "black",
      "Both" = "black"
    )) +
    
    labs(
      title = paste0("Distribution of Genes with PFAM hit(s) in GBOT | ", name),
      x = "",
      y = paste0("Number of genes (n = ", sum(summary_df$n), ")"),
      fill = NULL
    ) +
    
    theme_bw() +
    theme(
      text = element_text(size = 10, family = "Courier New"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(save, "PFAMDistrib.png"),
    plot = p,
    width = 6.5,
    height = 10,
    dpi = 600
  )
}

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
save_dir <- args[2]
input1 <- args[3]
input2 <- args[4]
input3 <- args[5]

PFAMNumber(name, save_dir, input1, input2, input3)