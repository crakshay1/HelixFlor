CDS_type <- function(dir, utr) {
  
  library(ggplot2)
  library(dplyr)
  
  directory <- "~/SUSHI/Work/GBOT/flagv5.0/"
  resume <- read.delim(paste0(directory, dir, "StatsSAMECDS.txt"))
  
  df_tab <- as.data.frame(table(resume))
  colnames(df_tab) <- c("Category", "Freq")
  
  df_tab$Category <- factor(df_tab$Category, levels = df_tab$Category)
  
  my_cols <- c(
    "=" = "#63003C",
    "~" = "#00a3a6",
    "c" = "#D6DBDF",
    "j" = "#2C3E50",
    "k" = "#8E44AD",
    "m" = "#E67E22",
    "n" = "#95A5A6"
  )
  
  df_tab$perc <- df_tab$Freq / sum(df_tab$Freq) * 100
  
  x <- ggplot(df_tab, aes(x = Category, y = Freq, fill = Category)) +
    
    geom_col() +
    
    geom_text(
      aes(label = paste0(round(perc, 1), "%\n(", Freq, ")")),
      fontface = "bold",
      vjust = -0.3,
      size = 3
    ) +
    
    scale_fill_manual(values = my_cols) +
    labs(
      title = "CDS Structural Comparison: Helixer vs Reference Annotation",
      x = "Type",
      y = "Count"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        face = "bold",
        size = 14
      )
    )
  
  ggsave(
    paste0(directory, dir, "figs/CDS_stats.png"),
    plot = x,
    dpi = 600,
    width = 8,
    height = 10
  )
  
  resume2 <- read.delim(paste0(directory, dir, utr))
  tmp <- table(resume2$Type.Egalité)
  df_tab2 <- data.frame(
    Category = names(tmp),
    Freq = as.numeric(tmp)
  )
  
  df_tab2$Category <- factor(df_tab2$Category, levels = df_tab2$Category)
  
  my_cols <- c(
    "=" = "#63003C",
    "~" = "#00a3a6",
    "c" = "#D6DBDF",
    "j" = "#2C3E50",
    "k" = "#8E44AD",
    "m" = "#E67E22",
    "n" = "#95A5A6"
  )
  
  df_tab2$perc <- df_tab2$Freq / sum(df_tab2$Freq) * 100
  
  y <- ggplot(df_tab2, aes(x = Category, y = Freq, fill = Category)) +
    
    geom_col() +
    
    geom_text(
      aes(label = paste0(round(perc, 1), "%\n(", Freq, ")")),
      fontface = "bold",
      vjust = -0.3,
      size = 3
    ) +
    
    scale_fill_manual(values = my_cols) +
    labs(
      title = "mRNA Structural Comparison for identical CDS structures",
      x = "Type",
      y = "Count"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        face = "bold",
        size = 14
      )
    )
  
  ggsave(
    paste0(directory, dir, "figs/mRNA_stats.png"),
    plot = y,
    dpi = 600,
    width = 8,
    height = 10
  )
}



args <- commandArgs(trailingOnly = TRUE)

hey <- as.numeric(args[1])
ut <- args[2]
hey <- "Projet/Species/ARAPORT11/"
ut <- "UTRAnalysis/Arabidopsis_thaliana_UTR_ForSameCDS.tsv"

result <- CDS_type(hey, ut)