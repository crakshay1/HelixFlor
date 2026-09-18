# Added Missed Commons stats 

AddedVSMissedVSCommon <- function(qcover, threshold_ident, name, save, Added ,Missed, Common, pfama, pfamm, pfamc, resC, resA, resM, iso) {
  suppressMessages(
    suppressWarnings(library(readr)))
  suppressMessages(
    suppressWarnings(library(dplyr)))
  suppressMessages(
    suppressWarnings(library(ggplot2)))
  suppressMessages(
    suppressWarnings(library(gridExtra)))
  
  
  add <- suppressMessages(
    suppressWarnings(
      readr::read_delim(Added, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
    )
  )
  
  miss1 <- suppressMessages(
    suppressWarnings(
      readr::read_delim(Missed, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
    )
  )
  
  comm1 <- suppressMessages(
    suppressWarnings(
      readr::read_delim(Common, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
    )
  )
  misses <- sub("\\.[0-9]+$", "", miss1$Reference)  # enlève .1, .2, .3, etc.
  comms <- sub("\\.[0-9]+$", "", comm1$Reference)  # enlève .1, .2, .3, etc.
  if (iso != "True") {
    comm <- comm1[!duplicated(comms),]  
    miss <- miss1[!duplicated(misses),] 
  } else {
    miss <- miss1
    comm <- comm1
  }
  

  # FIRST LEVEL : Median Length
  ALength <- median(add$Length)
  MLength <- median(miss$Length)
  CLength <- median(comm$Length)
  
  funnel_df <- data.frame(
    Step = factor(c("Added", "Covered", "Missed"),
                  levels = c("Added", "Covered", "Missed")),
    Count = c(ALength, CLength, MLength)
  )
  
  pA <- ggplot(funnel_df, aes(x = Step, y = Count, fill = Step)) +
    geom_col() +
    
    geom_text(
      aes(label = paste0("(", Count, ")")),
      vjust = -0.4,
      size = 3
    ) +
    
    scale_fill_manual(values = c(
      "Covered" = "#00a3a6",
      "Added" = "#63003C",
      "Missed" = "#D6DBDF"
    )) +
    
    theme_bw() +
    labs(
      title = paste0("A - Median Length of Helixer CDS | ", name),
      subtitle = paste0("Added genes : ", nrow(add) , " | ",
                        "Missed genes : ", nrow(miss), " | ",
                        "Common genes : ", nrow(comm)),
      x = "",
      y = "Median Length (bp)"
    ) +
    theme(
      text = element_text(size = 8),
      legend.position = "none"
    )
  
  
  # SECOND LEVEL : Exons 
  ALength <- round(median(add$NbExons),1)
  MLength <- round(median(miss$NbExons),1)
  CLength <- round(median(comm$NbExons),1)
  
  funnel_df <- data.frame(
    Step = factor(c("Added", "Covered", "Missed"),
                  levels = c("Added", "Covered", "Missed")),
    Count = c(ALength, CLength, MLength)
  )
  
  pB <- ggplot(funnel_df, aes(x = Step, y = Count, fill = Step)) +
    geom_col() +
    
    geom_text(
      aes(label = paste0("(", Count, ")")),
      vjust = -0.4,
      size = 3
    ) +
    
    scale_fill_manual(values = c(
      "Covered" = "#00a3a6",
      "Added" = "#63003C",
      "Missed" = "#D6DBDF"
    )) +
    
    theme_bw() +
    labs(
      title = paste0("B - Median Exons Number in Helixer CDS | ", name),
      subtitle = paste0("Added genes : ", nrow(add) , " | ",
                        "Missed genes : ", nrow(miss), " | ",
                        "Common genes : ", nrow(comm)),
      x = "",
      y = "Number of Exons"
    ) +
    theme(
      text = element_text(size = 8),
      legend.position = "none"
    )
  
  
  # THIRD LEVEL : PFAM
  summary_df <- data.frame(
    PFAM = c("Covered", "Added", "Missed"),
    n = as.numeric(c(pfamc, pfama, pfamm))
  )
  
  summary_df <- summary_df %>%
    mutate(
      freq = n / sum(n),
      label = paste0("(",n,")")
    )
  
  # We hate pie charts
  pC <- ggplot(summary_df, aes(x = PFAM, y = n, fill = PFAM)) +
    geom_col() +
    geom_text(
      aes(label = label, color = PFAM),
      vjust = -0.4,
      size = 3,
      show.legend = FALSE
    ) +
    
    scale_fill_manual(values = c(
      "Added" = "#63003C",
      "Missed" = "#D6DBDF",
      "Covered" = "#00a3a6"
    )) +
    
    scale_colour_manual(values = c(
      "Added" = "black",
      "Missed" = "black",
      "Covered" = "black"
    )) +
    
    labs(
      title = paste0("C - Distribution of Genes with PFAM hit(s) in GBOT | ", name),
      x = "",
      y = paste0("Number of genes with PFAM hits (n = ", sum(summary_df$n), " / ",nrow(comm)+nrow(add)+nrow(miss),")"),
      fill = NULL
    ) +
    
    theme_bw() +
    theme(
      text = element_text(size = 8),
      legend.position = "none"
    )
  
    # FOURTH LEVEL : BLAST
    evalue_threshold <- 1e-5
    machin <- function(file) {
      lines <- readLines(file)
      
      fixed <- character()
      current <- lines[1]
      
      for (i in 2:length(lines)) {
        
        ncols <- length(strsplit(current, "\t")[[1]])
        
        if (ncols < 12) {
          current <- paste0(current, " ", trimws(lines[i]))
        } else {
          fixed <- c(fixed, current)
          current <- lines[i]
        }
      }
      
      fixed <- c(fixed, current)
      bak_file <- paste0(file, ".bak")
      file.copy(file, bak_file, overwrite = TRUE)
      
      writeLines(fixed, file)
    }
    machin(resA)
    machin(resM)
    machin(resC)
    
    blastCommon <- suppressMessages(
      suppressWarnings(
        readr::read_delim(resC, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
      )
    )
    blastCommon <- blastCommon[
      (  blastCommon$Gene %in% comm$Helixer &
         blastCommon$PIdent >= threshold_ident &
         blastCommon$QCover >= threshold_ident &
         blastCommon$Evalue <= evalue_threshold) |
        is.na(blastCommon$PIdent),
    ]
    blastCommon$Uniprot <- ifelse(
    is.na(blastCommon$Evalue),
    "No hit",
    "Uniprot Hit"
    )

    blastAdded <- suppressMessages(
    suppressWarnings(
      readr::read_delim(resA, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
        )
    )
    
    blastAdded <- blastAdded[
      (blastAdded$PIdent >= threshold_ident &
         blastAdded$QCover >= threshold_ident &
         blastAdded$Evalue <= evalue_threshold) |
        is.na(blastAdded$PIdent),
    ]
    blastAdded$Uniprot <- ifelse(
    is.na(blastAdded$Evalue),
    "No hit",
    "Uniprot Hit"
    )

    blastMissed <- suppressMessages(
    suppressWarnings(
      readr::read_delim(resM, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
        )
    )
    
    blastMissed <- blastMissed[
      ( blastMissed$Gene %in% miss$Reference &
        blastMissed$PIdent >= threshold_ident &
         blastMissed$QCover >= threshold_ident &
         blastMissed$Evalue <= evalue_threshold) |
        is.na(blastMissed$PIdent),
    ]
    blastMissed$Uniprot <- ifelse(
    is.na(blastMissed$Evalue),
    "No hit",
    "Uniprot Hit"
    )

    summary_df <- data.frame(
    Category = c(
        rep("Covered", nrow(blastCommon)),
        rep("Added", nrow(blastAdded)),
        rep("Missed", nrow(blastMissed))
    ),
    Status = c(
        blastCommon$Uniprot,
        blastAdded$Uniprot,
        blastMissed$Uniprot
    )
    ) %>%
    count(Category, Status)

    summary_df <- summary_df %>%
    group_by(Category) %>%
    mutate(
        Total = sum(n),
        Label = paste0("(",n,")")
    ) %>%
    ungroup()

    pD <- ggplot(summary_df, aes(x = Category, y = n, fill = Status)) +
    geom_col(width = 0.7) +
    geom_text(
        aes(label = Label, color = Status),
        position = position_stack(vjust = 0.5),
        size = 2.5
    ) +
      scale_fill_manual(values = c(
        "Uniprot Hit" = "#0066CC",   # UniProt-associated blue
        "No hit" = "#D6DBDF"
      )) +
      scale_color_manual(values = c(
        "Uniprot Hit" = "white",
        "No hit" = "black"
      )) +
    guides(
        fill = guide_legend(title = "Uniprot Proteins"),
        color = "none"
      ) +

    labs(
        title = paste0(
        "D - Identification of Helixer proteins in Uniprot", " | ", name),
        subtitle = paste0("Identity ≥ ", threshold_ident , "% | ",
                          "Qcover ≥ ", qcover, "% | ",
                          "Evalue ≤ ", evalue_threshold),
        x = "",
        y = paste0("Number of proteins (n = ", sum(summary_df$n),")"),
        fill = NULL
    ) +

    theme_bw() +
    theme(
        text = element_text(size = 8),
        legend.position = "right"
    )

    
  # Final figure ahh
  suppressMessages(
    suppressWarnings(
          final_fig <- grid.arrange(
          pA, pB, pC, pD,
          ncol = 2
      )
    )
  )
  
  ggsave(file.path(save, "../figs/MissedAddedCommon.png"),
         final_fig,
         width = 10,
         height = 8.8,
         dpi = 600)
  
}
# Rscript Projet/Rapport/CodesFig/ProteinAnalysisSpecies.r TAIR12 Projet/Species/TAIR12 Projet/Species/TAIR12/Protein_Analysis.tsv Projet/Species/TAIR12/Arabidopsis_thaliana_TAIR12_pfamHelixer.txt 

args <- commandArgs(trailingOnly = TRUE)

qcover <- as.numeric(args[1])
pident <- as.numeric(args[2])
name <- args[3]
save_dir <- args[4]
added <- args[5]
missed <- args[6]
common <- args[7]
input1 <- args[8]
input2 <- args[9]
input3 <- args[10]
resC <- args[11]
resA <- args[12]
resM <- args[13]
isoform <- args[14]

result <- AddedVSMissedVSCommon(qcover, pident, name, save_dir,
                                added, missed, common,
                                input1, input2, input3,
                                resC, resA, resM, isoform)