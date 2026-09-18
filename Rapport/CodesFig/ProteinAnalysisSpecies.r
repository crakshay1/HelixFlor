# Goal here is filtering out proteins with no blastp hits, 
# and also checking if the Helixer proteins with PFAM motifs 
# are actually accurate... (Blud switched to english)

Protein_AnalysisSpecies <- function(name, save, threshold_ident, directory1, directory2 = NULL, directory3 = NULL) {
  suppressMessages(
    suppressWarnings(library(readr)))
  suppressMessages(
    suppressWarnings(library(dplyr)))
  suppressMessages(
    suppressWarnings(library(ggplot2)))
  suppressMessages(
    suppressWarnings(library(gridExtra)))
  
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
  machin(directory1)
  blast_results <- suppressMessages(
    suppressWarnings(
      readr::read_delim(directory1, delim = "\t", escape_double = FALSE, show_col_types = FALSE)
    )
  )
  # Species Name, Save Directory, Blast results, PFAM proteins -> Args 
  if (!is.null(directory2)) { # Given only if the file is not empty !
    pfam <- read.delim(directory2, header = FALSE)
    commons <- inner_join( # PFAM Proteins only
      blast_results,
      pfam,
      by = c("Gene" = "V2")
    )
    # Order by Evalue
    df_sorted1 <- commons[order(commons$Evalue), ]
    df_true1 <- df_sorted1 %>%
                filter(!is.na(df_sorted1$Evalue)) %>%
                distinct(Gene, .keep_all = TRUE) # Getting rid of proteins with no hit
    write.table(df_true1,
              file = file.path(save, "../BlastAnalysisADDED/FilteredPFAMProteinAnalysis.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  }
  
  df_sorted <- blast_results[order(blast_results$Evalue), ] # Order by Evalue
  df_true <- df_sorted[!is.na(df_sorted$Evalue), ] # Getting rid of proteins with no hit
  df_true <- df_sorted[df_sorted$PIdent >= threshold_ident,]

  # Merging EggNogs and Blast results if no Uniprot match
  if (!is.null(directory3)) {
    # Gene	CDSLength	ExonsNb	Protein	GO Annot	Function	ALength	QCover	Score	Evalue	PIdent	Species
    # #query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
    eggnog <- read.delim(directory3, header = TRUE, sep = "\t") %>%
      select(query, seed_ortholog, GOs, Description, score, evalue, max_annot_lvl) %>%
      group_by(query) %>%
      slice_min(evalue, n = 1, with_ties = FALSE) %>%
      ungroup()

    no_match <- df_sorted %>%
      filter(is.na(Evalue)) %>%
      select(Gene,CDSLength,ExonsNb)
    eggnog$clean_query <- sub("_[^_]*$", "", eggnog$query)
    eggnog <- eggnog %>% select(-query)
    giving_match <- inner_join(no_match, eggnog, by = c("Gene" = "clean_query"))
    giving_match <- giving_match[order(giving_match$evalue), ]
    write.table(giving_match,
              file = file.path(save, "../BlastAnalysisADDED/FilteredEggNogProtein.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  } 

  # Filtering 
  df_true2 <- df_true[df_true$QCover > 50, ]
  df_true3 <- df_true[df_true$QCover > 60, ]
  df_true4 <- df_true[df_true$QCover > 70, ]
  
  write.table(df_true,
              file = file.path(save, "../BlastAnalysisADDED/FilteredProteinAnalysisNoFilter.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  write.table(df_true2,
              file = file.path(save, "../BlastAnalysisADDED/FilteredProteinAnalysis70.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  write.table(df_true3,
              file = file.path(save, "../BlastAnalysisADDED/FilteredProteinAnalysis60.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
  write.table(df_true4,
              file = file.path(save, "../BlastAnalysisADDED/FilteredProteinAnalysis50.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  # ===================================================================================
  # Below is the figure part           
  
  # FIRST LEVEL : Uniprot
  blast_results$Uniprot <- ifelse(is.na(blast_results$Evalue), "No hit", "Uniprot Hit")
  total <- nrow(blast_results) 
  hit <- sum(blast_results$Uniprot == "Uniprot Hit") 
  no_hit <- sum(blast_results$Uniprot == "No hit")
  
  blast_filtered <- blast_results %>%
    filter(
      is.na(QCover) | QCover > threshold_ident
    )
  pA <- ggplot(blast_filtered, aes(x = PIdent)) +
    
    # zone <= 50% (background)
    annotate("rect",
             xmin = -Inf, xmax = 50,
             ymin = -Inf, ymax = Inf,
             fill = "#D6DBDF", alpha = 0.4) +
    
    geom_histogram(
      bins = 30,
      fill = "#63003C",
      color = "white"
    ) +
    
    geom_vline(xintercept = c(50, 70, 90),
               linetype = "dashed",
               color = "#2C3E50") +
    
    annotate("text",
             x = 25,
             y = Inf,
             label = "≤50% identity",
             vjust = 2,
             size = 3.5,
             color = "#2C3E50") +
    
    theme_bw(base_size = 12) +
    
    labs(
      title = "A - Identity distribution of Helixer proteins in Uniprot",
      subtitle = paste0(
        "QCover filtered > ", threshold_ident, "%",
        " | Total proteins: ", length(unique(blast_results$Gene)),
        " | Hits: ", sum(blast_filtered$Uniprot == "Uniprot Hit")
      ),
      x = "Identity %",
      y = "Number of proteins"
    ) +
    
    theme(
      panel.grid.minor = element_blank()
    )
  
  # SECOND LEVEL : PFAM
  if (!is.null(directory2)) {
    
    pfam_genes <- unique(pfam$V2)
    blast_filtered$PFAM <- ifelse(blast_filtered$Gene %in% pfam_genes,
                                 "Present", "Absent")
    
    dfB <- blast_filtered %>%
      count(Uniprot, PFAM)
    
    pB <- ggplot(dfB, aes(x = Uniprot, y = n, fill = PFAM)) +
      geom_col() +
      geom_text(
        aes(
          label = paste0("(", n, "/", total, ")"),
          color = PFAM
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
        fill = guide_legend(title = "PFAM Proteins"),
        color = "none"
      ) +
      theme_bw() +
      labs(
        title = "B - PFAM support across Helixer-specific proteins and UniProt mappings",
        subtitle = paste0(
        "QCover filtered > ", threshold_ident, "%",
        " | Total proteins: ", length(unique(blast_results$Gene)),
        " | Hits: ", sum(blast_filtered$Uniprot == "Uniprot Hit")
      ),
        x = "",
        y = ""
      )
  } else {
    pB <- ggplot() + theme_void() + labs(title = "B — PFAM not available")
  }
  
  
  # THIRD LEVEL : EggNOG
  if (!is.null(directory3)) {
    
    no_match <- df_sorted %>%
      filter(is.na(Evalue)) %>%
      select(Gene)
    
    no_match$EggNOG <- ifelse(no_match$Gene %in% giving_match$Gene,
                              "Recovered",
                              "Missing")
    
    dfC <- no_match %>% count(EggNOG)
    
    pC <- ggplot(dfC, aes(x = EggNOG, y = n, fill = EggNOG)) +
      geom_col() +
      geom_text(
        aes(label = paste0("(", n, "/", nrow(no_match), ")")),
        vjust = -0.4,
        size = 3
      ) +
      scale_fill_manual(values = c(
        "Recovered" = "#63003C",
        "Missing" = "#D6DBDF"
      )) +
      theme_bw() +
      labs(title = "C - EggNOG annotations for Uniprot unmatched proteins",
           x = "",
           y = "") +
      theme(
        legend.position = "none"
      )
  } else {
    pC <- ggplot() + theme_void() + labs(title = "C — EggNOG not available")
  }
  
  
  # Final figure ahh
  suppressMessages(
    suppressWarnings(
      final_fig <- grid.arrange(
        pA, pB, pC,
        ncol = 3
      )
    )
  )
  
  ggsave(file.path(save, "../figs/ProteinAnalysisEveryLevel.png"),
         final_fig,
         width = 20,
         height = 7,
         dpi = 600)
  
}
# Rscript Projet/Rapport/CodesFig/ProteinAnalysisSpecies.r TAIR12 Projet/Species/TAIR12 Projet/Species/TAIR12/Protein_Analysis.tsv Projet/Species/TAIR12/Arabidopsis_thaliana_TAIR12_pfamHelixer.txt 

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
save_dir <- args[2]
threshold_ident <- as.numeric(args[3])
input1 <- args[4]

input2 <- NULL
if (length(args) > 4) {
  input2 <- args[5]
}

input3 <- NULL
if (length(args) > 5) {
  input3 <- args[6]
}
result <- Protein_AnalysisSpecies(name, save_dir, threshold_ident, input1, input2, input3)