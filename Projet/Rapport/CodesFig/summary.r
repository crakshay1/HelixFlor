library(gt)
library(readr)

name <- "TAIR12" 
directory <- "~/SUSHI/Work/GBOT/flagv5.0/Projet/Species/TAIR12"


stats <- readr::read_delim(paste0(directory, "/../../Rapport/resume.txt"), delim = "\t")

first <- data.frame(
Annotation = c("Reference", "Helixer"),
CDS = c(stats[stats$Species == name, ]$NbExperts, stats[stats$Species == name, ]$NbHelixer),
mRNAs = c(stats[stats$Species == name, ]$NbExpertsRNA, stats[stats$Species == name, ]$NbHelixerRNA),
Isoforms = c(stats[stats$Species == name, ]$NbExpertsRNA > stats[stats$Species == name, ]$NbExperts,
             stats[stats$Species == name, ]$NbHelixerRNA > stats[stats$Species == name, ]$NbHelixer)
)

second <- data.frame(
  Added = stats[stats$Species == name, ]$Added,
  Missed = stats[stats$Species == name, ]$Missed,
  Shared = stats[stats$Species == name, ]$Common
)

file2 <- paste0(directory, "/UTRAnalysis/Arabidopsis_thaliana_UTR_ForAlmostSameCDS.tsv")
file <- paste0(directory, "/UTRAnalysis/Arabidopsis_thaliana_UTR_ForSameCDS.tsv")
same <- readr::read_delim(file, delim = "\t")
almost <- readr::read_delim(file2, delim = "\t")
same$length5h <- nchar(same$`5 UTR Helixer`)
same$length3h <- nchar(same$`3 UTR Helixer`)
almost$length5h <- nchar(almost$`5 UTR Helixer`)
almost$length3h <- nchar(almost$`3 UTR Helixer`)

same$length5e <- nchar(same$`5 UTR Expert`)
same$length3e <- nchar(same$`3 UTR Expert`)
almost$length5e <- nchar(almost$`5 UTR Expert`)
almost$length3e <- nchar(almost$`3 UTR Expert`)


third <- data.frame(
  CDS = c("Identical", "Almost Identical"),
  Number = c(length(same$Gène), length(almost$Gène)),
  Identity = c(mean(same$`Score Identité`), mean(almost$`Score Identité`)),
  #hUTR5 = c(median(same$length5h, na.rm = TRUE), median(almost$length5h, na.rm = TRUE)),
  #eUTR5 = c(median(same$length5e, na.rm = TRUE), median(almost$length5e, na.rm = TRUE)),
  #hUTR3 = c(median(same$length3h, na.rm = TRUE), median(almost$length3h, na.rm = TRUE)),
  #eUTR3 = c(median(same$length3e, na.rm = TRUE), median(almost$length3e, na.rm = TRUE)),
  type = c(names(which.max(table(unlist(strsplit(same$`Type Egalité`, ""))))),
           names(which.max(table(unlist(strsplit(almost$`Type Egalité`, "")))))),
  hLEN5 = c(mean(same$`5 UTR Helixer Len`, na.rm = TRUE), 
             mean(almost$`5 UTR Helixer Len`, na.rm = TRUE)),
  hLEN3 = c(mean(same$`3 UTR Helixer Len`, na.rm = TRUE), 
             mean(almost$`3 UTR Helixer Len`, na.rm = TRUE)),
  rLEN5 = c(mean(same$`5 UTR Expert Len`, na.rm = TRUE), 
             mean(almost$`5 UTR Expert Len`, na.rm = TRUE)),
  rLEN3 = c(mean(same$`3 UTR Expert Len`, na.rm = TRUE), 
             mean(almost$`3 UTR Expert Len`, na.rm = TRUE))
)

gt(first) %>%
    tab_header(
        title = paste0(name,": Genome annotations summary"),
        subtitle = "Comparative overview of Helixer and reference genome annotations"
    ) %>%
    fmt_number(
        columns = c(CDS, mRNAs),
        decimals = 0
    ) %>%
    tab_style(
        style = list(
        cell_text(weight = "bold")
        ), locations = cells_column_labels()
    ) %>%
    tab_options(
        table.font.size = px(12),
        heading.title.font.size = px(14),
        heading.subtitle.font.size = px(11)
    ) %>%
    cols_align(
        align = "center",
        columns = everything()
    )

gt(second) %>%
    tab_header(
        title = paste0(name,": Genome annotations summary"),
        subtitle = "CDS annotated by Helixer"
    ) %>%
    fmt_number(
        columns = c(Added, Missed, Shared),
        decimals = 0
    ) %>%
    tab_style(
        style = list(
        cell_text(weight = "bold")
        ), locations = cells_column_labels()) %>%
    tab_options(
        table.font.size = px(12),
        heading.title.font.size = px(14),
        heading.subtitle.font.size = px(11)
    ) %>%
    cols_align(
        align = "center",
        columns = everything()
    ) %>%
    tab_footnote(
      footnote = "Proteins generated from these CDS were analyzed.",
      locations = cells_body(columns = Added)
    )

gt(third) %>%
    tab_header(
        title = paste0(name,": Genome annotations summary"),
        subtitle = "mRNAs Structural analysis for specific CDS structures"
    ) %>%
    fmt_number(
        columns = c(Identity, hLEN5, hLEN3, rLEN5, rLEN3),
        decimals = 0
    ) %>%
    tab_style(
        style = list(
        cell_text(weight = "bold")
        ), locations = cells_column_labels()) %>%
    tab_options(
        table.font.size = px(12),
        heading.title.font.size = px(14),
        heading.subtitle.font.size = px(11)
    ) %>%
    cols_align(
        align = "center",
        columns = everything()
    ) %>%
    cols_label(
      CDS = "CDS Structure",
      Number = "Number",
      Identity = "Mean Identity %",
      type = "Most Frequent Type",
      hLEN5 = "5'UTR Helixer Length",
      #hUTR5 = "5'UTR Helixer ALD",
      hLEN3 = "3'UTR Helixer Length",
      #hUTR3 = "3'UTR Helixer ALD",
      rLEN5 = "5'UTR Reference Length",
      #eUTR5 = "5'UTR Reference ALD",
      rLEN3 = "3'UTR Reference Length",
      #eUTR3 = "3'UTR Reference ALD"
    ) %>%
    tab_footnote(
      footnote = "Only mRNAs linked to the CDS above were compared here.",
      locations = cells_body(columns = c(Identity, rLEN3, rLEN5, hLEN3, hLEN5))
    ) %>%
    tab_footnote(
        footnote = "Average Length of UTR regions",
        locations = cells_body(columns = c(rLEN5, rLEN3, hLEN3, hLEN5))
    )

if (name == "TAIR12") {
    VersusCDS <- readr::read_delim(paste0(directory, "/Analysis/TAIR12vsARAPORT11_CDS.txt"))
    VersusRNA <- readr::read_delim(paste0(directory, "/Analysis/TAIR12vsARAPORT11_mRNA.txt"))
    fourth <- data.frame(
        Number = nrow(VersusCDS),
        DiffLength = mean(abs(VersusCDS$TAIRLength - VersusCDS$ARALength)),
        DiffExons = mean(abs(VersusCDS$TAIRExons - VersusCDS$ARAExons))
    )
    valid_CDS <- VersusCDS[abs(VersusCDS$TAIRLength - VersusCDS$ARALength) <= 10 & 
                           abs(VersusCDS$TAIRExons - VersusCDS$ARAExons) <= 2, ]
    
    invalid_RNA <- VersusRNA[!VersusRNA$id_feat %in% valid_CDS$id_feat, ]
    
    fifth <- data.frame(
      Number     = nrow(invalid_RNA),
      DiffLength = mean(abs(invalid_RNA$TAIRLength - invalid_RNA$ARALength), na.rm = TRUE),
      DiffExons  = mean(abs(invalid_RNA$TAIRExons - invalid_RNA$ARAExons), na.rm = TRUE)   
    )
    
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    
    commons <- read_delim(
      paste0(directory, "/Stats/Commons_Arabidopsis_thaliana_TAIR12.txt"),
      delim = "\t"
    )
    commons1 <- read_delim(
      "~/SUSHI/Work/GBOT/flagv5.0/Projet/Species/ARAPORT11/Stats/Commons_Arabidopsis_thaliana.txt",
      delim = "\t"
    )
    df <- VersusCDS %>%
      left_join(
        commons %>%
          select(Reference, Length, NbExons) %>%
          rename(
            HelixerLength_TAIR = Length,
            HelixerExons_TAIR  = NbExons
          ),
        by = c("id_feat" = "Reference")
      ) %>%
      left_join(
        commons1 %>%
          select(Reference, Length, NbExons) %>%
          rename(
            HelixerLength_ARA = Length,
            HelixerExons_ARA  = NbExons
          ),
        by = c("id_feat" = "Reference")
      )
    
    
    df_plot <- df %>%
      mutate(
        Diff_Helixer_TAIR = abs(HelixerLength_TAIR - TAIRLength),
        Diff_Helixer_ARA  = abs(HelixerLength_ARA - ARALength)
      ) %>%
      select(Diff_Helixer_TAIR,
             Diff_Helixer_ARA) %>%
      pivot_longer(cols = everything(),
                   names_to = "Comparison",
                   values_to = "Value")
    
    df_plot$Comparison <- recode(df_plot$Comparison,
                                 Diff_Helixer_TAIR = "Helixer vs TAIR12",
                                 Diff_Helixer_ARA  = "Helixer vs ARAPORT11"
    )
    
    p4 <- ggplot(df_plot,
                 aes(x = Comparison,
                     y = Value,
                     fill = Comparison)) +
      scale_fill_manual(values = c(
        "Helixer vs ARAPORT11" = "#c1ff72",
        "Helixer vs TAIR12" = "#7ed957"
      )) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.15, outlier.alpha = 0.2) +
      
      scale_y_log10() +
      
      theme_bw() +
      
      labs(
        title = "Helixer structural agreement with reference annotations",
        x = "",
        y = "Absolute CDS length difference (log scale)"
      ) +
      
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      )
    
    ggsave(
      filename = paste0(directory, "/Figure4_Helixer_vs_refs.png"),
      plot = p4,
      width = 6.5,
      height = 5,
      dpi = 600
    )
    

    library(Biostrings)
    # Read the FASTA files
    
    fasta11 <- readAAStringSet(paste0(directory, "/all_proteins_Arabidopsis_thaliana_TAIR12_Versus_ARA.fasta"))
    fasta22 <- readAAStringSet(paste0(directory, "/all_proteins_Arabidopsis_thaliana_TAIR12_Versus_TAIR.fasta"))
    
    # Fonction : garde uniquement l'isoforme principale (suffixe .1)
    keep_primary_isoform <- function(fasta) {
      gene_id <- sub("\\.[0-9]+$", "", names(fasta))  # enlève .1, .2, .3, etc.
      fasta_dedup <- fasta[!duplicated(gene_id)]       # garde la première occurrence par gène
      names(fasta_dedup) <- gene_id[!duplicated(gene_id)]
      fasta_dedup
    }
    
    fasta1 <- keep_primary_isoform(fasta11)
    fasta2 <- keep_primary_isoform(fasta22)
    
    # Vérification
    length(fasta11)          # nombre total de transcrits avant
    length(fasta1)  # nombre de gènes (isoforme .1 uniquement)
    
    length(fasta22)
    length(fasta2)
    
    # Extract names (headers) and sequences
    headers1 <- names(fasta1)
    headers2 <- names(fasta2)
    
    seqs1 <- as.character(fasta1)
    seqs2 <- as.character(fasta2)
    
    # Create named vectors for easy lookup by header
    names(seqs1) <- headers1
    names(seqs2) <- headers2
    
    # Find common headers
    common_headers <- intersect(headers1, headers2)
    
    # Compare sequences for the headers that exist in both files
    if (length(common_headers) > 0) {
      # Extract sequences for common headers in the same order
      match_seq1 <- seqs1[common_headers]
      match_seq2 <- seqs2[common_headers]
      
      # Check where they are exactly equal
      same_sequence <- match_seq1 == match_seq2
      not_same_sequence <- !same_sequence
      num_same <- sum(same_sequence)
      num_diff <- sum(not_same_sequence)
    } else {
      num_same <- 0
      num_diff <- 0
    }
    
    I_decided <- FALSE
    if (I_decided) {
      library(ggplot2)
      # Select differing proteins
      diff_headers <- common_headers[not_same_sequence]
    
      # Build single FASTA input
      tmp_in <- tempfile(fileext = ".fa")
      tmp_out <- tempfile(fileext = ".aln.fa")
      
      seq_lines <- c()
      for (h in diff_headers) {
        seq_lines <- c(
          seq_lines,
          paste0(">", h, "_1"),
          match_seq1[h],
          paste0(">", h, "_2"),
          match_seq2[h]
        )
      }
      writeLines(seq_lines, tmp_in)
      
      # Run MAFFT once
      system(paste("mafft --auto", tmp_in, ">", tmp_out))
      
      # Read alignment
      aln <- readAAStringSet(tmp_out)
      aln_seq <- as.character(aln)
      
      # Compute identity per pair
      identity_vals <- sapply(diff_headers, function(h) {
        
        s1 <- as.character(aln_seq[[paste0(h, "_1")]])
        s2 <- as.character(aln_seq[[paste0(h, "_2")]])
        
        s1_vec <- strsplit(s1, "")[[1]]
        s2_vec <- strsplit(s2, "")[[1]]
        
        valid <- s1_vec != "-"
        
        100 * sum(s1_vec == s2_vec & valid) / sum(valid)
      })
      
      
      # Build dataframe
      identity_df <- data.frame(
        Header = diff_headers,
        Identity = as.numeric(identity_vals)
      )
      
      # Categorization
      identity_df$Category <- cut(
        identity_df$Identity,
        breaks = c(-Inf, 50, 70, 90, Inf),
        labels = c("< 50%", "50–70%", "70–90%", "> 90%")
      )
      
      identity_df$Category <- factor(
        identity_df$Category,
        levels = c("> 90%", "70–90%", "50–70%", "< 50%")
      )
      
      
      # Plot with ggplot
      ggplot(identity_df, aes(x = Category, fill = Category)) +
        geom_bar() +
        theme_bw() +
        scale_fill_manual(values = c(
          "> 90%" = "#00a3a6",
          "70–90%" = "#63003C"
        )) +
        labs(
          title = "Differing Protein identity distribution (MAFFT alignment)",
          x = "Identity% ranges",
          y = "Number of proteins"
        ) +
        theme(
          legend.position = "none"
        )
    }
    
    # Build the final summary table
    result_table <- data.frame(
      Metric = c("With Same Sequence", 
                 "With Different Sequence",
                 "Total"),
      Count = c(num_same, num_diff, num_diff+num_same)
    )
    print(
    gt(fourth) %>%
      tab_header(
          title = paste0(name,": Genome annotations summary"),
          subtitle = "TAIR12 VS ARAPORT11 : CDS"
      ) %>%
      fmt_number(
          columns = c(DiffExons),
          decimals = 1
      ) %>%
      tab_style(
          style = list(
          cell_text(weight = "bold")
          ), locations = cells_column_labels()) %>%
      tab_options(
          table.font.size = px(12),
          heading.title.font.size = px(14),
          heading.subtitle.font.size = px(11)
      ) %>%
      cols_align(
          align = "center",
          columns = everything()
      ) %>%
      cols_label(
        DiffExons = "Exons Number Difference",
        Number = "Number",
        DiffLength = "Length Difference (bp)",
      )
    )
  
    print(
    gt(fifth) %>%
      tab_header(
          title = paste0(name,": Genome annotations summary"),
          subtitle = "TAIR12 VS ARAPORT11 : mRNA for different CDS"
      ) %>%
      fmt_number(
          decimals = 0
      ) %>%
      tab_style(
          style = list(
          cell_text(weight = "bold")
          ), locations = cells_column_labels()) %>%
      tab_options(
          table.font.size = px(12),
          heading.title.font.size = px(14),
          heading.subtitle.font.size = px(11)
      ) %>%
      cols_align(
          align = "center",
          columns = everything()
      ) %>%
      cols_label(
        DiffExons = "Exons Number Difference",
        Number = "Number",
        DiffLength = "Length Difference (bp)",
      )
    )
    gt(result_table)%>%
      tab_header(
          title = paste0(name,": Genome annotations summary"),
          subtitle = "TAIR12 VS ARAPORT11 : Protein Comparison"
      ) %>%
      fmt_number(
          decimals = 0
      ) %>%
      tab_style(
          style = list(
          cell_text(weight = "bold")
          ), locations = cells_column_labels()) %>%
      tab_options(
          table.font.size = px(12),
          heading.title.font.size = px(14),
          heading.subtitle.font.size = px(11)
      ) %>%
      cols_align(
          align = "center",
          columns = everything()
      ) %>%
      cols_label(
          Metric = "Shared Proteins"
        )
}
