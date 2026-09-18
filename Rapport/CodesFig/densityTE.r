estimate_density <- function(dir) {
  suppressMessages(
    suppressWarnings(
      library(dplyr))
  )
  suppressMessages(
    suppressWarnings(
      library(ggplot2))
  )
  
  # Groupes à analyser
  groups <- c("added", "missed", "commonHEL", "commonEXP")
  
  # Stockage des résultats
  all_density <- list()
  
  for (group in groups) {
    windows <- read.table(
      paste0(dir, group, "_window.bed"),
      sep = "\t",
      header = FALSE
    )
    
    colnames(windows)[1:4] <- c("chr", "start", "end", "Gene")
    
    windows <- windows %>%
      mutate(
        Window_size = end - start
      )
    
    overlap <- read.table(
      paste0(dir, group, "_overlap.txt"),
      sep = "\t",
      header = FALSE
    )
    
    colnames(overlap)[4] <- "Gene"
    TE_per_gene <- overlap %>%
      group_by(Gene) %>%
      summarise(
        TE_bp = sum(V9),
        .groups = "drop"
      )
    
    density <- windows %>%
      select(Gene, Window_size) %>%
      left_join(
        TE_per_gene,
        by = "Gene"
      ) %>%
      mutate(
        TE_bp = ifelse(is.na(TE_bp), 0, TE_bp),
        Density = TE_bp / Window_size,
        Group = group
      )
    
    all_density[[group]] <- density
  }
  
  
  # Fusion des quatre groupes
  density_TE <- bind_rows(all_density)
  
  density_summary <- density_TE %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean_density = mean(Density),
      median_density = median(Density),
      sd_density = sd(Density)
    )
  
  write.csv(
    density_summary,
    file = paste0(dir, "../TES/", "density.csv"),
    row.names = FALSE
  )
  
  
  final_fig <- ggplot(
    density_TE,
    aes(
      x = Group,
      y = Density
    )
  ) +
    geom_boxplot() +
    ylab("TE density (TE bp / window size)") +
    xlab("Gene category") +
    theme_classic()
  
  
  kruskal.test(
    Density ~ Group,
    data = density_TE
  )
  
  ggsave(file.path(dir, paste0("../figs/DENSITY_TE.png")),
         final_fig,
         width = 10,
         height = 8.8,
         dpi = 600)
}

args <- commandArgs(trailingOnly = TRUE)

nom <- args[1]
result <- estimate_density(nom)

