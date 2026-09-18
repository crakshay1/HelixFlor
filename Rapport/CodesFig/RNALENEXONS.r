library(gridExtra)
library(ggplot2)
name <- "TAIR12" 
directory <- "~/SUSHI/Work/GBOT/flagv5.0/Projet/Species/TAIR12"
mRNAeX <- readr::read_delim(paste0(directory, "/ExpertsLENEXONS.txt"), delim = "\t")
mRNAHe <- readr::read_delim(paste0(directory, "/HelixerLENEXONS.txt"), delim = "\t")
mRNAeX$Type <- name
mRNAHe$Type <-"Helixer"
overall_stats <- data.frame(
                            "Type" = c("Helixer", name),
                            "Length" = c(median(mRNAHe$Length), median(mRNAeX$Length)),
                            "Exons" = c(median(mRNAHe$NbExons), median(mRNAeX$NbExons)),
                            "5UTR" = c(median(mRNAHe$UTR5), median(mRNAeX$UTR5)),
                            "3UTR" = c(median(mRNAHe$UTR3), median(mRNAeX$UTR3))
                            )
overall_BoxStats <- data.frame(
    "Type" = c(mRNAHe$Type, mRNAeX$Type),
    "Id_feat" = c(mRNAHe$id_feat, mRNAeX$id_feat),
    "Length" = c(mRNAHe$Length, mRNAeX$Length),
    "Exons" = c(mRNAHe$NbExons, mRNAeX$NbExons),
    "UTR5" = c(mRNAHe$UTR5, mRNAeX$UTR5),
    "UTR3" = c(mRNAHe$UTR3, mRNAeX$UTR3)
)
p1<-ggplot(overall_BoxStats, aes(x=Type, y=Length, fill = Type)) +
  geom_boxplot(outlier.shape = NA, staplewidth = 0.2)+
  coord_cartesian(ylim = c(0, 5000)) +
  labs(
    x = NULL,
    y = "Length (bp)",
    title = "A - Length of mRNA"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p2<-ggplot(overall_BoxStats, aes(x=Type, y=Exons, fill = Type)) +
  geom_boxplot(outlier.shape = NA, staplewidth = 0.2)+
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    x = NULL,
    y = "Number of Exons",
    title = "B - Number of exons in mRNA"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p3<-ggplot(overall_BoxStats, aes(x=Type, y=UTR5, fill = Type)) +
  geom_boxplot(outlier.shape = NA, staplewidth = 0.2)+
  coord_cartesian(ylim = c(0, 600)) +
  labs(
    x = NULL,
    y = "5UTR Length (bp)",
    title = "C - 5'UTR length in mRNA"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

p4<-ggplot(overall_BoxStats, aes(x=Type, y=UTR3, fill = Type)) +
  geom_boxplot(outlier.shape = NA, staplewidth = 0.2)+
  coord_cartesian(ylim = c(0, 600)) +
  labs(
    x = NULL,
    y = "3UTR Length (bp)",
    title = "D - 3'UTR length in mRNA"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

final_fig <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = paste0("Differences in mRNA between Helixer & ", name, " annotations")
)

final_fig

overall_stats |>
  gt() |>
  cols_label(
    Length = "Length (bp)",
    Exons = "Number of Exons",
    X5UTR = "5'UTR (bp)",
    X3UTR = "3'UTR (bp)"
  ) |>
  fmt_number(
    columns = -Type,
    decimals = 0,
  ) |>
  tab_header(
    title = paste0("Differences in mRNA between Helixer & ", name, " annotations"),
    subtitle = "Median for all metrics and without isoforms"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  cols_align(
    align = "center",
    columns = everything()
    )|>
  opt_stylize(style = 1)