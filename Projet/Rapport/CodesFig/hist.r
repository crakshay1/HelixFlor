# L'objectif ici est de pouvoir montrer graphiquement d'un point de vue général
# les différences entre l'annotation officielle et celle de Helixer.


library(ggplot2)
library(ggtext)
# Lecture des données
directory <- "~/SUSHI/Work/GBOT/flagv5.0"
resume <- read.delim(paste0(directory,'/Projet/Rapport/resume.txt'))
resume$Total <- resume$Missed + resume$Added + resume$Common

# Reshape
df <- data.frame(
  Species = rep(resume$Species, each = 3),
  Category = rep(c("Common", "Added", "Missed"), times = nrow(resume)),
  Value = as.vector(t(resume[, c("Common", "Added", "Missed")] / resume$Total))
)

# Renommer propre pour affichage
df$Category <- factor(df$Category,
                      levels = c("Common", "Added", "Missed"),
                      labels = c("In common", "Added by Helixer", "Missed by Helixer"))

# Tri
projects <- c(
  paste("Araly", "v1.0.12"),
  paste("ARAPORT", "11"),
  paste("Bd21", "v3.2"),
  paste("GDDH13", "v1.1"),
  paste("OBDH", "v1.0"),
  paste("Heinz1706", "ITAG4.0"),
  paste("TAIR", "12.1"),
  paste("NAM", "5.0.5")
)

species_levels <- unique(resume$Species[order(resume$Species)])

species_labels <- paste0(
  species_levels,
  "<br><span style='font-size:8pt; color:#666'>",
  projects,
  "</span>"
)

species_labels[species_levels == "ARAPORT 11"] <-
  paste0(
    "Arabidopsis thaliana",
    "<br><span style='font-size:8pt; color:#666'>ARAPORT 11</span>"
  )

species_labels[species_levels == "TAIR 12.1"] <-
  paste0(
    "Arabidopsis thaliana",
    "<br><span style='font-size:8pt; color:#666'>TAIR 12.1</span>"
  )

labels_sorted <- species_labels[order(
  sub("<br>.*", "", species_labels),
  species_labels
)]

species_level1 <- sub("<br>.*", "", labels_sorted)
species_level2 <- sub(".*<span.*>(.*)</span>.*", "\\1", labels_sorted)

species_key <- ifelse(
  grepl("Arabidopsis thaliana", species_level1),
  species_level2,
  species_level1
)

print(species_key)
df$Species <- factor(df$Species,
                     levels = species_key)
# Plot publication
p <- ggplot(df, aes(x = Species, y = Value * 100, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  scale_x_discrete(labels = labels_sorted) +
  
  # Valeurs au-dessus des barres
  geom_text(aes(label = paste0(round(Value*100,1) ,"%\n(", t(resume[, c("Common", "Added", "Missed")]),")")),
            lineheight = 0.8,
            fontface = "bold",
            family = "Courier New",
            size = 2.25,
            position = position_dodge(width = 0.75),
            vjust = -0.3) +
  
  scale_fill_manual(values = c(
    "In common" = "#00a3a6",
    "Missed by Helixer" = "#a2a5a8",
    "Added by Helixer" = "#63003C"
  )) +
  
  labs(
    y = "Number of Genes (%)",
    x = NULL,
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10, family = "Courier New"),
    axis.text.x = element_markdown(angle = 30, hjust = 1, family = "Courier New"),
    axis.title.y = element_text(face = "bold", family = "Courier New"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.6)
  )

# Export
ggsave(paste0(directory,"/Projet/Rapport/Figures/figure_CDS_comparison.png"),
       plot = p,
       width = 15,
       height = 5,
       dpi = 600)
