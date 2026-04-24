library(ggplot2)

# Données comme dans ton graphe
data <- data.frame(
  species = rep(c("Arabidopsis Thaliana", "Brachypodium distachyon", "Malus domestica","Rosa chinensis", "Solanum lycopersicum"), each = 3),
  tool = rep(c("En commun", "Ajoutés par Helixer", "Ratés par Helixer"), times = 5),
  walltime = c(45925/(45925+1479+2341), 1479/(45925+1479+2341), 2341/(45925+1479+2341), 
               51203/(51203+3621+5632), 3621/(51203+3621+5632), 5632/(51203+3621+5632), 
               38005/(38005+5910+7110), 5910/(38005+5910+7110), 7110/(38005+5910+7110),
               32306/(32306+12196+7362), 12196/(32306+12196+7362), 7362/(32306+12196+7362),
               26397/(26397+3229+7677), 3229/(26397+3229+7677), 7677/(26397+3229+7677)
               )
)

# Barplot groupé
ggplot(data, aes(x = species, y = walltime*100, fill = tool)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    y = "% CDS",
    x = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = -0.1),
    plot.subtitle = element_text(face = "bold")
  ) +
  ggtitle("")






# Arabidopsis added
at_ad <- unique(Arabidopsis.thaliana_added$V1)
ad <- sapply(at_ad, function(element) {
  sum(Arabidopsis.thaliana_added$V1 == element)
})

# Arabidopsis missed
at_mi <- unique(Arabidopsis.thaliana_missed$V1)
mi <- sapply(at_mi, function(element) {
  sum(Arabidopsis.thaliana_missed$V1 == element)
})

# Malus added
ma_ad <- unique(Malus.domestica_added$V1)
ad1 <- sapply(ma_ad, function(element) {
  sum(Malus.domestica_added$V1 == element)
})

# Malus missed
ma_mi <- unique(Malus.domestica_missed$V1)
mi1 <- sapply(ma_mi, function(element) {
  sum(Malus.domestica_missed$V1 == element)
})

barplot(ad, col="darkblue", xlab = "Chromosome", ylab = "CDS ajoutés")
barplot(ad1, col="darkorange",, xlab = "Chromosome", ylab = "CDS ajoutés")
barplot(mi, col="darkred", , xlab = "Chromosome", ylab = "CDS ratés")
barplot(mi1, col="darkgreen", , xlab = "Chromosome", ylab = "CDS ratés")


length(MDsame[MDsame$Score.Identité >= 95,1])/length(MDsame$Score.Identité) * 100
length(ATsame[ATsame$Score.Identité >= 95,1]) / length(ATsame$Score.Identité) * 100

length(ATalmost[ATalmost$Score.Identité >= 95,1]) /length(ATalmost$Score.Identité) * 100
length(MDalmost[MDalmost$Score.Identité >= 95,1]) /length(MDalmost$Score.Identité) * 100