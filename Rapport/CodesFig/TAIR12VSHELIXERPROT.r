library(Biostrings)
library(ggplot2)

fasta1 <- readAAStringSet("/home/crakshay/SUSHI/Work/GBOT/flagv5.0/Projet/Species/TAIR12/all_proteins_Arabidopsis_thaliana_TAIR12_CommonHELIXER.fasta")
fasta2 <- readAAStringSet("/home/crakshay/SUSHI/Work/GBOT/flagv5.0/Projet/Species/TAIR12/all_proteins_Arabidopsis_thaliana_TAIR12_CommonHELIXER.fasta")

identities <- mapply(function(s1, s2) {
  
  aln <- Biostrings::pairwiseAlignment(s1, s2, type = "global")
  
  pat <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
  sub <- strsplit(as.character(alignedSubject(aln)), "")[[1]]
  
  valid <- pat != "-" & sub != "-"
  
  sum(pat[valid] == sub[valid]) / sum(valid) * 100
  
}, fasta1, fasta2)

identity_class <- cut(
  identities,
  breaks = c(-Inf, 50, 70, 90, 99.999, 100),
  labels = c("<50%", "50–69%", "70–89%", "90–99%", "100%"),
  include.lowest = TRUE
)

df <- as.data.frame(table(identity_class))
colnames(df) <- c("Identity", "Count")

ggplot(df, aes(x = Identity, y = Count)) +
  geom_col(fill = "#2C7FB8") +
  geom_text(aes(label = Count), vjust = -0.3) +
  labs(
    title = "Protein Identity Distribution (Helixer vs Reference)",
    x = "Identity class",
    y = "Number of protein pairs"
  ) +
  theme_classic(base_size = 14)