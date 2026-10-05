# L'objectif ici est de retrouver tous les TEs prédits par Helixer qui ne sont pas sensés être des TEs...
# (Oui j'suis revenu au français)

directory <- "~/SUSHI/Work/GBOT/flagv5.0/Projet/Species/B73/"
input_file  <- paste0(directory,"found_TEs.tsv")  # Sacré BlastN
output_te    <- paste0(directory,"/realHelixer_TEs.tsv") # Gènes identifiés comme TE 


data <- read_delim(input_file, delim = "\t", col_types = cols())

# - QCover >= 90% 
# - Evalue <= 1e-10
# - PIdent >= 60%
# Je voulais faire Chr==TEChr mais les ET sont mobiles...

te_false_positives <- data %>%
  filter(
    QCover >= 90,
    Evalue <= 1e-10,
    PIdent >= 60
  ) %>%
  # Comme ça les meilleurs hits (proches de 0) se retrouvent en haut
  arrange(Evalue) %>%
  # On garde uniquement le meilleur hit par gène
  distinct(Gene, .keep_all = TRUE)

first <- data.frame(
  Helixer = c("TEs", "No TEs"),
  Number = c(nrow(te_false_positives), 7267-nrow(te_false_positives))
)


# On prend tous les gènes du fichier de départ QUI NE SONT PAS dans notre liste de faux positifs
gt(first) %>%
  tab_header(
    title = "Detection of transposable elements within Helixer-specific genes",
    subtitle = "Results obtained after BlastN against MAIZEGDB's TEs Database"
  ) %>%
  fmt_number(
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
  ) %>%
  tab_footnote(
    footnote = "Qcover ≥90% / Evalue ≤1e-10 / Identity ≥90%",
  )

write_delim(te_false_positives, output_te, delim = "\t")