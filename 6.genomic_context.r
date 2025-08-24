library(ggalluvial)
library(dplyr)
library(readxl)
library(grid)

if (!file.exists("arg_data.xlsx") || !file.exists("metadata.xlsx")) {
  message("...")
  quit(save = "no", status = 1) 
}



arg_data <- read_excel("arg_data.xlsx")
metadata  <- read_excel("metadata.xlsx")

depleted_args <- c(
  "ESCHERICHIA COLI EF-TU MUTANTS CONFERRING RESISTANCE TO PULVOMYCIN",
  "VANS","ERMB","MEFA","VANU","MEL","EFRA","EFRB","TOLC",
  "BCRA","TET34","MACB","TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR",
  "MPHE","LSAE","LNUB","DNA-BINDING_PROTEIN_H-NS",
  "LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE","TRIC","OXA-333"
)

enriched_args <- c(
  "CCRA","BRUCELLA SUIS MPRF",
  "STAPHYLOCOCCUS MUPA CONFERRING RESISTANCE TO MUPIROCIN"
)

all_args_of_interest <- c(depleted_args, enriched_args)

#simplifies long  names
arg_label_map <- c(
  "ESCHERICHIA COLI EF-TU MUTANTS CONFERRING RESISTANCE TO PULVOMYCIN" = "Ecol_EFTu_PLV",
  "VANS"  = "vanS",  "ERMB"  = "ermB",  "MEFA" = "mefA",
  "VANU"  = "vanU",  "MEL"   = "mel",   "EFRA" = "efrA",
  "EFRB"  = "efrB",  "TOLC"  = "tolC",  "BCRA" = "bcrA",
  "TET34" = "tet34", "MACB"  = "macB",  "TRUNCATED_PUTATIVE_RESPONSE_REGULATOR_ARLR" = "arlR",
  "MPHE"  = "mphE",  "LSAE"  = "lsaE",  "LNUB" = "lnuB",
  "DNA-BINDING_PROTEIN_H-NS" = "H-NS",
  "LLMA_23S_RIBOSOMAL_RNA_METHYLTRANSFERASE" = "LlmA_23S_CLI",
  "TRIC"  = "triC",  "OXA-333" = "OXA-333",
  "CCRA"  = "ccrA",  "BRUCELLA SUIS MPRF" = "mprF",
  "STAPHYLOCOCCUS MUPA CONFERRING RESISTANCE TO MUPIROCIN" = "mupA"
)

group_MGE <- function(pred) {
  if (is.na(pred)) return(NA)
  clean_pred <- tolower(trimws(as.character(pred)))
  
  if (grepl("unclassified", clean_pred)) {
    return(NA)
  } else if (clean_pred == "chromosome") {
    return("Chromosome")
  } else if (clean_pred %in% c("phage","plasmid","ambiguous (plasmid/phage)")) {
    return("MGE")
  } else if (clean_pred == "ambiguous (phage/chromosome)") {
    return("Ambiguous Phage/Chromosome")
  } else {
    return("Other")
  }
}

arg_data <- arg_data %>%
  mutate(
    sample_id   = as.character(sample_id),
    ARG_label   = as.character(arg_label_map[ARG]),
    MGE_grouped = sapply(MGE_prediction, group_MGE)
  )


plot_prep <- arg_data %>%
  inner_join(metadata, by = "sample_id") %>%
  filter(group == "influent", ARG %in% all_args_of_interest) %>%
  count(ARG_label, MGE_grouped) %>%
  filter(!is.na(MGE_grouped) & MGE_grouped != "Other") %>%
  mutate(MGE_grouped = if_else(MGE_grouped == "Ambiguous Phage/Chromosome",
                               "Ambiguous",
                               MGE_grouped))

# rare args
rare_cutoff <- 3
rare_ARGs <- plot_prep %>%
  group_by(ARG_label) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  filter(total_n < rare_cutoff) %>%
  pull(ARG_label)

plot_prep <- plot_prep %>%
  mutate(ARG_grouped = if_else(ARG_label %in% rare_ARGs, "Other ARGs", ARG_label)) %>%
  group_by(ARG_grouped, MGE_grouped) %>%
  summarise(n = sum(n), .groups = "drop")


plot_prep <- plot_prep %>%
  mutate(Enrichment = if_else(ARG_grouped %in% c("ccrA","mprF","mupA"), "Enriched", "Depleted"))

plot_prep$MGE_grouped <- droplevels(as.factor(plot_prep$MGE_grouped))
plot_prep$Enrichment  <- factor(plot_prep$Enrichment, levels = c("Depleted","Enriched"))

alluvial_plot <- ggplot(plot_prep,
                        aes(axis1 = ARG_grouped,
                            axis2 = Enrichment,
                            axis3 = MGE_grouped,
                            y     = n)) +
  geom_alluvium(aes(fill = MGE_grouped),
                width = 0.25,
                color = "black",
                size = 0.3) +
  geom_stratum(width = 0.25, fill = "grey90", color = "black", size = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5, check_overlap = TRUE) +
  scale_x_discrete(
    limits = c("ARG_grouped","Enrichment","MGE_grouped"),
    labels = c("ARG","Enrichment Status","Genomic Context"),
    expand = c(0.1,0.2)
  ) +
  scale_fill_manual(
    name = "Genomic Context",
    values = c("Chromosome"="#008080","MGE"="#FF4500","Ambiguous"="#8000FF")
  ) +
  labs(x = "Category", y = "Number of Influent ARGs") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    panel.border         = element_rect(fill=NA,color="black",size=0.3),
    axis.line            = element_line(color="black",size=0.3),
    axis.ticks           = element_line(size=0.3),
    axis.title.x         = element_text(size=12,face="bold"),
    axis.title.y         = element_text(size=12,face="bold"),
    axis.text.x          = element_text(angle=45,hjust=1,size=8,color="black"),
    axis.text.y          = element_text(size=8),
    panel.grid.major.y   = element_line(linetype="dashed",size=0.3,color="grey70"),
    panel.grid.major.x   = element_blank(),
    panel.grid.minor     = element_blank(),
    legend.position      = c(1.02,0.5),
    legend.justification = c(0,0.5),
    legend.background    = element_blank(),
    legend.key           = element_rect(fill=NA,color=NA),
    legend.title         = element_text(size=8),
    legend.text          = element_text(size=8),
    legend.key.size      = unit(0.8,"lines"),
    legend.box.margin    = margin(l=5,unit="pt"),
    plot.margin          = margin(t=5,r=100,b=5,l=5,unit="pt")
  )

ggsave("genomic_contexts_ARGs.pdf",
       plot   = alluvial_plot,
       width  = 9,
       height = 6,
       units  = "in",
