# This script was originally used with older version of adhesiomeR (commit 0f675a9).
# Results obtained with newer version may differ

library(adhesiomeR)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggplot2)
library(plotly)
library(pals)
library(biogram)

source("functions/process_results.R")
source("functions/colocalization_functions.R")

if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/adhesiomeR/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/adhesiomeR/"
}

if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Projects/adhesiomeR/"
}

all_blast_res <- read_in_blast_results(paste0(data_path, "Pathotypes/adhesiomeR_results/"))
adhesins_lengths <- adhesiomeR::adhesins_lengths
gene_groups <- adhesiomeR::gene_groups

blast_res_lens <- left_join(all_blast_res, adhesins_lengths, by = c("Subject" = "Gene"))
blast_res_lens[["Subject coverage"]] <- blast_res_lens[["Alignment length"]]/blast_res_lens[["Length"]]*100
all_res <- filter(blast_res_lens, `% identity` > 80 & `Subject coverage` > 80)

multigene_systems <- adhesins_df %>% 
  group_by(System) %>% 
  summarise(n = n()) %>% 
  filter(n > 2) %>% 
  .[["System"]]

system_results <- read_in_system_results(paste0(data_path, "Pathotypes/adhesiomeR_results/")) %>% 
  pivot_longer(2:(ncol(.)-1), names_to = "System", values_to = "Presence") %>% 
  filter(Presence == "Present") %>% 
  left_join(adhesiomeR::adhesins_df)

# Systems that have too many full operons on a single contig and will need
# addition of operons on multiple contigs
problematic_systems <- c("AA/II", "AA/III", "Afa-III", "CS14", "CS28B", "F1845", "Sfp", "P_2", "F4/K88")

### Full operons on the same contig or multiple contigs in the same genome
lapply(multigene_systems, function(ith_system) {
  print(paste0(ith_system))
  genes <- filter(adhesins_df, System == ith_system)[["Gene"]]
  x <- left_join(filter(system_results, Gene %in% genes),
                 filter(all_res, Subject %in% genes),
                 by = c("File", "pathotype", "Gene" = "Subject")) %>% 
    left_join(adhesins_lengths, by = "Gene")
  dat_contigs <- get_operon_fullnes_info(x, genes, type = "contigs")
  dat_genomes <- get_operon_fullnes_info(x, genes, type = "genomes")
  dat_contigs_summ <- dat_contigs %>% 
    group_by(pathotype, full) %>% 
    summarise(count = n()) 
  dat_genomes_summ <- dat_genomes %>% 
    group_by(pathotype, full) %>% 
    summarise(count = n()) 
  plot_dat <- bind_rows(mutate(dat_contigs_summ, type = "Contig"), 
                        mutate(dat_genomes_summ, type = "Genome"))
  if(nrow(plot_dat) > 0) {
    p <- ggplot(plot_dat, aes(x = count, y = pathotype, fill = full)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = count), size = 3, position = position_dodge(width = 0.9), hjust = -0.05) +
      ggtitle(ith_system) +
      theme_bw() +
      facet_wrap(~type) +
      scale_fill_manual("Full operon", values = c("TRUE" = "#e42b24", "FALSE" = "#85c1ff")) +
      theme(legend.position = "bottom")
    ggsave(paste0(data_path, "Gene_localization_operons/Contigs/", gsub("/", "-", ith_system), ".png"), p, width = 12, height = 6)
  }
  write.csv(dat_contigs, paste0(data_path, "Gene_localization_operons/Contigs/", gsub("/", "-", ith_system), "_contigs.csv"), row.names = FALSE)
  write.csv(dat_genomes, paste0(data_path, "Gene_localization_operons/Contigs/", gsub("/", "-", ith_system), "_genomes.csv"), row.names = FALSE)
  
  # Plot colocalization
  plot_data_contigs <- get_colocalization_data(dat_contigs, x, type = "contigs")
  plot_data_genomes <- get_colocalization_data(dat_genomes, x, type = "genomes")
  
  lapply(unique(plot_data_contigs[["pathotype"]]), function(ith_pathotype) {
    plot_subset_c <- filter(plot_data_contigs, pathotype == ith_pathotype)
    pc <- plot_operon_colocalization(plot_subset_c, pathotype = ith_pathotype, system = ith_system)
    if(nrow(plot_subset_c) > 5000) {
      htmlwidgets::saveWidget(ggplotly(pc), paste0(data_path, "Gene_localization_operons/Colocalization/Full_contigs_", gsub("/", "-", ith_system), "_", ith_pathotype, ".html"), selfcontained = TRUE)
    } else {
      ggsave(paste0(data_path, "Gene_localization_operons/Colocalization/Full_contigs_", gsub("/", "-", ith_system), "_", ith_pathotype, ".png"), 
             pc, width = 3+length(genes)/3, height = 3+nrow(plot_subset_c)/65, limitsize = FALSE)
    }
    plot_subset_g <- filter(plot_data_genomes, pathotype == ith_pathotype)
    pg <- plot_operon_colocalization(plot_subset_g, pathotype = ith_pathotype, system = ith_system)
    if(nrow(plot_subset_g) > 5000) {
      htmlwidgets::saveWidget(ggplotly(pg), paste0(data_path, "Gene_localization_operons/Colocalization/Full_genomes_", gsub("/", "-", ith_system), "_", ith_pathotype, ".html"), selfcontained = TRUE)
    } else {
      ggsave(paste0(data_path, "Gene_localization_operons/Colocalization/Full_genomes_", gsub("/", "-", ith_system), "_", ith_pathotype, ".png"), 
             pg, width = 5+length(genes)/3, height = 3+nrow(plot_subset_g)/65, limitsize = FALSE)
    }
  })
  ref_data <- mutate(plot_data_contigs, Type = "Contig", `Operon length` = end-start)
  write.csv(ref_data, paste0(data_path, "Gene_localization_operons/Reference_sequences/", gsub("/", "-", ith_system), ".csv"), row.names = FALSE)
  if(ith_system %in% problematic_systems) {
    ref_data2 <- mutate(plot_data_genomes, Type = "Genome") %>% 
      filter(!(File %in% unique(plot_data_contigs[["File"]])))
    write.csv(ref_data2, paste0(data_path, "Gene_localization_operons/Reference_sequences/", gsub("/", "-", ith_system), "_2.csv"), row.names = FALSE)
  }
})

# Example on type 1 fimbriae:
#---
genes <- filter(adhesins_df, System == "Type_1")[["Gene"]]
x <- filter(all_res, Subject %in% genes)
plot_dat <- x %>% 
  group_by(pathotype, File, Query) %>% 
  summarise(n_found = n(),
            n_unique = n_distinct(Subject),
            full = ifelse(n_unique == length(genes), TRUE, FALSE))
#---

# Check if similar system may have hits to the same contigs:
contig_dat_files <- list.files(paste0(data_path, "Gene_localization_operons/Contigs/"), 
                               pattern = "_contigs.csv",
                               full.names = TRUE)
all_contigs <- lapply(contig_dat_files, function(ith_file) {
  x <- read.csv(ith_file) 
  if(nrow(x) > 1) {
    mutate(x, system = gsub("_contigs.csv", "", last(strsplit(ith_file, "/")[[1]])))
  }
}) %>% bind_rows()

contigs_with_many_systems <- all_contigs %>% 
  filter(full == TRUE) %>% 
  group_by(Query) %>% 
  summarise(n = n(),
            systems = paste0(system, collapse = ", "))



