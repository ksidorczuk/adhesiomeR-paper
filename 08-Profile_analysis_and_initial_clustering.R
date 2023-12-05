library(rmarkdown)
library(dplyr)
library(cluster)
library(factoextra)
library(DT)
library(tidyr)
library(ggplot2)
library(adhesiomeR)
library(class)
library(patchwork)
library(htmltools)
library(ggdendro)

set.seed(73607254)

source("functions/process_results.R")
source("functions/cluster_plots.R")

if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Users/kar23yal/Dropbox/adhesiomeR/"
} else if(Sys.info()[["nodename"]] == "adhesiomer-testing") {
  data_path <- "~/hpc-home/adhesiomeR/"
} else {
  data_path <- "/hpc-home/kar23yal/adhesiomeR/"
}

colors <- adhesiomeR::pathotype_colors

# All data
gene_data <- read_in_gene_results(paste0(data_path, "Pathotypes/adhesiomeR_results_strict/"))

unique_gene_data <- unique(select(gene_data, -c(File, pathotype)))
x <- unique_gene_data[which(colSums(as.matrix(unique_gene_data)) != 0)]

dist_jaccard <- dist(x = as.matrix(x), method = "binary")
dist_man <- dist(x = as.matrix(x), method = "manhattan")

meta <- mutate(unique_gene_data, profile = 1:nrow(unique_gene_data)) %>% 
  left_join(gene_data) %>% 
  select(c("profile", "File", "pathotype"))

profiles_genomes_df <- meta %>%
  group_by(profile) %>% 
  summarise(genomes = n()) %>% 
  arrange(desc(genomes)) %>% 
  mutate(Adhesins_profile = paste0("A-", 1:nrow(.)))

mutate(unique_gene_data, profile = 1:nrow(unique_gene_data)) %>% 
  left_join(meta) %>% 
  left_join(profiles_genomes_df) %>% 
  select(-c("profile", "genomes")) %>% 
  write.csv(paste0(data_path, "Clustering/Profiles_all_genes.csv"), row.names = FALSE)

a_profile_df <- mutate(unique_gene_data, profile = 1:nrow(unique_gene_data)) %>% 
  left_join(meta) %>% 
  left_join(profiles_genomes_df) %>% 
  select(-c("profile", "genomes", "File", "pathotype")) %>%
  unique()
write.csv(a_profile_df, paste0(data_path, "Clustering/Profiles_table_all_genes.csv"), row.names = FALSE)

write.csv(select(profiles_genomes_df, -profile), paste0(data_path, "Clustering/Profiles_all_genes_summary.csv"), row.names = FALSE)



# Fimbrial adhesins
fimbrial_systems <- adhesiomeR::adhesins_df_grouped %>%
  group_by(System) %>%
  summarise(n = n()) %>% 
  filter(n > 1) %>% 
  .[["System"]]
fimbrial_genes <- unique(filter(adhesiomeR::adhesins_df_grouped, System %in% fimbrial_systems)[["Gene"]])

fimbrial_x <- unique(x[, which(colnames(x) %in% fimbrial_genes)])

dist_fimbrial <- dist(x = as.matrix(fimbrial_x), method = "binary")
dist_fimbrial_man <- dist(x = as.matrix(fimbrial_x), method = "manhattan")

fimbrial_meta <- mutate(fimbrial_x, profile = 1:nrow(fimbrial_x)) %>% 
  left_join(select(gene_data,  c(fimbrial_genes, "File", "pathotype"))) %>% 
  select(c("profile", "File", "pathotype"))

fimbrial_profiles_genomes_df <- fimbrial_meta %>%
  group_by(profile) %>% 
  summarise(genomes = n()) %>% 
  arrange(desc(genomes)) %>% 
  mutate(Fimbrial_profile = paste0("F-", 1:nrow(.)))

mutate(unique(fimbrial_x), profile = 1:nrow(unique(fimbrial_x))) %>% 
  left_join(select(gene_data, c(fimbrial_genes, "File", "pathotype"))) %>% 
  left_join(fimbrial_profiles_genomes_df) %>% 
  select(-c("profile", "genomes")) %>% 
  write.csv(paste0(data_path, "Clustering/Profiles_fimbrial_genes.csv"), row.names = FALSE)

f_profile_df <- mutate(unique(fimbrial_x), profile = 1:nrow(unique(fimbrial_x))) %>% 
  left_join(select(gene_data, c(fimbrial_genes, "File", "pathotype"))) %>% 
  left_join(fimbrial_profiles_genomes_df) %>% 
  select(-c("profile", "genomes", "File", "pathotype")) %>%
  unique() 
write.csv(f_profile_df, paste0(data_path, "Clustering/Profiles_table_fimbrial_genes.csv"), row.names = FALSE)

write.csv(select(fimbrial_profiles_genomes_df, -profile), paste0(data_path, "Clustering/Profiles_fimbrial_genes_summary.csv"), row.names = FALSE)


# Nonfimbrial adhesins
nonfimbrial_systems <- adhesiomeR::adhesins_df[["System"]][which(!(adhesiomeR::adhesins_df[["System"]] %in% fimbrial_systems))]
nonfimbrial_genes <- unique(filter(adhesiomeR::adhesins_df_grouped, System %in% nonfimbrial_systems)[["Gene"]])

nonfimbrial_x <- unique(x[, which(colnames(x) %in% nonfimbrial_genes)])

dist_nonfimbrial <- dist(x = as.matrix(nonfimbrial_x), method = "binary")
dist_nonfimbrial_man <- dist(x = as.matrix(nonfimbrial_x), method = "manhattan")

nonfimbrial_meta <- mutate(nonfimbrial_x, profile = 1:nrow(nonfimbrial_x)) %>% 
  left_join(gene_data) %>% 
  select(c("profile", "File", "pathotype"))

nonfimbrial_profiles_genomes_df <- nonfimbrial_meta %>%
  group_by(profile) %>% 
  summarise(genomes = n()) %>% 
  arrange(desc(genomes)) %>% 
  mutate(Nonfimbrial_profile = paste0("N-", 1:nrow(.)))

mutate(unique(nonfimbrial_x), profile = 1:nrow(unique(nonfimbrial_x))) %>% 
  left_join(select(gene_data, all_of(c(nonfimbrial_genes, "File", "pathotype")))) %>% 
  left_join(nonfimbrial_profiles_genomes_df) %>% 
  select(-c("profile", "genomes")) %>% 
  write.csv(paste0(data_path, "Clustering/Profiles_nonfimbrial_genes.csv"), row.names = FALSE)

n_profile_df <- mutate(unique(nonfimbrial_x), profile = 1:nrow(unique(nonfimbrial_x))) %>% 
  left_join(select(gene_data, c(nonfimbrial_genes, "File", "pathotype"))) %>% 
  left_join(nonfimbrial_profiles_genomes_df) %>% 
  select(-c("profile", "genomes", "File", "pathotype")) %>%
  unique() 
write.csv(n_profile_df, paste0(data_path, "Clustering/Profiles_table_nonfimbrial_genes.csv"), row.names = FALSE)

write.csv(select(nonfimbrial_profiles_genomes_df, -profile), paste0(data_path, "Clustering/Profiles_nonfimbrial_genes_summary.csv"), row.names = FALSE)

profiles <- list("A" = a_profile_df,
                 "F" = f_profile_df,
                 "N" = n_profile_df)
save(profiles, file = paste0(data_path, "adhesiomeR_files/profiles.rda"), compress = "xz")

# Save data

save(list = ls(), file = paste0(data_path, "Clustering/clustering.RData"))


render("Clara_clustering.Rmd", output_file = paste0("Clara_k=10_all_man.html"), output_dir = paste0(data_path, "Clustering/Tests/"),
       params = list(k = 10, type = "all_man"))


render("Clara_clustering.Rmd", output_file = paste0("Clara_k=8_fimbrial_man.html"), output_dir = paste0(data_path, "Clustering/Tests/"),
       params = list(k = 8, type = "fimbrial_man"))


render("Clara_clustering.Rmd", output_file = paste0("Clara_k=5_nonfimbrial_man.html"), output_dir = paste0(data_path, "Clustering/Tests/"),
       params = list(k = 5, type = "nonfimbrial_man"))

