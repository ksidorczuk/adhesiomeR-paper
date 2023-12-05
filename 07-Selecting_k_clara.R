library(dplyr)
library(cluster)
library(factoextra)
library(tidyr)
library(ggplot2)
library(adhesiomeR)
library(class)

set.seed(73607254)

read_in_gene_results <- function(adhesiomer_res_path, copies = FALSE) {
  lf <- if(copies == TRUE) {
    list.files(adhesiomer_res_path, pattern = "gene_copies", full.names = TRUE)
  } else {
    list.files(adhesiomer_res_path, pattern = "gene_presence", full.names = TRUE)
  }
  lapply(lf, function(ith_file) {
    read.csv(ith_file, check.names = FALSE) %>% 
      mutate(pathotype = strsplit(last(strsplit(ith_file, "/")[[1]]), "_")[[1]][1])
  }) %>% bind_rows()
}

if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Projects/adhesiomeR/"
} else {
  data_path <- "/hpc-home/kar23yal/adhesiomeR/"
}


colors <- adhesiomeR::pathotype_colors


gene_data <- read_in_gene_results(paste0(data_path, "Pathotypes/adhesiomeR_results_strict/"))

unique_gene_data <- unique(select(gene_data, -c(File, pathotype)))
x <- unique_gene_data[which(colSums(as.matrix(unique_gene_data)) != 0)]


dist_jaccard <- dist(x = as.matrix(x), method = "manhattan")


## Clustering - gene presence/absence
kmax = 200


clara_gap <- fviz_nbclust(x, clara, method = "gap_stat", k.max = kmax, diss = dist_jaccard, nboot = 10)
ggsave(paste0(data_path, "Clustering/SelectingK/clara_gapstat_man.png"), (clara_gap + theme(axis.text.x = element_text(angle = 90))), height = 6, width = 35)


save(list = ls(), file = paste0(data_path, "Clustering/SelectingK/clustering_data_clara_man.RData"))
