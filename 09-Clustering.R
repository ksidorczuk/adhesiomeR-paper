library(dplyr)
library(cluster)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(grid)
library(gridExtra)
library(adhesiomeR)
library(class)
library(patchwork)
library(ranger)
library(tidytext)

set.seed(73607254)

source("functions/process_results.R")
source("functions/cluster_plots.R")

get_size_df <- function(clustering_res, cluster_names) {
  data.frame(cluster = cluster_names,
             cl_cl[["clusinfo"]])
}


if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Users/kar23yal/Dropbox/adhesiomeR/"
} else {
  data_path <- "/hpc-home/kar23yal/adhesiomeR/"
}

load(paste0(data_path, "Clustering/clustering.RData"))


colors <- adhesiomeR::pathotype_colors

# Perform final clustering
clustering_all <- clara(x, k = 10, metric = "manhattan")
clustering_fimbrial <- clara(fimbrial_x, k = 8, metric = "manhattan")
clustering_nonfimbrial  <- clara(nonfimbrial_x, k = 5, metric = "manhattan")

# Combine clustering with data
all_clusters_df <- left_join(cbind(x, Cluster = clustering_all[["clustering"]]),
                             gene_data)
fimbrial_clusters_df <- left_join(cbind(fimbrial_x, Cluster = clustering_fimbrial[["clustering"]]),
                                  gene_data)
nonfimbrial_clusters_df <- left_join(cbind(nonfimbrial_x, Cluster = clustering_nonfimbrial[["clustering"]]),
                                     gene_data)

# Define cluster names
all_cluster_names <- c("A-G", "A-I", "A-E", "A-C", "A-A",
                       "A-D", "A-J", "A-B", "A-F", "A-H")
fimbrial_cluster_names <- c("F-F", "F-G", "F-A", "F-D", 
                            "F-C", "F-E", "F-H", "F-B")
nonfimbrial_cluster_names <- c("N-A", "N-D", "N-C", "N-B", "N-E")

# Summary with numbers of genomes in clusters
all_clusters_df %>% 
  group_by(Cluster) %>% 
  summarise(Genomes = n()) %>% 
  mutate(Cluster_name = all_cluster_names) %>% 
  write.csv(paste0(data_path, "Clustering/Clusters_all_genomes.csv"), row.names = FALSE)

fimbrial_clusters_df %>% 
  group_by(Cluster) %>% 
  summarise(Genomes = n()) %>% 
  mutate(Cluster_name = fimbrial_cluster_names) %>% 
  write.csv(paste0(data_path, "Clustering/Clusters_fimbrial_genomes.csv"), row.names = FALSE)

nonfimbrial_clusters_df %>% 
  group_by(Cluster) %>% 
  summarise(Genomes = n()) %>% 
  mutate(Cluster_name = nonfimbrial_cluster_names) %>% 
  write.csv(paste0(data_path, "Clustering/Clusters_nonfimbrial_genomes.csv"), row.names = FALSE)

# Change cluster names
clustering_all[["clustering"]] <- sapply(clustering_all[["clustering"]], function(i) all_cluster_names[i])
clustering_fimbrial[["clustering"]] <- sapply(clustering_fimbrial[["clustering"]], function(i) fimbrial_cluster_names[i])
clustering_nonfimbrial[["clustering"]] <- sapply(clustering_nonfimbrial[["clustering"]], function(i) nonfimbrial_cluster_names[i])

# Save clustering for use in R package
save(clustering_all, file = paste0(data_path, "adhesiomeR_files/clustering_all.rda"), compress = "xz")
save(clustering_fimbrial, file = paste0(data_path, "adhesiomeR_files/clustering_fimbrial.rda"), compress = "xz")
save(clustering_nonfimbrial, file = paste0(data_path, "adhesiomeR_files/clustering_nonfimbrial.rda"), compress = "xz")


### Publication plots ###
# Pathotype and cluster size plots
cl_path_plot_all <- get_pathotype_cluster_plot(all_clusters_df, clustering_all, all_cluster_names, 200)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_all.png"), cl_path_plot_all, width = 12, height = 6)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_all.eps"), cl_path_plot_all, width = 12, height = 6)
cl_path_plot_fimbrial <- get_pathotype_cluster_plot(fimbrial_clusters_df, clustering_fimbrial, fimbrial_cluster_names, 70)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_fimbrial.png"), cl_path_plot_fimbrial, width = 12, height = 6)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_fimbrial.eps"), cl_path_plot_fimbrial, width = 12, height = 6)
cl_path_plot_nonfimbrial <- get_pathotype_cluster_plot(nonfimbrial_clusters_df, clustering_nonfimbrial, nonfimbrial_cluster_names, 160)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_nonfimbrial.png"), cl_path_plot_nonfimbrial, width = 12, height = 5)
ggsave(paste0(data_path, "Clustering/Clusters_pathotypes_nonfimbrial.eps"), cl_path_plot_nonfimbrial, width = 12, height = 5)



imp_plot_all <- get_importance_plot(mutate(data.frame(clustering_all[["data"]], check.names = FALSE), 
                                           cluster = clustering_all[["clustering"]])) +
  ggtitle("All adhesin clusters")
imp_plot_fimbrial <- get_importance_plot(mutate(data.frame(clustering_fimbrial[["data"]], check.names = FALSE), 
                                                cluster = clustering_fimbrial[["clustering"]])) +
  ggtitle("Fimbrial adhesin clusters")
imp_plot_nonfimbrial <- get_importance_plot(mutate(data.frame(clustering_nonfimbrial[["data"]], check.names = FALSE), 
                                                   cluster = clustering_nonfimbrial[["clustering"]]), rows = 1) +
  ggtitle("Nonfimbrial adhesin clusters")

imp_combined <- wrap_plots(list(imp_plot_all, imp_plot_fimbrial, imp_plot_nonfimbrial), ncol = 1, heights = c(2,2,1))
ggsave(paste0(data_path, "Clustering/Clusters_gene_importance.png"), imp_combined, width = 12, height = 16)
ggsave(paste0(data_path, "Clustering/Clusters_gene_importance.eps"), imp_combined, width = 12, height = 16)
