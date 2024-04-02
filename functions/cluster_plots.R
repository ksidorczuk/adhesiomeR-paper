
my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")

options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

get_pathotype_plot <- function(cl_path) {
  cl_path %>% 
    ggplot(aes(x = cluster, y = n, fill = pathotype)) +
    geom_col(position = "fill") +
    theme_bw() +
    scale_fill_manual("Pathotype", values = colors) +
    coord_flip() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = size_df[["cluster"]])
}

get_size_plot <- function(size_df) {
  size_df %>% 
    ggplot(aes(x = cluster, y = size, label = size)) +
    geom_col() +
    theme_bw() +
    geom_text(nudge_y = 50) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
}


get_dendro_path_plot <- function(x_cl, cl_cl_path, size_df) {
  means_df <- lapply(unique(x_cl[["cluster"]]), function(i) {
    x_cl %>% 
      filter(cluster == i) %>% 
      select(-c("cluster", "genome")) %>% 
      colMeans()
  }) %>% bind_rows() %>% 
    mutate(cluster = unique(x_cl[["cluster"]]))
  
  clustering_clusters <- as.dendrogram(hclust(dist(as.matrix(select(means_df, -"cluster"))), method = "ward.D2"))
  cluster_order <- order.dendrogram(clustering_clusters)
  
  dendro <- clustering_clusters %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_clusters))) + 
    theme_void() + 
    coord_flip() + 
    scale_y_reverse()
  
  cl_cl_path2 <- cl_cl_path %>% 
    mutate(cluster = factor(cluster, levels = unique(cl_cl_path[["cluster"]])[cluster_order]))
  
  path_plot <- cl_cl_path2 %>% 
    ggplot(aes(x = cluster, y = n, fill = pathotype)) +
    geom_col(position = "fill") +
    theme_bw() +
    scale_fill_manual("Pathotype", values = colors) +
    coord_flip() +
    theme(legend.position = "bottom",
          axis.title.y = element_blank())
  
  size_plot <- size_df %>% 
    mutate(cluster = factor(cluster, levels = unique(cl_cl_path[["cluster"]])[cluster_order])) %>% 
    get_size_plot()
  
  wrap_plots(list(dendro, path_plot, size_plot), nrow = 1, widths = c(0.1, 0.6, 0.3), guides = "collect") &
    theme(legend.position = "bottom")
}

get_dendro <- function(x_cl) {
  means_df <- lapply(unique(x_cl[["cluster"]]), function(i) {
    x_cl %>% 
      filter(cluster == i) %>% 
      select(-c("cluster", "genome")) %>% 
      colMeans()
  }) %>% bind_rows() %>% 
    mutate(cluster = unique(x_cl[["cluster"]]))
  
  clustering_clusters <- as.dendrogram(hclust(dist(as.matrix(select(means_df, -"cluster"))), method = "ward.D2"))
  cluster_order <- order.dendrogram(clustering_clusters)
  
  dendro <- clustering_clusters %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("",
                       expand = c(0, 0)) +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_clusters)),
                     labels = dendro_data(clustering_clusters)[["labels"]][["label"]]) + 
    theme_minimal() +
    theme(axis.line = element_blank(), axis.text.y = element_blank(),
          panel.grid = element_blank())
  return(list("dendro" = dendro,
              "cl_order" = cluster_order))
}

get_pathotype_cluster_plot <- function(clusters_df, clustering, cluster_names, offset, type = "pathotypes", colors = "all") {
  if(colors == "all") {
    colors <- adhesiomeR::pathotype_colors
  }
  colors <- c("#b74786", "#d490b4", "#d73127", "#e7827b", "#ab615c",
              "#fc8d61", "#cda190", "#fcb598", "#ffff7d", "#ffffc9",
              "#e6e6b0", "#a6cee3", "#949494")
  cols_to_remove <- colnames(clusters_df)[which(!(colnames(clusters_df) %in% colnames(clustering[["data"]])))]
  x_cl <- clusters_df %>% 
    mutate(cluster = cluster_names[Cluster]) %>% 
    select(-all_of(cols_to_remove)) %>% 
    unique() %>% 
    mutate(genome = 1:nrow(.))
  
  means_df <- lapply(unique(x_cl[["cluster"]]), function(i) {
    x_cl %>% 
      filter(cluster == i) %>% 
      select(-c("cluster", "genome")) %>% 
      colMeans() %>%  
      t() %>% 
      data.frame(check.names = FALSE) %>% 
      mutate(cluster = i)
  }) %>% bind_rows()
  
  clustering_clusters <- as.dendrogram(hclust(dist(as.matrix(select(means_df, -"cluster"))), method = "ward.D2"))
  cluster_order <- order.dendrogram(clustering_clusters)
  genomes_df <- clusters_df %>% 
    mutate(cluster = cluster_names[Cluster]) %>% 
    group_by(cluster) %>% 
    summarise(genomes = n())
  size_df <- data.frame(cluster = cluster_names, clustering[["clusinfo"]]) %>% 
    left_join(genomes_df, by = "cluster")
  
  dendro <- clustering_clusters %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_clusters))) + 
    theme_void() + 
    coord_flip() + 
    scale_y_reverse()
  
  cl_cl_path <- summarise(group_by(mutate(clusters_df, cluster = cluster_names[Cluster]), cluster, pathotype), n = n())
  cl_cl_path[["pathotype"]] <- factor(cl_cl_path[["pathotype"]], levels = c("aEPEC", "tEPEC", "EAEC", "EHEC", "EIEC",
                                                                            "ETEC", "DAEC", "STEC", "UPEC", "APEC",
                                                                            "NMEC", "Nonpathogenic", "NA"))
  
  path_plot <- if(type == "pathotypes") {
    cl_cl_path %>% 
      ggplot(aes(x = cluster, y = n, fill = pathotype)) +
      geom_col(position = "fill") +
      theme_bw(base_size = 16) +
      scale_fill_manual("Pathotype", values = colors) +
      coord_flip() +
      theme(legend.position = "bottom",
            axis.title.y = element_blank()) +
      ylab("Pathotype fraction")
  } else {
    cl_cl_path %>% 
      mutate(pathotype = case_when(pathotype %in% c("APEC", "NMEC", "UPEC") ~ "InPEC",
                                   pathotype %in% c("aEPEC", "EAEC", "ETEC", "EHEC", 
                                                    "STEc", "DAEC", "EIEC", "tEPEC") ~ "ExPEC",
                                   pathotype == NA ~ "NA",
                                   pathotype == "Nonpathogenic" ~ "Nonpathogenic"))  %>% 
      ggplot(aes(x = cluster, y = n, fill = pathotype)) +
      geom_col(position = "fill") +
      theme_bw(base_size = 16) +
      scale_fill_manual("Pathotype", values = c("InPEC" = "#FFFF99", "ExPEC" = "#d73027", 
                                                "Nonpathogenic" = "#A6CEE3", `NA` = "#949494")) +
      coord_flip() +
      theme(legend.position = "bottom",
            axis.title.y = element_blank()) +
      ylab("Pathotype fraction")
  }
  
  size_plot <- size_df %>% 
    select(-all_of(c("max_diss", "av_diss", "isolation"))) %>% 
    pivot_longer(c("size", "genomes"), names_to = "type", values_to = "value") %>% 
    mutate(type = ifelse(type == "size", "Profiles", "Genomes")) %>% 
    ggplot(aes(x = cluster, y = value, label = value, fill = type)) +
    geom_col(position = "dodge") +
    theme_bw(base_size = 16) +
    geom_text(hjust = -0.05, size = 4, position = position_dodge(width = .9)) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ylab("Count") +
    ylim(c(0, max(size_df[["genomes"]]+offset))) +
    scale_fill_manual("Number of", values = c("Genomes" = "grey30", "Profiles" = "grey70"))
  
  wrap_plots(list(dendro, path_plot, size_plot), nrow = 1, widths = c(0.03, 0.67, 0.3)) &
    theme(legend.position = "bottom")
}


get_cluster_pathotype_plot <- function(clusters_df, clustering, cluster_names) {
  colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
  
  cols_to_remove <- colnames(clusters_df)[which(!(colnames(clusters_df) %in% colnames(clustering[["data"]])))]
  x_cl <- clusters_df %>% 
    mutate(cluster = cluster_names[Cluster]) %>% 
    select(-all_of(cols_to_remove)) %>% 
    unique() %>% 
    mutate(genome = 1:nrow(.))
  
  cl_cl_path <- summarise(group_by(mutate(clusters_df, cluster = cluster_names[Cluster]), cluster, pathotype), n = n())
  cl_cl_path[["pathotype"]] <- factor(cl_cl_path[["pathotype"]], levels = c("aEPEC", "tEPEC", "EAEC", "EHEC", "EIEC",
                                                                            "ETEC", "DAEC", "STEC", "UPEC", "APEC",
                                                                            "NMEC", "Nonpathogenic", "NA"))
  x <- cl_cl_path %>% 
    filter(!(pathotype %in% c("NMEC", "EIEC"))) %>% 
    pivot_wider(names_from = "cluster", values_from = "n", values_fill = 0) 
  
  clustering_paths <- as.dendrogram(hclust(dist(as.matrix(x)[,2:ncol(x)]), method = "ward.D2"))
  path_order <- order.dendrogram(clustering_paths)
  
  dendro <- clustering_paths %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("", position = "right") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_clusters))) + 
    theme_void() + 
    coord_flip() 
  
  
  cl_cl_path2 <- cl_cl_path %>% 
    filter(!(pathotype %in% c("NMEC", "EIEC")))
  cl_cl_path2[["pathotype"]] <- factor(cl_cl_path2[["pathotype"]], levels = x[["pathotype"]][path_order])
  path_plot <- cl_cl_path2 %>% 
    ggplot(aes(x = pathotype, y = n, fill = cluster)) +
    geom_col(position = "fill") +
    theme_bw(base_size = 16) +
    scale_fill_manual("Cluster", values = colors) +
    coord_flip() +
    theme(legend.position = "bottom",
          axis.title.y = element_blank()) +
    ylab("Frequency in cluster")
  
  wrap_plots(list(path_plot, dendro), nrow = 1, widths = c(0.9, 0.1)) &
    theme(legend.position = "bottom")
  
}



get_pathotype_normalized_cluster_plot <- function(clusters_df, clustering, cluster_names, offset) {
  colors <- adhesiomeR::pathotype_colors
  
  cols_to_remove <- colnames(clusters_df)[which(!(colnames(clusters_df) %in% colnames(clustering[["data"]])))]
  x_cl <- clusters_df %>% 
    mutate(cluster = cluster_names[Cluster]) %>% 
    select(-all_of(cols_to_remove)) %>% 
    unique() %>% 
    mutate(genome = 1:nrow(.))
  
  means_df <- lapply(unique(x_cl[["cluster"]]), function(i) {
    x_cl %>% 
      filter(cluster == i) %>% 
      select(-c("cluster", "genome")) %>% 
      colMeans() %>%  
      t() %>% 
      data.frame(check.names = FALSE) %>% 
      mutate(cluster = i)
  }) %>% bind_rows()
  
  clustering_clusters <- as.dendrogram(hclust(dist(as.matrix(select(means_df, -"cluster"))), method = "ward.D2"))
  cluster_order <- order.dendrogram(clustering_clusters)
  genomes_df <- clusters_df %>% 
    mutate(cluster = cluster_names[Cluster]) %>% 
    group_by(cluster) %>% 
    summarise(genomes = n())
  size_df <- data.frame(cluster = cluster_names, clustering[["clusinfo"]]) %>% 
    left_join(genomes_df, by = "cluster")
  
  dendro <- clustering_clusters %>%
    dendro_data %>%
    segment %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_segment() +
    scale_y_continuous("") +
    scale_x_discrete("",
                     limits = factor(1L:nobs(clustering_clusters))) + 
    theme_void() + 
    coord_flip() + 
    scale_y_reverse()
  
  cl_cl_path <- summarise(group_by(mutate(clusters_df, cluster = cluster_names[Cluster]), cluster, pathotype), n = n())
  cl_cl_path[["pathotype"]] <- factor(cl_cl_path[["pathotype"]], levels = c("aEPEC", "tEPEC", "EAEC", "EHEC", "EIEC",
                                                                            "ETEC", "DAEC", "STEC", "UPEC", "APEC",
                                                                            "NMEC", "Nonpathogenic", "NA"))
  path_sizes <- summarise(group_by(cl_cl_path, pathotype), count = sum(n)) %>% 
    mutate(fctr = 1000/count)
  path_plot <- cl_cl_path %>% 
    left_join(path_sizes, by = "pathotype") %>% 
    mutate(size = n*fctr) %>% 
    filter(!(pathotype %in% c("EIEC", "NMEC"))) %>% 
    ggplot(aes(x = cluster, y = size, fill = pathotype)) +
    geom_col(position = "fill") +
    theme_bw(base_size = 16) +
    scale_fill_manual("Pathotype", values = colors) +
    coord_flip() +
    theme(legend.position = "bottom",
          axis.title.y = element_blank()) +
    ylab("")
  
  size_plot <- size_df %>% 
    select(-all_of(c("max_diss", "av_diss", "isolation"))) %>% 
    pivot_longer(c("size", "genomes"), names_to = "type", values_to = "value") %>% 
    mutate(type = ifelse(type == "size", "Profiles", "Genomes")) %>% 
    ggplot(aes(x = cluster, y = value, label = value, fill = type)) +
    geom_col(position = "dodge") +
    theme_bw(base_size = 16) +
    geom_text(hjust = -0.05, size = 4, position = position_dodge(width = .9)) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ylab("Count") +
    ylim(c(0, max(size_df[["genomes"]]+offset))) +
    scale_fill_manual("Number of", values = c("Genomes" = "grey30", "Profiles" = "grey70"))
  
  wrap_plots(list(dendro, path_plot, size_plot), nrow = 1, widths = c(0.03, 0.67, 0.3)) &
    theme(legend.position = "bottom")
}


get_importance_plot <- function(cl_df, rows = 2) {
  importance_all <- lapply(unique(cl_df[["cluster"]]), function(i) {
    dat <- mutate(cl_df, cluster = ifelse(cluster == i, TRUE, FALSE))
    rf <- ranger(data = data.frame(dat), 
                 formula = cluster ~ ., 
                 write.forest = TRUE, 
                 classification = TRUE, 
                 importance = "impurity",
                 seed = 73607254)
    print(paste0(i, ": ", rf[["prediction.error"]]))
    imp_sorted <- sort(rf[["variable.importance"]], decreasing = TRUE)[1:20]
    data.frame(Cluster = i,
               Gene = names(imp_sorted),
               Importance = imp_sorted,
               row.names = NULL)
  }) %>% bind_rows() %>% 
    mutate(Cluster = as.factor(Cluster))
  importance_all[["Gene"]] <- sapply(importance_all[["Gene"]], function(i) gsub(".", "/", i, fixed = TRUE))
  
  ggplot(importance_all, aes(Importance, reorder_within(Gene, Importance, Cluster))) +
    geom_col() +
    facet_wrap(~Cluster, scales = "free_y", nrow = rows) +
    scale_y_reordered() +
    theme_bw() +
    ylab("Gene") +
    xlab("Gini importance")
}