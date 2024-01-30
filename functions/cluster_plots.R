
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
  colors <- if(colors == "all") {
    adhesiomeR::pathotype_colors
  } else {
    
  }
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
  