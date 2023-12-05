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

read_in_gene_results_identities <- function(adhesiomer_res_path, copies = FALSE) {
  lf <- if(copies == TRUE) {
    list.files(adhesiomer_res_path, pattern = "gene_copies", full.names = TRUE)
  } else {
    list.files(adhesiomer_res_path, pattern = "gene_presence", full.names = TRUE)
  }
  lapply(lf, function(ith_file) {
    read.csv(ith_file, check.names = FALSE) %>% 
      mutate(pathotype = strsplit(last(strsplit(ith_file, "/")[[1]]), "_")[[1]][1],
             threshold = strsplit(last(strsplit(ith_file, "/")[[1]]), "_")[[1]][2])
  }) %>% bind_rows()
}

read_in_system_results <- function(adhesiomer_res_path) {
  lapply(list.files(adhesiomer_res_path, pattern = "system", full.names = TRUE), function(ith_file) {
    read.csv(ith_file, check.names = FALSE) %>% 
      mutate(pathotype = strsplit(last(strsplit(ith_file, "/")[[1]]), "_")[[1]][1])
  }) %>% bind_rows()
}

read_in_blast_results <- function(adhesiomer_res_path) {
  lapply(list.files(adhesiomer_res_path, pattern = "blast", full.names = TRUE), function(ith_file) {
    read.csv(ith_file, check.names = FALSE) %>% 
      mutate(pathotype = strsplit(last(strsplit(ith_file, "/")[[1]]), "_")[[1]][1])
  }) %>% bind_rows()
}
