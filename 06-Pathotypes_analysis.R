library(adhesiomeR)

genomes_path <- "/home/ubuntu/hpc-home/RefSeq_Ecoli/genomes/"
results_path <- "/home/ubuntu/adhesiomeR_results_strict/"
pathotypes <- c("aEPEC", "DAEC", "EAEC", "EHEC", "EIEC", "ETEC", "NMEC", "NA", "Nonpathogenic", "tEPEC", "STEC", "UPEC", "APEC")

lapply(pathotypes, function(ith_pathotype) {
  print(paste0("Starting ", ith_pathotype, " pathotype"))
  genome_files <- list.files(path = paste0(genomes_path, ith_pathotype), 
                             full.names = TRUE)
  blast_results <- get_blast_res(genome_files, 8)
  write.csv(blast_results, paste0(results_path, ith_pathotype, "_blast_results.csv"), row.names = FALSE)
  presence_tab <- get_presence_table_strict(blast_results, n_threads = 8)
  write.csv(presence_tab, paste0(results_path, ith_pathotype, "_gene_presence.csv"), row.names = FALSE)
  presence_tab_counts <- get_presence_table_strict(blast_results, count_copies = TRUE, n_threads = 8)
  write.csv(presence_tab_counts, paste0(results_path, ith_pathotype, "_gene_copies.csv"), row.names = FALSE)
  summary_tab <- get_summary_table(presence_tab)
  write.csv(summary_tab, paste0(results_path, ith_pathotype, "_system_presence.csv"), row.names = FALSE)
})

results_path2 <- "/home/ubuntu/adhesiomeR_results_relaxed/"

lapply(pathotypes, function(ith_pathotype) {
  print(paste0("Starting ", ith_pathotype, " pathotype"))
  genome_files <- list.files(path = paste0(genomes_path, ith_pathotype), 
                             full.names = TRUE)
  blast_results <- get_blast_res(genome_files, 8)
  write.csv(blast_results, paste0(results_path2, ith_pathotype, "_blast_results.csv"), row.names = FALSE)
  presence_tab <- get_presence_table_relaxed(blast_results, n_threads = 8)
  write.csv(presence_tab, paste0(results_path2, ith_pathotype, "_gene_presence.csv"), row.names = FALSE)
  presence_tab_counts <- get_presence_table_relaxed(blast_results, count_copies = TRUE, n_threads = 8)
  write.csv(presence_tab_counts, paste0(results_path2, ith_pathotype, "_gene_copies.csv"), row.names = FALSE)
  summary_tab <- get_summary_table(presence_tab)
  write.csv(summary_tab, paste0(results_path2, ith_pathotype, "_system_presence.csv"), row.names = FALSE)
})
