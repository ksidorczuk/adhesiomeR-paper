# This script was originally used with older version of adhesiomeR (commit 0f675a9).
# Should correspond to current relaxed setting with thresholds of coverage and identity set to 80%
# To run it with current version of adhesiomeR, please use get_presence_table_relaxed function. 
library(adhesiomeR)

genomes_path <- "~/RefSeq_Ecoli/genomes/"
results_path <- "~/Dropbox/adhesiomeR/Pathotypes/adhesiomeR_results/"
pathotypes <- c("aEPEC", "DAEC", "EAEC", "EHEC", "EIEC", "ETEC", "NMEC", "NA", "Nonpathogenic", "tEPEC", "STEC", "UPEC", "APEC")

lapply(pathotypes, function(ith_pathotype) {
  print(paste0("Starting ", ith_pathotype, " pathotype"))
  genome_files <- list.files(path = paste0(genomes_path, ith_pathotype), 
                             full.names = TRUE)
  blast_results <- get_blast_res(genome_files, 12, blast_dir = "~/Programy/ncbi-blast-2.11.0+/bin/")
  presence_tab <- get_presence_table(blast_results, n_threads = 12, identity = 80, coverage = 80)
  presence_tab_counts <- get_presence_table(blast_results, count_copies = TRUE, n_threads = 12, identity = 80, coverage = 80)
  summary_tab <- get_summary_table(presence_tab)
  write.csv(blast_results, paste0(results_path, ith_pathotype, "_blast_results.csv"), row.names = FALSE)
  write.csv(presence_tab, paste0(results_path, ith_pathotype, "_gene_presence.csv"), row.names = FALSE)
  write.csv(presence_tab_counts, paste0(results_path, ith_pathotype, "_gene_copies.csv"), row.names = FALSE)
  write.csv(summary_tab, paste0(results_path, ith_pathotype, "_system_presence.csv"), row.names = FALSE)
})
