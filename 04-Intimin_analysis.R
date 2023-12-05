library(dplyr)
library(parallel)
library(future)
library(future.apply)
library(pbapply)
library(biogram)
library(tidyr)
library(ggplot2)
library(pals)
library(plotly)

source("functions/colocalization_functions.R")
source("functions/process_results.R")

do_blast_single <- function(input_file, blast_dir = Sys.which("blastn")) {
  input_file_name <- last(strsplit(input_file, "/")[[1]])
  system(paste0(blast_dir, " -db ~/Programy/ncbi-blast-2.11.0+/db/intimin -query ", gsub(" ", "\\ ", input_file, fixed = TRUE), " -out ", gsub(" ", "\\ ", input_file_name, fixed = TRUE), ".blast -outfmt 6"))
  res <- tryCatch(
    read.delim(paste0(input_file_name, ".blast"), header = FALSE), 
    error = function(e) {
      msg <- conditionMessage(e)
      if(msg == "no lines available in input") as.data.frame(matrix(nrow = 1, ncol = 12))
    })
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(paste0(input_file_name, ".blast"))
  mutate(res, File = input_file_name)
}

get_blast_res <- function(input_file_list, nt = 1, blast_dir = Sys.which("blastn")) {
  res <- pblapply(1:length(input_file_list), cl = nt, function(i) {
    do_blast_single(input_file_list[[i]], blast_dir)
  }) %>% bind_rows() 
}

genomes_path <- "~/RefSeq_Ecoli/genomes/"
results_path <- "~/Dropbox/adhesiomeR/VAF_variants/intimin/Blast_results/"
pathotypes <- c("APEC", "aEPEC", "DAEC", "EAEC", "EHEC", "EIEC", "ETEC", "NMEC", "NA", "Nonpathogenic", "tEPEC", "STEC", "UPEC")

lapply(pathotypes, function(ith_pathotype) {
  print(paste0("Starting ", ith_pathotype, " pathotype"))
  genome_files <- list.files(path = paste0(genomes_path, ith_pathotype), 
                             full.names = TRUE)
  blast_results <- get_blast_res(genome_files, 20, blast_dir = "~/Programy/ncbi-blast-2.11.0+/bin/blastn")
  write.csv(blast_results, paste0(results_path, ith_pathotype, "_blast_results.csv"), row.names = FALSE)
})

# Combine results into a single summary file
gene_seqs <- read_fasta("~/Dropbox/adhesiomeR/VAF_variants/intimin/eae.txt")
gene_lengths <- data.frame(Subject = sapply(names(gene_seqs), function(i) strsplit(i, " ")[[1]][1], USE.NAMES = FALSE),
                           Length = lengths(gene_seqs))

res_files <- list.files(results_path, full.names = TRUE)
res_summary <- lapply(res_files, function(ith_file) {
  x <- read.csv(ith_file, check.names = FALSE) %>% 
    left_join(gene_lengths) %>% 
    mutate(coverage = `Alignment length`/Length*100,
           Pathotype = gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]]))) %>% 
    filter(`% identity` > 50 & coverage > 50 & !is.na(Query)) %>% 
    mutate(Query = as.character(Query))
}) %>% bind_rows()

write.csv(res_summary, "~/Dropbox/adhesiomeR/VAF_variants/intimin/eae_blast_results_all.csv", row.names = FALSE)
