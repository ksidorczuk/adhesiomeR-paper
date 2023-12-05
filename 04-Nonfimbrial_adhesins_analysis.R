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
  system(paste0(blast_dir, " -db ~/Programy/ncbi-blast-2.11.0+/db/overlaps -query ", gsub(" ", "\\ ", input_file, fixed = TRUE), " -out ", gsub(" ", "\\ ", input_file_name, fixed = TRUE), ".blast -outfmt 6"))
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
  mutate(res, Subject = sapply(Subject, function(i) strsplit(i, "~~~")[[1]][2]))
}

genomes_path <- "~/RefSeq_Ecoli/genomes/"
results_path <- "~/Dropbox/adhesiomeR/Overlapping_genes/Blast_results/"
plot_path <- "~/Dropbox/adhesiomeR/Overlapping_genes/Colocalization/"
pathotypes <- c("APEC", "aEPEC", "DAEC", "EAEC", "EHEC", "EIEC", "ETEC", "NMEC", "NA", "Nonpathogenic", "tEPEC", "STEC", "UPEC")

lapply(pathotypes, function(ith_pathotype) {
  print(paste0("Starting ", ith_pathotype, " pathotype"))
  genome_files <- list.files(path = paste0(genomes_path, ith_pathotype), 
                             full.names = TRUE)
  blast_results <- get_blast_res(genome_files, 20, blast_dir = "~/Programy/ncbi-blast-2.11.0+/bin/blastn")
  write.csv(blast_results, paste0(results_path, ith_pathotype, "_blast_results.csv"), row.names = FALSE)
})

# Combine results into a single summary file
#gene_seqs <- read_fasta("~/Dropbox/adhesiomeR/Overlapping_genes/17.03.2023_autotransporters+overlapping_genes.txt")
#gene_seqs <- read_fasta("~/Dropbox/adhesiomeR/Overlapping_genes/19.03.2023_autotransporters+overlapping_genes_ver2.txt")
gene_seqs <- read_fasta("~/Dropbox/adhesiomeR/Overlapping_genes/22.04.2023_autotransporters+overlapping_genes_ver3")
gene_lengths <- data.frame(Subject = sapply(names(gene_seqs), function(i) strsplit(i, "~~~")[[1]][2]),
                           Length = lengths(gene_seqs))
gene_names <- c("tia", "iha", "fdeC", "etpA", "saa", "efa1/lifA", "toxB", "paa", "tibA", "aidA", "ycgV", "flu", 
                "yfaL", "yeeJ", "ypjA", "upaG", "cah", "espP", "sab", "ehaG", "ehaA", "ehaB", "aatA", "aatB", "tsh")

res_files <- list.files(results_path, full.names = TRUE)
res_summary <- lapply(res_files, function(ith_file) {
  df <- read.csv(ith_file, check.names = FALSE) %>% 
    left_join(gene_lengths) %>% 
    mutate(coverage = `Alignment length`/Length*100) %>% 
    filter(`% identity` > 80 & coverage > 80) %>%
    mutate(gene = sapply(.[["Subject"]], function(i) strsplit(i, "-")[[1]][1]),
           type = sapply(.[["Subject"]], function(i) strsplit(i, "-")[[1]][2]),
           version = sapply(.[["Subject"]], function(i) strsplit(i, "-")[[1]][3])) %>% 
    mutate(number = sapply(.[["version"]], function(i) strsplit(i, ".", fixed = TRUE)[[1]][1]),
           type = ifelse(is.na(type), "Gene", paste0(type, "-", number))) 
  summ <- df %>% 
    group_by(File, Query, gene, type) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "type", values_from = n) %>% 
    filter(!is.na(Gene)) %>% 
    mutate(Pathotype = gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])))
  
  summ_identity <- df %>% 
    group_by(File, Query, gene, type) %>% 
    summarise(identity = paste0(`% identity`, collapse = ";")) %>% 
    pivot_wider(names_from = "type", values_from = identity, names_prefix = "identity_")
  summ_coverage <- df %>% 
    group_by(File, Query, gene, type) %>% 
    summarise(identity = paste0(`coverage`, collapse = ";")) %>% 
    pivot_wider(names_from = "type", values_from = identity, names_prefix = "coverage_")
  summ_versions <-  df %>% 
    group_by(File, Query, gene, type) %>% 
    summarise(identity = paste0(`version`, collapse = ";")) %>% 
    pivot_wider(names_from = "type", values_from = identity, names_prefix = "version_")
  
  lapply(unique(summ[["gene"]]), function(ith_gene) {
    limits <- df %>%
      filter(gene == ith_gene) %>% 
      group_by(File, Query) %>%
      summarise(start = min(c(`Query end`, `Query start`)),
                end = max(c(`Query end`, `Query start`)),
                good = ifelse(any(Subject == ith_gene), TRUE, FALSE)) %>% 
      filter(good == TRUE)
    
    colocalization_dat <- df %>%
      filter(gene == ith_gene, Query %in% limits[["Query"]]) %>% 
      left_join(limits) %>% 
      mutate(Gene = type,
             Position = `Query start` - start + 1,
             Position2 = `Query end` - start + 1,
             Name = ifelse(type == "contigs", Query, paste0(File, "-", Query))) 
    bad <- unique(filter(colocalization_dat, Position > 10000 | Position2 > 10000)[["Query"]])
    cd_good <- filter(colocalization_dat, !(Query %in% bad))
    cd_bad <- filter(colocalization_dat, Query %in% bad)
    if(nrow(cd_good) > 0) {
      p1 <- plot_operon_colocalization(cd_good, gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), ith_gene)
      if(nrow(cd_good) > 2500) {
        htmlwidgets::saveWidget(ggplotly(p1), paste0(plot_path, gsub("/", "-", ith_gene), "_", gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), "_good.html"), selfcontained = TRUE)
      } else {
        ggsave(paste0(plot_path, gsub("/", "-", ith_gene), "_", gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), "_good.png"), 
               p1, width = log2(max(c(colocalization_dat[["Position"]], colocalization_dat[["Position2"]]))), height = 2+nrow(cd_good)/75, limitsize = FALSE)
      }
    }
    if(nrow(cd_bad) > 0) {
      p2 <- plot_operon_colocalization(cd_bad, gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), ith_gene)
      htmlwidgets::saveWidget(ggplotly(p2), paste0(plot_path, gsub("/", "-", ith_gene), "_", gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), "_bad.html"), selfcontained = TRUE)
      # ggsave(paste0(data_path, "Overlapping_genes/Colocalization/", gsub("/", "-", ith_gene), "_", gsub("_blast_results.csv", "", last(strsplit(ith_file, "/")[[1]])), "_bad.png"), 
      #        p2, width = log2(max(c(colocalization_dat[["Position"]], colocalization_dat[["Position2"]])))*2, height = 2+nrow(cd_bad)/70, limitsize = FALSE)
    }
  })
  left_join(left_join(left_join(summ, summ_versions), summ_identity), summ_coverage)
}) %>% bind_rows()

write.csv(res_summary[, c("Pathotype", "File", "Query", "gene", "Up-2", "Up-1", "Gene", "Down-1", "Down-2", "Down-3", 
                          "version_Up-2", "version_Up-1", "version_Down-1", "version_Down-2", "version_Down-3",
                          "identity_Up-2", "identity_Up-1", "identity_Gene", "identity_Down-1", "identity_Down-2", "identity_Down-3",
                          "coverage_Up-2", "coverage_Up-1", "coverage_Gene", "coverage_Down-1", "coverage_Down-2", "coverage_Down-3")], 
          "~/Dropbox/adhesiomeR/Overlapping_genes/Blast_results_summary_80.3.csv", row.names = FALSE)

adhesiomer_res <- read_in_gene_results("~/Dropbox/adhesiomeR/Pathotypes/adhesiomeR_results/") %>% 
  pivot_longer(2:466, names_to = "gene", values_to = "adhesiomer") %>% 
  setNames(c("File", "Pathotype", "gene", "adhesiomer"))

write.csv(left_join(res_summary[, c("Pathotype", "File", "Query", "gene", "Up-2", "Up-1", "Gene", "Down-1", "Down-2", "Down-3", 
                                    "version_Up-2", "version_Up-1", "version_Down-1", "version_Down-2", "version_Down-3",
                                    "identity_Up-2", "identity_Up-1", "identity_Gene", "identity_Down-1", "identity_Down-2", "identity_Down-3",
                                    "coverage_Up-2", "coverage_Up-1", "coverage_Gene", "coverage_Down-1", "coverage_Down-2", "coverage_Down-3")], 
                    adhesiomer_res, by = c("Pathotype", "File", "gene")),
          "~/Dropbox/adhesiomeR/Overlapping_genes/Blast_results_summary_80.3+adhesiomer.csv", row.names = FALSE)
