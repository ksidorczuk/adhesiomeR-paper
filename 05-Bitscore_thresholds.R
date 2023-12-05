library(dplyr)
library(ggplot2)
library(adhesiomeR)

if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Users/kar23yal/Dropbox/adhesiomeR/"
}

source("functions/process_results.R")

autotransporters <- c("tia", "iha", "fdeC", "etpA", "saa", "efa1/lifA", "toxB", 
                      "paa", "tibA", "aidA", "ycgV", "flu", "yfaL", "yeeJ", 
                      "ypjA", "upaG", "cah", "espP", "sab", "ehaG", "ehaA", 
                      "ehaB", "aatA", "aatB", "tsh", "eibG")

# Fimbrial adhesins
fimbrial_ref_files <- list.files(paste0(data_path, "Gene_localization_operons/Reference_sequences"),
                                 pattern = "_curated.csv", full.names = TRUE)

fimbrial_ref <- lapply(fimbrial_ref_files, function(ith_file) {
  read.csv(ith_file, check.names = FALSE)
}) %>% bind_rows() %>% 
  filter(Category == "Good")

fimbrial_thresholds <- lapply(unique(fimbrial_ref[["Gene"]]), function(ith_gene) {
  thr <- min(filter(fimbrial_ref, Gene == ith_gene)[["Bit score"]])
  data.frame(Gene = ith_gene,
             Threshold = thr)
}) %>% bind_rows()

# Autotransporters
at_blast <- read_in_blast_results(paste0(data_path, "Pathotypes/adhesiomeR_results/"))  %>% 
  left_join(adhesins_lengths, by = c("Subject" = "Gene")) %>% 
  mutate(`Subject coverage` = `Alignment length`/Length*100) %>% 
  filter(`Subject coverage` > 80 & Subject %in% autotransporters)
at_info <- read.csv(paste0(data_path, "Overlapping_genes/autotransporters_thresholds.csv"))

at_ref <- left_join(at_blast, at_info, by = c("Subject" = "Gene")) %>% 
  filter(`% identity` >= identity)

at_thresholds <- lapply(unique(at_ref[["Subject"]]), function(ith_gene) {
  thr <- min(filter(at_ref, Subject == ith_gene)[["Bit score"]])
  data.frame(Gene = ith_gene,
             Threshold = thr)
}) %>% bind_rows()

# Intimin
int_blast <- read.csv(paste0(data_path, "VAF_variants/intimin/eae_blast_results_all.csv"),
                      check.names = FALSE) 

int_ref <- filter(int_blast, 
                  Subject == "M58154.1:1644-4463_eae-alpha1",
                  `% identity` > 70, coverage > 70)

int_thresholds <- data.frame(Gene = "eae",
                             Threshold = min(int_ref[["Bit score"]]))


# Saa
saa_blast <- read.csv(paste0(data_path, "Overlapping_genes/25.08.2023_saa_variants_blast_filter_saa_5.csv"),
                      check.names = FALSE) %>% 
  setNames(c("Query", "Subject", "% identity", "Alignment length", "Mismatches", 
             "Gap opens", "Query start", "Query end", "Subject start", "Subject end", 
             "Evalue", "Bit score", "Subject coverage"))

saa_ref <- filter(saa_blast,
                  `% identity` > 95,
                  `Subject coverage` > 78)

saa_thresholds <- data.frame(Gene = "saa",
                             Threshold = min(saa_ref[["Bit score"]]))


# Combine all thresholds and save 

all_thresholds <- bind_rows(list(fimbrial_thresholds, at_thresholds, 
                                 int_thresholds, saa_thresholds))
write.csv(all_thresholds, paste0(data_path, "adhesiomeR_database/Bitscore_thresholds.csv"),
          row.names = FALSE)



# Table with all reference hits for supplements
all_blast_res <- read_in_blast_results(paste0(data_path, "Pathotypes/adhesiomeR_results_strict/"))
all_refs <- bind_rows(list(select(fimbrial_ref, -c("Category", "Comment", "start", "end", "Position", "Position2", "Name", "pathotype", "System")), 
                           select(mutate(at_ref, Gene = Subject), -c("Subject", "Mismatches", "Gap opens", "Subject start", "Subject end",
                                                                     "pathotype", "Length", "identity", "Subject coverage")), 
                           select(mutate(int_ref, Gene = "eae"), -c("Subject", "Mismatches", "Gap opens", "Subject start", "Subject end",
                                                                    "Pathotype", "Length",  "coverage")),
                           select(mutate(saa_ref, Gene = "saa"), -c("Subject", "Mismatches", "Gap opens", "Subject start", "Subject end", "Subject coverage"))))

write.csv(all_refs, paste0(data_path, "Publication/SupplementaryFiles/Table S5.csv"), row.names = FALSE)
