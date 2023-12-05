library(dplyr)
library(tidyr)
library(purrr)
library(ggdendro)
library(umap)

if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/adhesiomeR/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/adhesiomeR/"
}

source("functions/assembly_info_qc.R")
source("functions/get_ftp_links.R")
source("functions/pathotyping.R")

# Read in necessary files
RefSeq_assembly_summary_file <- paste0(data_path, "RefSeq_Ecoli_database/RefSeq_Ecoli_assembly_summary.txt")
RefSeq_assembly_summary <- read.delim(RefSeq_assembly_summary_file, skip = 1)

RefSeq_assembly_info_file <- paste0(data_path, "RefSeq_Ecoli_database/RefSeq_Ecoli_assembly_info.csv")
RefSeq_assembly_info <- read.csv(RefSeq_assembly_info_file)

platforms_file <- paste0(data_path, "RefSeq_Ecoli_database/Platforms_to_filter_out.txt")


# Filter RefSeq genomes
RefSeq_assembly_info_after_qc <- do_assembly_info_qc(RefSeq_assembly_info = RefSeq_assembly_info, 
                                                     platforms_file = platforms_file)
write.csv(RefSeq_assembly_info_after_qc, 
          paste0(data_path, "RefSeq_Ecoli_database/RefSeq_Ecoli_assembly_info_filtered.csv"),
          row.names = FALSE)

# Get ftp links for filtered RefSeq
get_ftp_links(assembly_info_filtered = RefSeq_assembly_info_after_qc,
              assembly_summary = RefSeq_assembly_summary, 
              outfile = paste0(data_path, "RefSeq_Ecoli_database/RefSeq_Ecoli_filtered_ftp_links.txt"))

