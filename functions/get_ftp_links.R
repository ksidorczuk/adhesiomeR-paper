extend_ftp_links <- function(partial_links, type = "_genomic.fna.gz") {
  lapply(partial_links, function(ith_link) {
    paste0(ith_link, "/", last(strsplit(ith_link, "/")[[1]]), type)
  }) %>% unlist() 
}

get_ftp_links <- function(assembly_info_filtered, assembly_summary, outfile) {
  filtered_summary <- filter(assembly_summary, `X..assembly_accession` %in% assembly_info_filtered[["assembly"]])
  ftp_links <- extend_ftp_links(filtered_summary[["ftp_path"]])
  writeLines(ftp_links, outfile)
  outfile
}

