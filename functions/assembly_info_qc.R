do_assembly_info_qc <- function(RefSeq_assembly_info, platforms_file) {
  x <- readLines(platforms_file)
  RefSeq_assembly_info %>% 
    mutate(coverage = as.numeric(coverage)) %>% 
    filter(!((coverage < 50 & contign50 < 1000000) | platform %in% x | (is.na(platform) & contign50 < 1000000))) %>% 
    filter(contigs <= 400 | is.na(contigs))
}

