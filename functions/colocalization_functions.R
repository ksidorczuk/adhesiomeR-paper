get_operon_fullnes_info <- function(filtered_blast_system_res, genes, type = "contigs") {
  x <- if(type == "contigs") {
    group_by(filtered_blast_system_res, pathotype, File, Query)
  } else {
    group_by(filtered_blast_system_res, pathotype, File)
  }
  summarise(x, 
            n_found = n(),
            n_unique = n_distinct(Gene),
            full = ifelse(n_unique == length(genes), TRUE, FALSE)) 
}

plot_operon_colocalization <- function(plot_data, pathotype, system) {
  ggplot(plot_data) +
    geom_segment(aes(x = Position, xend = Position2, y = Name, yend = Name, color = Gene)) +
    theme_bw(base_size = 3) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 2)) +
    ggtitle(paste(pathotype, system)) +
    scale_color_manual(values = unname(alphabet()))
}



#' @param full_operons_dat a data frame with information on number of found and unique genes
#' in each contig, and information if it is a full operon
#' @param filtered_blast_system_res system presence results filtered to only include current
#' system, joined with blast results (also filtered)
#' @param type 'contigs' for data on operons encoded on a single contig, or 'genomes' for data
#' on operons that may be encoded on multiple contigs
get_colocalization_data <- function(full_operons_dat, filtered_blast_system_res, type = "contigs") {
  df <- if(type == "contigs") {
    filter(filtered_blast_system_res, Query %in% filter(full_operons_dat, full == TRUE)[["Query"]])
  } else {
    filter(filtered_blast_system_res, File %in% filter(full_operons_dat, full == TRUE)[["File"]])
  }
  colocalization_dat <- df %>%
    group_by(File, pathotype, System, Query, Gene) %>% 
    summarise(`% identity` = `% identity`[which.max(`Bit score`)],
              `Bit score` = `Bit score`[which.max(`Bit score`)],
              `Alignment length` = `Alignment length`[which.max(`Bit score`)],
              Evalue = Evalue[which.max(`Bit score`)],
              `Query start` = `Query start`[which.max(`Bit score`)],
              `Query end` = `Query end`[which.max(`Bit score`)])
  
  if(any(is.na(colocalization_dat[["Query start"]]))) {
    print(paste0("here!"))
  }
  if(any(is.na(colocalization_dat[["Query end"]]))) {
    print(paste0("here2!"))
  }
  limits <- colocalization_dat %>%
    group_by(pathotype, File, Query) %>%
    summarise(start = min(c(`Query end`, `Query start`)),
              end = max(c(`Query end`, `Query start`)))
  
  colocalization_dat %>% 
    left_join(limits) %>% 
    mutate(Position = `Query start` - start + 1,
           Position2 = `Query end` - start + 1,
           Name = ifelse(type == "contigs", Query, paste0(File, "-", Query))) 
}
