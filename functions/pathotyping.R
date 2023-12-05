create_blast_res_matrix <- function(blast_res_dir) {
  files <- list.files(blast_res_dir)
  lapply(files, function(i) {
    gene <- strsplit(i, '.', fixed = TRUE)[[1]][1]
    if(length(readLines(paste0(blast_res_dir, i))) == 0) {
      df <- data.frame("accession" = as.character(NA), gene = NA) %>%
        setNames(., c("accession", gene)) %>%
        drop_na()
    } else {
      temp_df <- read.delim(paste0(blast_res_dir, i), header = FALSE, stringsAsFactors = FALSE) %>%
        select(2,11) %>%
        filter(V11 < 10^-100) %>%
        setNames(c("accession", gene))
      df <- filter(temp_df, !duplicated(temp_df[["accession"]]))
    }
    df
  }) %>% reduce(full_join, by = "accession") %>%
    rowwise() %>%
    mutate(pathotype = case_when(
      !is.na(eae) & !is.na(bfpA) & (is.na(stx1a) & is.na(stx1b)) & (is.na(stx2a) & is.na(stx2b)) ~ "tEPEC",
      !is.na(eae) & is.na(bfpA) & (is.na(stx1a) & is.na(stx1b)) & (is.na(stx2a) & is.na(stx2b)) ~ "aEPEC",
      (!is.na(stx1a) & !is.na(stx1b) & !is.na(stx2a) & !is.na(stx2b) & !is.na(eae)) |
        (!is.na(stx1a) & !is.na(stx1b) & !is.na(eae)) | (!is.na(stx2a) & !is.na(stx2b) & !is.na(eae)) ~ "EHEC",
      (!is.na(ial) | !is.na(ipaH)) & is.na(stx1a) & is.na(stx2a)  ~ "EIEC",
      (!is.na(eltA) | !is.na(eltB) | !is.na(sta1) | !is.na(sta2) | !is.na(stb)) & is.na(aggR) ~ "ETEC",
      !is.na(aggR) & !is.na(aatA) ~ "EAEC",
      !is.na(afaA) & (!is.na(afaE_I) | !is.na(afaE_III) | !is.na(daaE) | !is.na(draE)) & !is.na(sat) ~ "DAEC",
      !is.na(fyuA) & !is.na(fimH) &
        ((!is.na(chuA) & !is.na(yfcV)) | (!is.na(chuA) & !is.na(vat)) | (!is.na(yfcV) & !is.na(vat))) ~ "UPEC",
      all(is.na(c(aatA, afaA, eae, ial, ipaH, eltA, eltB, sta1, sta2, stb, aggR, afaE_I, afaE_III, daaE, draE, sat, vat))) == TRUE &
        ((!is.na(stx1a) & !is.na(stx1b)) | (!is.na(stx2a) & !is.na(stx2b))) ~ "STEC",
      all(is.na(c(aatA, afaA, afaE_I, afaE_III, aggR, bfpA, daaE, draE, eae, eltA, eltB, ial, ipaH, 
                  pet, sat, sta1, sta2, stb, stx1a, stx1b, stx2a, stx2b, vat))) == TRUE &
        (all(is.na(c(iucC, neuC, sitA, yfcV))) == TRUE & !is.na(fyuA) |
            all(is.na(c(fyuA, neuC, sitA, yfcV))) == TRUE & !is.na(iucC) |
            all(is.na(c(fyuA, iucC, sitA, yfcV))) == TRUE & !is.na(neuC) |
            all(is.na(c(fyuA, iucC, neuC, yfcV))) == TRUE & !is.na(sitA) |
            all(is.na(c(fyuA, iucC, neuC, sitA))) == TRUE & !is.na(yfcV)) |
        all(is.na(c(aatA, afaA, afaE_I, afaE_III, aggR, bfpA, daaE, draE, eae, eltA, eltB, ial, ipaH, pet, sat, sta1, 
                    sta2, stb, stx1a, stx1b, stx2a, stx2b, vat, fyuA, iucC, neuC, sitA, yfcV)) == TRUE) ~ "nonpathogenic",
      all(is.na(c(sitA, vat, neuC, iucC, neuA)) == FALSE) ~ "NMEC"
    ))
}


plot_genes_for_pathotypes <- function(blast_res_matrix, outdir) {
  blast_res_matrix <- mutate(blast_res_matrix, pathotype = ifelse(is.na(pathotype), "unknown", pathotype))
  lapply(unique(blast_res_matrix[["pathotype"]]), function(ith_pathotype) {
    plot_dat <- filter(blast_res_matrix, pathotype == ith_pathotype) %>% 
      mutate(across(2:32, ~ !is.na(.)))
    dendro_accessions <- as.dendrogram(hclust(d = dist(x = as.matrix(plot_dat[, 2:(ncol(plot_dat)-1)]))))
    accessions_order <- order.dendrogram(dendro_accessions)
    dendro_genes <- as.dendrogram(hclust(d = dist(t(x = as.matrix(plot_dat[, 2:(ncol(plot_dat)-1)])))))
    genes_order <- order.dendrogram(dendro_genes)

    long_dat <- pivot_longer(plot_dat, 2:32, names_to = "Gene", values_to = "E-value") 
    
    long_dat[["accession"]] <- factor(long_dat[["accession"]], 
                                   levels = plot_dat[["accession"]][accessions_order],
                                   ordered = TRUE)
    long_dat[["Gene"]] <- factor(long_dat[["Gene"]], 
                                   levels = colnames(plot_dat)[2:(ncol(plot_dat)-1)][genes_order],
                                   ordered = TRUE)
    heatmap <- ggplot(long_dat, aes(x = Gene, y = accession, fill = `E-value`)) +
      geom_tile() +
      scale_fill_manual(values = c(`TRUE` = "darkgreen", `FALSE` = "grey90")) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45))
    x <- nrow(plot_dat)
    h <- case_when(x >= 1000 ~ x/60, 
                   x >= 100 & x < 1000 ~ x*0.0175+1,
                   x < 100 ~ x*0.1+1)
    
    ggsave(filename = paste0(outdir, "/", ith_pathotype, "_blast_found_genes.png"),
           plot = heatmap,
           height = h,
           width = 10,
           limitsize = FALSE)  
  })
  outdir
}


add_filenames <- function(blast_res_matrix, assembly_summary, outfile) {
  filename_dat <- select(assembly_summary, c(`X..assembly_accession`, asm_name)) %>% 
    mutate(acc = sapply(assembly_summary[["X..assembly_accession"]], function(i) strsplit(i, ".", fixed = TRUE)[[1]][1]),
           filename = gsub(" ", "_", paste0(`X..assembly_accession`, "_", asm_name, "_genomic.fna")))
  with_names <- left_join(blast_res_matrix, select(filename_dat, c(acc, filename)), by = c("accession" = "acc"))
  write.csv(with_names, outfile, row.names = FALSE)
  with_names
}