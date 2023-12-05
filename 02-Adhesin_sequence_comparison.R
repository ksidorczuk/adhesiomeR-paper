library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)


if(Sys.info()[["nodename"]] == "kasia-MACH-WX9") {
  data_path <- "/media/kasia/Data/Dropbox/adhesiomeR/"
}

if(Sys.info()[["nodename"]] == "ryzen") {
  data_path <- "~/Dropbox/adhesiomeR/"
}

if(Sys.info()[["nodename"]] == "N160246") {
  data_path <- "C:/Projects/adhesiomeR/"
}

plot_selected_systems <- function(dat) {
  ggplot(dat, aes(x = query, y = subject, fill = identity_percent)) +
    geom_tile(color = "white") +
    geom_text(aes(label = identity_percent), size = 2) +
    facet_wrap(~System, scales = "free") +
    scale_fill_gradient(low = "white", high = "#e00e00") +
    theme_bw() +
    theme(axis.text = element_text(size = 6))
}    

adhesin_dat <- read_xlsx(paste0(data_path, "adhesiomeR_database/Adhesins.xlsx"), sheet = "Genes")

blast_res <- read.delim(paste0(data_path, "adhesiomeR_database/BLAST/adhesiomeR_sequences_blast_res.tsv"), header = FALSE) %>% 
  setNames(c("query", "subject", "identity_percent", "alignment_length", "mismatches", "gap_opens", 
             "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score")) %>% 
  mutate(query = sapply(.[["query"]], function(i) strsplit(i, "~~~")[[1]][2]),
         subject = sapply(.[["subject"]], function(i) strsplit(i, "~~~")[[1]][2]))

# Large plots with comparisons - full and reduced 
clus_dat <- select(blast_res, c(query, subject, identity_percent)) 
clus_dat_pivoted <- clus_dat %>% 
  pivot_wider(query, names_from = subject, values_from = identity_percent, values_fn = max)

clustering_all <- clus_dat[!duplicated(select(clus_dat, c(query, subject))),] %>% 
  left_join(select(adhesin_dat, c(System, Gene, `Function/Class`)), by = c("query" = "Gene")) %>% 
  ggplot(aes(x = query, y = subject, fill = identity_percent)) +
  geom_tile(color = "white") +
  facet_wrap(~System, scales = "free") +
  scale_fill_gradient(low = "white", high = "#e00e00") +
  theme_bw() +
  theme(axis.text = element_text(size = 6))

clustering_reduced <- clus_dat[!duplicated(select(clus_dat, c(query, subject))),] %>% 
  filter(query != subject) %>% 
  left_join(select(adhesin_dat, c(System, Gene, `Function/Class`)), by = c("query" = "Gene")) %>% 
  ggplot(aes(x = query, y = subject, fill = identity_percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = identity_percent), size = 3) +
  facet_wrap(~System, scales = "free") +
  scale_fill_gradient(low = "white", high = "#e00e00") +
  theme_bw() +
  theme(axis.text = element_text(size = 8))

ggsave(filename = "Adhesin_blast_comparison_all.png", clustering_all, width = 60, height = 45, limitsize = FALSE, path = paste0(data_path, "adhesiomeR_database/BLAST/"))
ggsave(filename = "Adhesin_blast_comparison_reduced.png", clustering_reduced, width = 40, height = 30, limitsize = FALSE, path = paste0(data_path, "adhesiomeR_database/BLAST/"))


# Comparisons of selected systems for easier analysis
clus_dat2 <- select(blast_res, c(query, subject, identity_percent, e_value)) 
dat <- clus_dat2[!duplicated(select(clus_dat2, c(query, subject))),] %>% 
  filter(e_value < 10^-100) %>% 
  left_join(select(adhesin_dat, c(System, Gene, `Function/Class`)), by = c("query" = "Gene"))

afas <- dat %>% 
  filter(System %in% c("Afa-I", "Afa-III", "Afa-VIII", "F1845", "Dr")) %>% 
  plot_selected_systems()
ggsave(filename = "Afa_systems.png", plot = afas, width = 10, height = 6, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs5_cs7 <- dat %>% 
  filter(System %in% c("CS5", "CS7")) %>% 
  plot_selected_systems()
ggsave(filename = "CS5_CS7_systems.png", plot = cs5_cs7, width = 10, height = 6, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

k88 <- dat %>% 
  filter(System %in% c("F41", "K88", "Lda", "CS31A", "CS23", "CS13")) %>% 
  plot_selected_systems()
ggsave(filename = "K88_systems.png", plot = k88, width = 14, height = 9, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs8_cs21 <- dat %>% 
  filter(System %in% c("CS8", "CS21")) %>% 
  plot_selected_systems()
ggsave(filename = "CS8_CS21_systems.png", plot = cs8_cs21, width = 10, height = 6, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs12_cs20 <- dat %>% 
  filter(System %in% c("CS12", "CS20", "CS28A", "CS28B")) %>% 
  plot_selected_systems()
ggsave(filename = "CS12_CS20_systems.png", plot = cs12_cs20, width = 12, height = 8, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs1_cs17 <- dat %>% 
  filter(System %in% c("CS1", "CS17", "CS19", "PCF071")) %>% 
  plot_selected_systems()
ggsave(filename = "CS1_CS17_systems.png", plot = cs1_cs17, width = 9, height = 6, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs18_cs26 <- dat %>% 
  filter(System %in% c("CS18", "CS26", "CS27A", "CS27B", "CS30")) %>% 
  plot_selected_systems()
ggsave(filename = "CS18_CS26_systems.png", plot = cs18_cs26, width = 12, height = 9, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

cs4_cs14 <- dat %>% 
  filter(System %in% c("CS4", "CS14", "CFA/I")) %>% 
  plot_selected_systems()
ggsave(filename = "CS4_CS14_systems.png", plot = cs4_cs14, width = 10, height = 4, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

f17 <- dat %>% 
  filter(System %in% c("F17a", "F17b", "F17d")) %>% 
  plot_selected_systems()
ggsave(filename = "F17_systems.png", plot = f17, width = 10, height = 4, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

p_s <- dat %>% 
  filter(System %in% c("F1C", "S_1", "S_2", "P_1", "P_2")) %>% 
  plot_selected_systems()
ggsave(filename = "P_S_systems.png", plot = p_s, width = 14, height = 9, path = paste0(data_path, "adhesiomeR_database/BLAST/"))

