library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readxl)

source("functions/process_results.R")
source("functions/benchmark_functions.R")

outpath <- "C:/Users/kar23yal/Dropbox/adhesiomeR/Benchmarks/"
metadata_path <- "C:/Users/kar23yal/Dropbox/adhesiomeR/ETEC_human/metadata/"
etec_results_system_level <- "C:/Users/kar23yal/Dropbox/adhesiomeR/ETEC_human/adhesiomer strict results/ETEC_human_system_presence.csv"

sra_df <- read.csv(paste0(metadata_path, "SraRunInfo.csv"))
annot_df <- read_excel(paste0(metadata_path, "Table_S2_original.xlsx"))
annot_df[["Isolates (UG designation)"]] <- sapply(annot_df[["Isolates (UG designation)"]], function(i) {gsub("sc", "", i)})
nr_df <- read_excel(paste0(metadata_path, "Table_S3.xlsx"))

annot_nr <- left_join(annot_df, nr_df, by = c("Isolates (UG designation)" = "Isolate (UG designation)"))
annot_all <- left_join(annot_nr, select(unique(sra_df), c("Run", "Sample")), by = c("ERS-number" = "Sample")) %>% 
  mutate(File = paste0(Run, ".fa")) %>% 
  filter(!is.na(Run))

adhesiomer_res_etec <- read.csv(etec_results_system_level, 
                                check.names = FALSE)


plot_df <- left_join(select(unique(annot_all), all_of(c("CF profile", "File"))),
                     adhesiomer_res_etec, by = "File")

pcr_info <- lapply(1:nrow(plot_df), function(i) {
  s <- strsplit(plot_df[["CF profile"]][i], "+", fixed = TRUE)[[1]]
  data.frame(System = s,
             Presence = "Present",
             File = plot_df[["File"]][i],
             Results = "Experimental")
}) %>% bind_rows() %>% 
  pivot_wider(names_from = "System", values_from = "Presence", values_fill = "Absent") %>% 
  pivot_longer(3:ncol(.), values_to = "Presence", names_to = "System") %>% 
  filter(System != "CF NEG") 


# Confusion matrix - select systems analysed in the paper
test_df <- select(plot_df, all_of(c("File", "CS1", "CS3", "CS21", "CS2", "CFA/I", "CS7", "CS6", 
                                    "CS8", "CS5", "CS17", "CS19", "CS23", "CS12", "CS14", "CS13", 
                                    "CS28A", "CS20", "CS4", "CS18"))) %>% 
  pivot_longer(2:ncol(.), names_to = "System", values_to = "Presence") %>% 
  mutate(Results = "adhesiomeR",
         System = ifelse(System == "CS28A", "CS28a", System)) %>% 
  bind_rows(pcr_info) %>% #remove one file without annotations
  filter(File != "ERR119532.fa") 

table("Experimental" = filter(test_df, Results == "Experimental")[["Presence"]],
      "adhesiomeR" = filter(test_df, Results == "adhesiomeR")[["Presence"]])


cm_plot <- lapply(unique(test_df[["System"]]), function(ith_system) {
  lapply(unique(test_df[["Results"]]), function(ith_results) {
    x <- filter(test_df, System == ith_system, Results == ith_results) %>%  
      group_by(Presence) %>% 
      summarise(n = n())
    data.frame(Results = ith_results,
               System = ith_system,
               x)
  }) %>% bind_rows()
}) %>% bind_rows() %>% 
  mutate(Presence = factor(Presence, levels = c("Present", "Partial", "Absent"))) %>% 
  ggplot(aes(x = Presence, y = Results, fill = n, label = n)) +
  geom_tile(color = "grey40") +
  facet_wrap(~System, nrow = 4) +
  scale_fill_gradient2("Number\nof strains") +
  geom_text() +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90)) 

ggsave("C:/Users/kar23yal/Dropbox/adhesiomeR/Benchmarks/ETEC/Confusion_matrix.png", cm_plot, width = 10, height = 7)
ggsave("C:/Users/kar23yal/Dropbox/adhesiomeR/Benchmarks/ETEC/Confusion_matrix.eps", cm_plot, width = 10, height = 7)

cm_plot_all <- lapply(unique(test_df[["Results"]]), function(ith_results) {
  x <- filter(test_df, Results == ith_results) %>%  
    group_by(Presence) %>% 
    summarise(n = n())
  data.frame(Results = ith_results,
             x)
}) %>% bind_rows() %>% 
  mutate(Presence = factor(Presence, levels = c("Present", "Partial", "Absent"))) %>% 
  ggplot(aes(x = Presence, y = Results, fill = n, label = n)) +
  geom_tile(color = "grey40") +
  scale_fill_gradient2("Number\nof strains") +
  geom_text() +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90)) 
ggsave(paste0(outpath, "ETEC/Confusion_matrix_all.png"), cm_plot_all, width = 8, height = 4)
ggsave(paste0(outpath, "ETEC/Confusion_matrix_all.eps"), cm_plot_all, width = 8, height = 4)


cm_plot_all_simplified <- lapply(unique(test_df[["Results"]]), function(ith_results) {
  x <- filter(test_df, Results == ith_results) %>% 
    mutate(Presence_simplified = ifelse(Presence == "Present", "Present", "Not present")) %>% 
    group_by(Presence_simplified) %>% 
    summarise(n = n())
  data.frame(Results = ith_results,
             x)
}) %>% bind_rows() %>% 
  mutate(Presence_simplified = factor(Presence_simplified, levels = c("Present", "Not present"))) %>% 
  ggplot(aes(x = Presence_simplified, y = Results, fill = n, label = n)) +
  geom_tile(color = "grey40") +
  scale_fill_gradient2("Number\nof strains") +
  geom_text() +
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 
ggsave(paste0(outpath, "ETEC/Confusion_matrix_all_simplified.png"), cm_plot_all_simplified, width = 4, height = 4)
ggsave(paste0(outpath, "ETEC/Confusion_matrix_all_simplified.eps"), cm_plot_all_simplified, width = 4, height = 4)
ggsave(paste0(outpath, "ETEC/Confusion_matrix_all_simplified.svg"), cm_plot_all_simplified, width = 4, height = 4)


# CF neg
# remove one file without annotations
cf_neg_df <- filter(plot_df, `CF profile` == "CF NEG", File != "ERR119532.fa") %>% 
  select(all_of(c("File", "CFA/I", "CS1", "CS12", "CS13", "CS14", "CS17", "CS18", "CS19", "CS2", 
                  "CS20", "CS21", "CS23", "CS26", "CS28A", "CS28B", "CS3", "CS30", "CS31A", 
                  "CS4", "CS5", "CS6", "CS7", "CS8", "CS27A/CS27B", "CS28A")))
cf_neg_df[["Present"]] <- sapply(cf_neg_df[["File"]], function(i) {
  sum(filter(plot_df, File == i)[2:ncol(cf_neg_df)] == "Present")
})  
cf_neg_df[["Partial"]] <- sapply(cf_neg_df[["File"]], function(i) {
  sum(filter(plot_df, File == i)[2:ncol(cf_neg_df)] == "Partial")
})  

sum(cf_neg_df[["Present"]] > 0, na.rm = TRUE)
sum(cf_neg_df[["Present"]] > 0, na.rm = TRUE)/nrow(cf_neg_df)
sum(cf_neg_df[["Partial"]] > 0, na.rm = TRUE)
sum(cf_neg_df[["Partial"]] > 0, na.rm = TRUE)/nrow(cf_neg_df)


# CF neg - only systems they checked
cf_neg_df2 <- filter(plot_df, `CF profile` == "CF NEG", File != "ERR119532.fa") %>% 
  select(all_of(c("File", "CFA/I", "CS1", "CS12", "CS13", "CS14", "CS17", "CS18", "CS19", "CS2", 
                  "CS20", "CS21", "CS23", "CS28A",  "CS3", 
                  "CS4", "CS5", "CS6", "CS7", "CS8",  "CS28A")))
cf_neg_df2[["Present"]] <- sapply(cf_neg_df2[["File"]], function(i) {
  sum(filter(plot_df, File == i)[2:ncol(cf_neg_df2)] == "Present")
})  
cf_neg_df2[["Partial"]] <- sapply(cf_neg_df2[["File"]], function(i) {
  sum(filter(plot_df, File == i)[2:ncol(cf_neg_df2)] == "Partial")
})  

sum(cf_neg_df2[["Present"]] > 0)
sum(cf_neg_df2[["Present"]] > 0)/nrow(cf_neg_df2)
sum(cf_neg_df2[["Partial"]] > 0)
sum(cf_neg_df2[["Partial"]] > 0)/nrow(cf_neg_df2)
