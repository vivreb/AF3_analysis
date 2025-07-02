library(protti)
library(bio3d)
library(cry)
library(dplyr)
library(msa)
library(stringr)
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(RColorBrewer)

source("functions_for_AF3_structural_analysis.R")

files <- list.files(path = "all_top_models/")

ranked_0_files <- data.frame(filename = files) %>% 
  filter(grepl("*\\.pdb", filename)) %>% 
  mutate(uniprot_id = str_extract(filename, "[a-z0-9]*")) %>% 
  mutate(ligand = str_extract(filename, "\\_[a-z]*\\_(ligand\\_)?"))

complete_files <- data.frame(uniprot_id = c(rep(tolower(kinases_to_predict_first %>% pull(uniprot_id) %>% unique()), 3)),
                             ligand = c(rep("_atp_", 410), rep("_stu_", 410), rep("_no_ligand_", 410))) %>% 
  left_join(ranked_0_files, by = c("uniprot_id", "ligand"))

usable_kinases_atp <- complete_files %>% 
  filter(!is.na(filename)) %>%
  filter(ligand %in% c("_atp_", "_no_ligand_")) %>% 
  group_by(uniprot_id) %>% 
  filter(n() == 2) %>% 
  ungroup()

usable_kinases_stu <- complete_files %>% 
  filter(!is.na(filename)) %>%
  filter(ligand %in% c("_stu_", "_no_ligand_")) %>% 
  group_by(uniprot_id) %>% 
  filter(n() == 2) %>% 
  ungroup()

prediction_error_stu <- calculate_reg_domain_prediction_error(data = usable_kinases_stu,
                                                              protein_id_column = uniprot_id,
                                                              filename_column = filename,
                                                              path = "Z:/Viviane/Computational/AF3/all_top_models")

prediction_error_stu %>% write.csv(paste(getwd(), "/stu_to_no_ligand_comparison_prediction_error.csv", sep = ""))


prediction_error_stu %>% 
  ggplot(aes(x = rmssd -rmssd_cat)) +
  geom_histogram(binwidth = 0.25)


change_in_sasa_stu <- calculate_change_in_surface_accessibility(data = usable_kinases_stu, 
                                                                protein_id_column = uniprot_id,
                                                                filename_column = filename,
                                                                all_domain_annotations = kinases_to_predict_first,
                                                                path = "Z:/Viviane/Computational/AF3/all_top_models")

change_in_sasa_stu %>% write.csv(paste(getwd(), "/stu_to_no_ligand_comparison_change_in_sasa.csv", sep = ""))




rmsd_and_surface_stu <- prediction_error_stu %>%
  mutate(kinase_domain_filename = str_replace(kinase_domain_filename, "\\_model.pdb", "")) %>% 
  left_join((change_in_sasa_stu %>% 
               mutate(kinase_domain_file = str_replace(kinase_domain_file, "\\_model.csv", ""))), 
            by = c("uniprot_id", "kinase_domain_filename" = "kinase_domain_file"))



annotated_test_rmsd_and_surface <- rmsd_and_surface_stu %>% 
  left_join((human_kinome %>% 
               mutate(pg_protein_accessions = tolower(pg_protein_accessions)) %>% 
               distinct(pg_protein_accessions, gene_name, kinase_group, kinase_family)), by = c("uniprot_id" = "pg_protein_accessions")) %>% 
  group_by(kinase_group) %>% 
  mutate(label_group = ifelse(n() <= 2, "Other", kinase_group)) %>% 
  ungroup() %>% 
  mutate(changing = ifelse(pae_weighted_rmssd > 2, "yes", "no")) %>% 
  mutate(changing = ifelse(abs(change_in_residue_surface) > 75, "yes", changing)) 


annotated_test_rmsd_and_surface %>% write.csv("Z:/Viviane/Computational/AF3/all_kinases_pae_weighted_rmssd_and_surface_stu.csv")