library(readxl)
library(tidyr)
library(dplyr)
library(protti)
library(stringr)
library(featurefinder)
library(rjson)

# Load all human kinases (only catalytic domains)

human_kinome <- read_excel("231129_human_kinome.xlsx")

all_kinases <- human_kinome %>%
  pull(pg_protein_accessions)

# Fetch uniprot annotations of domains

uniprot_domains <-
  fetch_uniprot(
    c(all_kinases),
    columns = c(
      "gene_primary",
      "sequence",
      "length",
      "ft_transmem",
      "ft_binding",
      "ft_domain",
      "ft_region",
      "ft_motif",
      "ft_zn_fing"
    ),
  ) %>%
  dplyr::rename(
    length_protein = length,
    uniprot_id = accession
  ) 

#Make list of transmembrane kinases to filter later

transmembrane_proteins <- uniprot_domains %>% 
  filter(!is.na(ft_transmem)) %>% 
  distinct(uniprot_id, gene_primary)

# Extract domain annotations 

domain_annotations <- uniprot_domains %>% 
  mutate(regions = paste(ft_domain, ft_region, ft_motif, ft_zn_fing, sep = "; ")) %>% 
  mutate(regions = str_replace_all(regions, '["]', "")) %>% 
  mutate(regions = str_replace_all(regions, '; \\/', "_")) %>% 
  mutate(regions = str_extract_all(regions, "[A-Z\\_]*\\s[0-9]*(\\.\\.[0-9]*)?\\_[a-z]*=[A-Za-z\\s-\\/0-9\\(\\)]*\\_?")) %>% 
  unnest(regions) %>% 
  distinct(uniprot_id, gene_primary, sequence, length_protein, regions) %>% 
  mutate(region_type = str_extract(regions, "[A-Z\\_]*")) %>% 
  mutate(region_start = str_extract(regions, " [0-9]*")) %>% 
  mutate(region_start = str_replace(region_start, "\\.", "")) %>% 
  mutate(region_start = as.numeric(region_start)) %>% 
  mutate(region_end = str_extract(regions, "\\.\\.[0-9]*")) %>% 
  mutate(region_end = str_replace(region_end, "\\.\\.", "")) %>% 
  mutate(region_end = as.numeric(region_end)) %>% 
  mutate(region_end = ifelse(is.na(region_end), region_start, region_end)) %>% 
  mutate(region_name = str_extract(regions, "=[A-Za-z\\s-\\/0-9\\(\\)]*")) %>% 
  mutate(region_name = str_replace(region_name, "=", "")) %>% 
  dplyr::select(-c(regions))

# Load and prepare data from pubmed on pseudosubstrate domains in kinases

pseudosubstrate_annotations_pubmed <- read_excel("Z:/Viviane/Lists/kinase_pseudosubstrate_domains_literature_search.xlsx") %>% 
  dplyr::select(-c(source)) %>% 
  left_join((domain_annotations %>% 
               dplyr::select(-c(region_type, region_start, region_end, region_name)) %>% 
               distinct()), by = "uniprot_id") %>% 
  distinct(uniprot_id, gene_primary, sequence, length_protein, region_type, region_start, region_end, region_name)

filtered_domain_annotations <- domain_annotations %>% 
  rbind(pseudosubstrate_annotations_pubmed) %>% 
  group_by(uniprot_id) %>% 
  mutate(kinase_domain = ifelse(grepl("kinase", region_name), TRUE, FALSE)) %>%
  mutate(kinase_domain = ifelse(grepl("Kinase", region_name), TRUE, kinase_domain)) %>%
  mutate(kinase_domain = ifelse(grepl("catalytic", region_name), TRUE, kinase_domain)) %>%
  mutate(kinase_domain = max(kinase_domain)) %>% 
  filter(kinase_domain == 1) %>% 
  dplyr::select(-c(kinase_domain)) %>%
  arrange(region_start) %>% 
  mutate(region_start_no = row_number()) %>% 
  arrange(region_end) %>% 
  mutate(region_end_no = row_number()) %>% 
  ungroup() %>% 
  arrange(gene_primary, region_start_no) %>% 
  mutate(keep = ifelse(region_end_no >= region_start_no, TRUE, FALSE)) %>% 
  filter(keep == TRUE) %>% 
  group_by(uniprot_id) %>% 
  arrange(region_start) %>% 
  mutate(region_start_no = row_number()) %>% 
  arrange(region_end) %>% 
  mutate(region_end_no = row_number()) %>% 
  ungroup() %>% 
  arrange(gene_primary, region_start_no) %>% 
  mutate(keep = ifelse(region_end_no >= region_start_no, TRUE, FALSE)) %>% 
  filter(keep == TRUE) %>% 
  mutate(keep = ifelse(region_name == "Disordered", FALSE, keep)) %>% ### Exclude any disordered region, linker, or region important for any interaction
  mutate(keep = ifelse(region_name == "Linker", FALSE, keep)) %>% 
  mutate(keep = ifelse(grepl("Essential", region_name), FALSE, keep)) %>%
  mutate(keep = ifelse(grepl("Interaction", region_name), FALSE, keep)) %>%
  mutate(keep = ifelse(grepl("Important", region_name), FALSE, keep)) %>% 
  mutate(keep = ifelse(grepl("Necessary", region_name), FALSE, keep)) %>% 
  mutate(keep = ifelse(grepl("N-terminal", region_name), FALSE, keep)) %>% 
  mutate(keep = ifelse(grepl("Required", region_name), FALSE, keep)) %>% 
  mutate(keep = ifelse(grepl("Sufficient", region_name), FALSE, keep)) %>%
  mutate(keep = ifelse(grepl("Mediates", region_name), FALSE, keep)) %>%
  mutate(keep = ifelse(grepl("May", region_name), FALSE, keep)) %>%
  filter(keep == TRUE) %>% 
  dplyr::select(-c(region_start_no, region_end_no))

# Find overlapping domains and group them

to_be_concatenated_domain_annotations <- filtered_domain_annotations %>% 
  arrange(region_start) %>% 
  group_by (uniprot_id) %>% 
  mutate(domain_number = row_number()) %>% 
  ungroup() %>% 
  pivot_longer(region_start:region_end, names_to = "region_annotation", values_to = "aa_number") %>% 
  group_by(uniprot_id) %>% 
  arrange(aa_number) %>% 
  mutate(region_annotation_no = row_number()) %>% 
  ungroup() %>% 
  dplyr::select(-c(aa_number)) %>% 
  pivot_wider(names_from = c(region_annotation), values_from = c(region_annotation_no)) %>% 
  mutate(corr_domain = ifelse(region_end == region_start + 1, TRUE, FALSE)) %>% 
  left_join(filtered_domain_annotations %>% 
              dplyr::select(uniprot_id, region_name, region_start, region_end) %>% 
              dplyr::rename(region_start_aa = region_start,
                     region_end_aa = region_end), 
            by = c("uniprot_id", "region_name"))

overlapping_domains <- to_be_concatenated_domain_annotations %>% 
  filter(corr_domain == FALSE) %>% 
  arrange(gene_primary)

domain_numbers <- overlapping_domains %>% 
  pull(domain_number)

group <- cumsum(c(1, (domain_numbers[-1] - domain_numbers[-length(domain_numbers)]) != 1))

overlapping_domains$group <- group

# If regulatory domains overlap, they are concatenated into one domain. If the PK domain overlaps with other domains, the domains
# are split so that the kinase domain is undistrubed, and the regulatory domain does not overlap with the kinase domain

concatenated_domains <- overlapping_domains %>% 
  group_by(group) %>% 
  mutate(new_region_name = ifelse(grepl("rotein kinase", region_name), 1, 0)) %>%
  mutate(new_region_name = max(new_region_name) + new_region_name) %>% 
  mutate(new_domain_start = ifelse(new_region_name >= 1, region_start_aa, min(region_start_aa))) %>% 
  mutate(new_domain_start = ifelse(domain_number == max(domain_number) & new_region_name == 1, min(region_end_aa), new_domain_start)) %>% 
  mutate(new_domain_end = ifelse(new_region_name >= 1, region_end_aa, max(region_end_aa))) %>% 
  mutate(new_domain_end = ifelse(domain_number == min(domain_number) & new_region_name == 1, max(region_start_aa), new_domain_end)) %>% 
  mutate(new_region_name = ifelse(new_region_name >= 1, region_name, "Regulatory domain")) %>% 
  ungroup() %>% 
  dplyr::select(-c(region_name, region_start_aa, region_end_aa)) %>% 
  dplyr::rename(region_name = new_region_name,
         region_start_aa = new_domain_start,
         region_end_aa = new_domain_end) %>% 
  distinct(uniprot_id, 
           gene_primary, 
           sequence, 
           length_protein, 
           region_name, 
           region_start_aa, 
           region_end_aa)

complete_domain_annotations <- to_be_concatenated_domain_annotations %>% 
  filter(corr_domain == TRUE) %>% 
  dplyr::select(-c(keep, domain_number, region_start, region_end, corr_domain, region_type)) %>% 
  rbind(concatenated_domains) %>% 
  arrange(gene_primary, region_start_aa) %>% 
  distinct() %>% 
  mutate(seq_length = region_end_aa - region_start_aa) %>% 
  filter(seq_length > 6) %>%   
  group_by(uniprot_id) %>% 
  arrange(gene_primary, region_start_aa) %>% 
  mutate(domain_number = row_number()) %>% 
  ungroup() %>% 
  dplyr::select(-c(seq_length))
  
complete_domain_annotations %>% pull(gene_primary) %>% unique() %>% length()

# 
# multi_domain_kinases <- complete_domain_annotations %>% 
#   filter(!gene_primary %in% c("TTN", "OBSCN")) %>% # Exclude titin and obscurin, very big proteins!
#   group_by(uniprot_id) %>% 
#   filter(n() > 1) %>% 
#   ungroup()
# 
# multi_domain_kinases %>% pull(gene_primary) %>% unique() %>% length()

kinases_to_predict_first <- complete_domain_annotations %>% 
  filter(!gene_primary %in% c("TTN", "OBSCN")) %>% # Exclude titin and obscurin, very big proteins!
  filter(!uniprot_id %in% (transmembrane_proteins %>% pull(uniprot_id)))

kinases_to_predict_first %>% pull(gene_primary) %>% unique() %>% length()

# Make fasta format annotation

alphabet_vector <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")

fasta_annotation_domain <- kinases_to_predict_first %>%
  group_by(uniprot_id) %>% 
  mutate(ligand_number = max(domain_number) + 1) %>% 
  ungroup() %>% 
  mutate(domain_name_final = gsub(" ", "_", tolower(gsub("(.)(A-Z)", "\\1 \\2", region_name)))) %>%
  #mutate(fasta_annotation = paste(">sp|", uniprot_id, "|", gene_primary, "_", domain_name_final, " OS=Homo sapiens OX=9606", sep = "")) %>%
  mutate(domain_sequence = str_sub(sequence, region_start_aa, region_end_aa)) %>%
  mutate(domain_id = alphabet_vector[domain_number]) %>% 
  mutate(ligand_id = alphabet_vector[ligand_number]) %>% 
  distinct(uniprot_id, gene_primary, domain_number, domain_id, domain_name_final, domain_sequence, ligand_id) %>% 
  arrange(uniprot_id, domain_id)

# Make fasta files for all domains

index_vector_df <- fasta_annotation_domain %>%
  group_by(uniprot_id) %>%
  mutate(domain_number = max(domain_number)) %>%
  ungroup() %>%
  distinct(uniprot_id, domain_number)

index_list <- list()
index_list[index_vector_df$uniprot_id] <- index_vector_df$domain_number


for(protein_id in index_vector_df$uniprot_id){
  
  name <- paste(tolower(protein_id), "_no_ligand", sep = "")
  
  sequence_block <- c()
  
  max_domain <- as.numeric(index_list[protein_id])
  
  for(domain in 1:as.numeric(index_list[protein_id])){
    
    sequence_ids <- fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_id)
    sequences <- fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_sequence)
    
    sequence_block <- paste(sequence_block, "\n    {\n      \"protein\": {\n        \"id\": [\"", paste(sequence_ids, collapse='", "'), "\"],\n        \"sequence\": \"", paste(sequences, collapse='", "'), "\"\n      }\n    }", sep = "")
    
    if(domain < max_domain){
      sequence_block <- paste(sequence_block, ",", sep = "")
    }
    
  }
  
  test <- paste("{\n  \"name\": \"", name, "\",\n  \"sequences\": [", sequence_block, "\n  ],\n  \"modelSeeds\": [1],\n  \"dialect\": \"alphafold3\",\n  \"version\": 1\n}", sep = "")

  dir.create(paste("json_files_kinases/", tolower(protein_id), "_no_ligand", sep = ""))
  dir.create(paste("json_files_kinases/", tolower(protein_id), "_no_ligand", "/alphafold3_jsons", sep = ""))
  
  write(test, paste("json_files_kinases/", tolower(protein_id), "_no_ligand", "/alphafold3_jsons/", tolower(protein_id), "_no_ligand.json", sep = ""))

  filename <- file(paste("json_files_kinases/", tolower(protein_id), "_no_ligand/config.yaml", sep = ""))
  
  text <- c("alphafold3_databases: /cluster/project/alphafold/alphafold3\nalphafold3_models: /cluster/home/reberv/alphafold3_params\nalphafold3_docker: /cluster/project/beltrao/alphafold3-44e1fd5-lowspec.sif")
  
  writeLines(text, filename)
  
  close(filename)
  

}

for(protein_id in index_vector_df$uniprot_id){
  
  for(ligand in c("ATP", "STU")){
    
    
    name <- paste(tolower(protein_id), "_", tolower(ligand), sep = "")
    
    ligand_ids <- fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% pull(ligand_id) %>% unique()
    ligands <- ligand
    
    sequence_block <- c()
    
    for(domain in 1:as.numeric(index_list[protein_id])){
      
      sequence_ids <- fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_id)
      sequences <- fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_sequence)
      
      sequence_block <- paste(sequence_block, "\n    {\n      \"protein\": {\n        \"id\": [\"", paste(sequence_ids, collapse='", "'), "\"],\n        \"sequence\": \"", paste(sequences, collapse='", "'), "\"\n      }\n    },", sep = "")
      
    }
      
    
    test <- paste("{ \n  \"name\": \"", name, "\", \n  \"sequences\": [", sequence_block, "\n    { \n      \"ligand\": { \n        \"id\": [\"", ligand_ids, "\"], \n        \"ccdCodes\": [\"", ligands, "\"] \n      } \n    } \n  ], \n  \"modelSeeds\": [1], \n  \"dialect\": \"alphafold3\", \n  \"version\": 1 \n} ", sep = "")
    
    dir.create(paste("json_files_kinases/", tolower(protein_id), "_", tolower(ligand), sep = ""))
    dir.create(paste("json_files_kinases/", tolower(protein_id), "_", tolower(ligand), "/alphafold3_jsons", sep = ""))
    
    
    write(test, paste("json_files_kinases/", tolower(protein_id), "_", tolower(ligand), "/alphafold3_jsons/", tolower(protein_id), "_", tolower(ligand), ".json", sep = ""))
    
    
    filename <- file(paste("json_files_kinases/", tolower(protein_id), "_", tolower(ligand), "/config.yaml", sep = ""))

      text <- c("alphafold3_databases: /cluster/project/alphafold/alphafold3\nalphafold3_models: /cluster/home/reberv/alphafold3_params\nalphafold3_docker: /cluster/project/beltrao/alphafold3-44e1fd5-lowspec.sif")

      writeLines(text, filename)

      close(filename)


  }
  
}
