### Calculate QC parameters ###

get_qcs <- function(data, 
                    protein_id_column,
                    filename_column,
                    path) {
  
  protein_id_vector <- data %>% 
    pull({{protein_id_column}}) %>% 
    unique()
  
  output_data <- data.frame("uniprot_id" = character(0), 
                            "filename" = character(0),
                            "change_in_median_pae" = numeric(0),
                            "fraction_disordered" = numeric(0),
                            "has_clash" = numeric(0), 
                            "iptm" = numeric(0), 
                            "ptm" = numeric(0), 
                            "ranking_score" = numeric(0),
                            "pae_to_all_domains" = numeric(0),
                            "reference_domain" = numeric(0),
                            "target_domain" = numeric(0))
  
  for(protein_id in protein_id_vector){
    
    print(protein_id)
    
    file_name_with_ligand <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      filter(!grepl("ligand", {{filename_column}})) %>% 
      pull({{filename_column}})
    
    levels <- c(file_name_with_ligand, paste(protein_id, "_no_ligand.pdb", sep = ""))
    
    files <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      dplyr::mutate(sort_by_col = factor({{filename_column}}, levels = levels)) %>% 
      arrange(sort_by_col) %>% 
      dplyr::mutate(json_file_names = str_replace( {{ filename_column }} , "model.pdb", "summary_confidences.json")) %>% 
      pull(json_file_names)
    
    # read in two json files
    
    json_1 <- fromJSON(file = paste(path, "/all_summary_confidences/", files[1], sep = ""))
    
    json_2 <- fromJSON(file = paste(path, "/all_summary_confidences/", files[2], sep = ""))
    

    
    # Use median PAE
    
    files <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      dplyr::mutate(sort_by_col = factor({{filename_column}}, levels = levels)) %>% 
      arrange(sort_by_col) %>% 
      dplyr::mutate(json_file_names = str_replace( {{ filename_column }} , "model.pdb", "confidences.json")) %>% 
      pull(json_file_names)
    
    # read in two json files
    
    json_3 <- fromJSON(file = paste(path, "/all_confidences/", files[1], sep = ""))
    
    json_4 <- fromJSON(file = paste(path, "/all_confidences/", files[2], sep = ""))
    
    chains_json_4 <- json_4[["token_chain_ids"]]
    
    all_chains <- chains_json_4 %>% unique()
    
    json_3_paes <- json_3[["pae"]]
    
    json_4_paes <- json_4[["pae"]]
    
    json_3_median_pae_matrix <- matrix(nrow = length(all_chains), ncol = length(all_chains))
    rownames(json_3_median_pae_matrix) <- all_chains
    colnames(json_3_median_pae_matrix) <- all_chains
    
    json_4_median_pae_matrix <- matrix(nrow = length(all_chains), ncol = length(all_chains))
    rownames(json_4_median_pae_matrix) <- all_chains
    colnames(json_4_median_pae_matrix) <- all_chains
    
    for(chain in all_chains){
      
      json_3_chain_pae <- c()
      json_4_chain_pae <- c()
      
      json_3_current_chain <- json_3_paes[which(chains_json_4 == chain)]
      json_4_current_chain <- json_4_paes[which(chains_json_4 == chain)]
      
      for(chain_2 in all_chains){
        
        for(item in 1:length(json_4_current_chain)){
          
          json_3_chain_pae <- append(json_3_chain_pae, json_3_current_chain[[item]][which(chains_json_4 == chain_2)])
          json_4_chain_pae <- append(json_4_chain_pae, json_4_current_chain[[item]][which(chains_json_4 == chain_2)])
          
          
        }
        
        json_3_median_pae_matrix[chain, chain_2] <- median(json_3_chain_pae)
        json_4_median_pae_matrix[chain, chain_2] <- median(json_4_chain_pae)
        
      }
      
    }
    
    abs_change_in_pae <- abs(json_3_median_pae_matrix - json_4_median_pae_matrix) %>% mean()    
    
    
    
    catalytic_domain <- fasta_annotation_domain %>%
      filter(uniprot_id == toupper(protein_id)) %>%
      filter(
        domain_name_final == "kinase" |
          grepl("protein_kinase", domain_name_final) |
          grepl("pseudokinase", domain_name_final) |
          grepl("histidine_kinase", domain_name_final) |
          grepl("catalytic", domain_name_final)
      ) %>%
      pull(domain_number)
    
    for(ref_domain in c(catalytic_domain)){
      
      tmp_pae <- json_1[["chain_pair_pae_min"]]

       
    output_data <- output_data %>%
      rbind(data.frame("uniprot_id" = protein_id,
                       "filename" = files[1],
                       "change_in_median_pae" = abs_change_in_pae,
                       "fraction_disordered" = json_1[["fraction_disordered"]],
                       "has_clash" = json_1[["has_clash"]],
                       "iptm" = json_1[["iptm"]],
                       "ptm" = json_1[["ptm"]],
                       "ranking_score" = json_1[["ranking_score"]],
                       "pae_to_all_domains" = tmp_pae[[ref_domain]],
                       "reference_domain" = ref_domain) %>% 
              mutate(target_domain = row_number()))
    
    tmp_pae <- json_2[["chain_pair_pae_min"]]
    
    output_data <- output_data %>%
      rbind(data.frame("uniprot_id" = protein_id,
                       "filename" = files[2],
                       "change_in_median_pae" = NA,
                       "fraction_disordered" = json_2[["fraction_disordered"]],
                       "has_clash" = json_2[["has_clash"]],
                       "iptm" = ifelse(is.null(json_2[["iptm"]]), NA, json_2[["iptm"]]) ,
                       "ptm" = json_2[["ptm"]],
                       "ranking_score" = json_2[["ranking_score"]],
                       "pae_to_all_domains" = tmp_pae[[ref_domain]],
                       "reference_domain" = ref_domain) %>% 
              mutate(target_domain = row_number()))
      
    }
   
    
  }
  
  return(output_data)
  
}


### Calculate prediction error ###

calculate_reg_domain_prediction_error <- function(data, 
                                                  protein_id_column,
                                                  filename_column,
                                                  path) {
  
  
  protein_id_vector <- data %>% 
    pull({{protein_id_column}}) %>% 
    unique()
  
  output_data <- data.frame("uniprot_id" = character(0), 
                            "kinase_domain_filename" = character(0),
                            "rmsd" = numeric(0), 
                            "rmssd_cat" = numeric(0), 
                            "rmssd" = numeric(0), 
                            "pae_weighted_rmssd" = numeric(0),
                            "max_distance" = numeric(0), 
                            "no_changing_aa" = numeric(0), 
                            "no_aligned_residues" = numeric(0))
  
  for(protein_id in protein_id_vector){
    
    print(protein_id)
    
    file_name_with_ligand <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      filter(!grepl("ligand", {{filename_column}})) %>% 
      pull({{filename_column}})
    
    levels <- c(file_name_with_ligand, paste(protein_id, "_no_ligand.pdb", sep = ""))
    
    files <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      dplyr::mutate(sort_by_col = factor({{filename_column}}, levels = levels)) %>% 
      arrange(sort_by_col) %>% 
      pull( {{filename_column}} )
    
    # read in two PDB files
    
    pdb_1 <- read.pdb(file = paste(path, "/", files[1], sep = ""))
    pdb_2 <- read.pdb(file = paste(path, "/", files[2], sep = ""))
    
    catalytic_domain <- fasta_annotation_domain %>% 
      filter(uniprot_id == toupper(protein_id)) %>% 
      filter( domain_name_final == "kinase" | grepl("protein_kinase", domain_name_final) | grepl("pseudokinase", domain_name_final) | grepl("histidine_kinase", domain_name_final) | grepl("catalytic", domain_name_final)) %>% 
      pull(domain_id)
    
    
    high_confidence_domains <- fasta_annotation_domain %>% 
      filter(uniprot_id == toupper(protein_id)) %>% 
      filter(domain_number %in% (qcs_stu %>% 
                                   filter(uniprot_id == protein_id) %>% 
                                   filter(pae_to_all_domains < 10) %>% 
                                   pull(target_domain) %>% 
                                   unique())) %>% 
      filter(!domain_id %in% catalytic_domain) %>% 
      pull(domain_id)
      
    # align the files
    alignment <- struct.aln(fixed = pdb_1, 
                            mobile = pdb_2, 
                            fixed.inds = atom.select(pdb_1, "protein", chain=catalytic_domain), 
                            mobile.inds = atom.select(pdb_2, "protein", chain=catalytic_domain), 
                            #max.cycles = 0,
                            exefile = "msa")
       
    
    
    coordinate <- c()
    
    for(i in 1:length(alignment$xyz)){
      
      coordinate <- append(coordinate, pdb_1$xyz[[i]])
      coordinate <- append(coordinate, alignment$xyz[[i]])
      
    }
    
    
    reference_paes <- qcs_stu %>% 
      filter(uniprot_id == protein_id) %>% 
      group_by(target_domain) %>% 
      mutate(diff_pae = max(pae_to_all_domains) - min(pae_to_all_domains)) %>% 
      mutate(pae = min(pae_to_all_domains)) %>% 
      ungroup() %>% 
      dplyr::distinct(target_domain, pae, diff_pae) %>% 
      mutate(chain = alphabet_vector[target_domain]) %>% 
      dplyr::select(-c(target_domain))
    
    
    aligned_coordinates <- data.frame("coordinate" = coordinate) %>% 
      mutate(index = row_number() - 1) %>% 
      mutate(dimension = ifelse(floor((index %% 6) / 2) == 0, "x", "y")) %>% 
      mutate(dimension = ifelse(floor((index %% 6) / 2) == 2, "z", dimension)) %>% 
      mutate(amino_acid = floor(index/6) + 1) %>% 
      mutate(pdb_id = index %% 2 + 1) %>% 
      distinct(amino_acid, pdb_id, dimension, coordinate) %>% 
      pivot_wider(names_from = c(dimension, pdb_id), values_from = coordinate) %>% 
      mutate(b_1 = pdb_1$atom$b[c(1:length(pdb_2$atom$b))] ) %>% 
      mutate(b_2 = pdb_2$atom$b) %>% 
      mutate(chain = pdb_2$atom$chain) %>% 
      left_join(reference_paes, by = c("chain")) %>% 
      mutate(distance_squared = (x_1 - x_2) ^ 2 + (y_1 - y_2) ^ 2 + (z_1 - z_2) ^ 2) %>% 
      group_by(amino_acid) %>% 
      mutate(flexibility_score_cat_domain = ifelse((max(b_1, b_2) < 70), 0, 1)) %>% 
      mutate(flexibility_score_cat_domain = ifelse(!chain %in% catalytic_domain, 0, flexibility_score_cat_domain)) %>% 
      #mutate(flexibility_score_cat_domain = ((min(b_1, b_2) / 100) ** 2)) %>% 
      mutate(flexibility_score_reg_and_cat_domain = ifelse((max(b_1, b_2) < 70), 0, 1)) %>% 
      mutate(flexibility_score_reg_and_cat_domain = ifelse(!chain %in% c(high_confidence_domains, catalytic_domain), 0, flexibility_score_reg_and_cat_domain)) %>% 
      mutate(flexibility_score_reg_and_cat_domain = ifelse(chain %in% high_confidence_domains, 1, flexibility_score_reg_and_cat_domain)) %>% 
      mutate(flexibility_score_pae_weighted = ifelse((max(b_1, b_2) < 70), 0, 1)) %>% 
      mutate(flexibility_score_pae_weighted = ifelse(!chain %in% catalytic_domain, min((40 - pae) / 30, 1), flexibility_score_pae_weighted)) %>% 
      #mutate(flexibility_score_reg_and_cat_domain = ((min(b_1, b_2) / 100) ** 2)) %>% 
      mutate(scaled_square_distance_cat_domain = distance_squared * flexibility_score_cat_domain * flexibility_score_cat_domain) %>% 
      mutate(scaled_square_distance_reg_and_cat_domain = distance_squared * flexibility_score_reg_and_cat_domain * flexibility_score_reg_and_cat_domain) %>% 
      mutate(scaled_square_distance_pae_weighted = (distance_squared) * flexibility_score_pae_weighted * flexibility_score_pae_weighted) %>% 
      ungroup()
    

    #aligned_coordinates %>% pull(scaled_square_distance) %>% mean() %>% sqrt()
    
    number_aa_with_dist_1 <- aligned_coordinates %>%
      filter(distance_squared > 1) %>%
      dplyr::pull(distance_squared) %>%
      length()
    
    max_distance <- aligned_coordinates %>% 
      filter(!is.na(distance_squared)) %>% 
      pull(distance_squared) %>% 
      max()
    
    no_aligned_residues <- aligned_coordinates %>% 
      filter(!is.na(distance_squared)) %>% 
      pull(distance_squared) %>% 
      length()
    
    
    # RMSD and max distance
    
    output_data <- output_data %>% rbind(data.frame("uniprot_id" = protein_id, 
                                                    "kinase_domain_filename" = files[1],
                                                    "rmsd" = rmsd(pdb_1$xyz, alignment$xyz), 
                                                    "rmssd_cat" = aligned_coordinates %>% pull(scaled_square_distance_cat_domain) %>% mean() %>% sqrt(),
                                                    "rmssd" = aligned_coordinates %>% pull(scaled_square_distance_reg_and_cat_domain) %>% mean() %>% sqrt(), 
                                                    "pae_weighted_rmssd" = aligned_coordinates %>% pull(scaled_square_distance_pae_weighted) %>% mean() %>% sqrt(),
                                                    "max_distance" = sqrt(max_distance), 
                                                    "no_changing_aa" = number_aa_with_dist_1,
                                                    "no_aligned_residues" = no_aligned_residues))
    
    
  }
  
  return(output_data)
  
}



### Calculate SASA ###

calculate_change_in_surface_accessibility <- function(data,
                                                      protein_id_column,
                                                      filename_column,
                                                      all_domain_annotations,
                                                      path) {
  
  
  protein_id_vector <- data %>% 
    pull({{protein_id_column}}) %>% 
    unique()
  
  output_data <- data.frame("uniprot_id" = character(0), 
                            "kinase_domain_file" = character(0),
                            "change_in_main_chain_surface" = numeric(0), 
                            "change_in_side_chain_surface" = numeric(0), 
                            "change_in_residue_surface" = numeric(0))
  
  for(protein_id in protein_id_vector){
    
    print(protein_id)
    
    file_name_with_ligand <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      filter(!grepl("ligand", {{filename_column}})) %>% 
      pull({{filename_column}})
    
    levels <- c(file_name_with_ligand, paste(protein_id, "_no_ligand.pdb", sep = ""))
    
    files <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      dplyr::mutate(sort_by_col = factor({{filename_column}}, levels = levels)) %>% 
      arrange(sort_by_col) %>% 
      dplyr::mutate(csv_file_name = str_replace( {{ filename_column }} , ".pdb", ".csv")) %>% 
      pull(csv_file_name)
    
    
    # read in two csv files
    
    csv_1 <- read.csv(file = paste(path, "/", files[1], sep = "")) %>% 
      dplyr::select(-c(X))
    
    csv_2 <- read.csv(file = paste(path, "/", files[2], sep = "")) %>% 
      dplyr::select(-c(X))
    
    
    
    # align the files
    
    surface_comparison <- csv_1 %>% 
      dplyr::rename(chain_index_1 = chain_index, 
                    main_chain_surface_1 = main_chain_surface,
                    side_chain_surface_1 = side_chain_surface) %>% 
      mutate(kinase_domain_file = files[1]) %>% 
      left_join((csv_2 %>% 
                   dplyr::rename(chain_index_2 = chain_index, 
                                 main_chain_surface_2 = main_chain_surface,
                                 side_chain_surface_2 = side_chain_surface)),
                by = c("aa_index", "aa"), relationship = "many-to-many") %>% 
      mutate(change_in_main_chain_surface = main_chain_surface_2 - main_chain_surface_1) %>% 
      mutate(change_in_side_chain_surface = side_chain_surface_2 - side_chain_surface_1) %>% 
      mutate(change_in_residue_surface = change_in_main_chain_surface + change_in_side_chain_surface)
    
    # get binding site
    
    binding_site_name <- data %>% 
      arrange({{protein_id_column}}, {{filename_column}}) %>% 
      dplyr::filter( {{ protein_id_column }} == protein_id) %>% 
      dplyr::filter( !grepl( "ligand", {{filename_column}}) ) %>% 
      dplyr::mutate(binding_site_name = str_replace( {{ filename_column }} , "model.pdb", "binding_site.pdb")) %>% 
      pull(binding_site_name)
    
    pdb_binding_site <- read.pdb(file = paste(path, "/binding_sites/", binding_site_name, sep = ""))
    
    binding_site_residues <- pdb_binding_site$atom$resno
    
    catalytic_domain <- fasta_annotation_domain %>% 
      filter(uniprot_id == toupper(protein_id)) %>% 
      filter(
        domain_name_final == "kinase" |
          grepl("protein_kinase", domain_name_final) |
          grepl("pseudokinase", domain_name_final) |
          grepl("histidine_kinase", domain_name_final) |
          grepl("catalytic", domain_name_final)
      ) %>%
      pull(domain_id)
    
    
    
    change_in_surface_at_ATP_binding_site <- surface_comparison %>% 
      filter(aa_index %in% binding_site_residues) %>% 
      filter(chain_index_1 %in% catalytic_domain) %>% 
      mutate(change_in_main_chain_surface = sum(change_in_main_chain_surface)) %>% 
      mutate(change_in_side_chain_surface = sum(change_in_side_chain_surface)) %>% 
      mutate(change_in_residue_surface = sum(change_in_residue_surface)) %>% 
      distinct(kinase_domain_file, change_in_main_chain_surface, change_in_side_chain_surface, change_in_residue_surface)
    
    
    output_data <- output_data %>% rbind(change_in_surface_at_ATP_binding_site %>% 
                                           mutate(uniprot_id = protein_id) %>% 
                                           distinct(uniprot_id,
                                                    kinase_domain_file,
                                                    change_in_main_chain_surface,
                                                    change_in_side_chain_surface,
                                                    change_in_residue_surface))
    
  }
  
  return(output_data)
  
}

