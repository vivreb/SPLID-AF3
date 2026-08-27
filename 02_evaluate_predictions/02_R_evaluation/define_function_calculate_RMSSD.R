# This function calculates RMSSD for apo and staurosporine bound predictions

calculate_RMSSD <- function(data, 
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
      mutate(flexibility_score_reg_and_cat_domain = ifelse((max(b_1, b_2) < 70), 0, 1)) %>% 
      mutate(flexibility_score_reg_and_cat_domain = ifelse(!chain %in% c(high_confidence_domains, catalytic_domain), 0, flexibility_score_reg_and_cat_domain)) %>% 
      mutate(flexibility_score_reg_and_cat_domain = ifelse(chain %in% high_confidence_domains, 1, flexibility_score_reg_and_cat_domain)) %>% 
      mutate(flexibility_score_pae_weighted = ifelse((max(b_1, b_2) < 70), 0, 1)) %>% 
      mutate(flexibility_score_pae_weighted = ifelse(!chain %in% catalytic_domain, min((31.75 - pae) / 21.75, 1), flexibility_score_pae_weighted)) %>% 
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