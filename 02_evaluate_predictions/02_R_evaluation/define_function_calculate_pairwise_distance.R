# This function calculates the distance changes between residues in apo and staurosporine bound structures


calculate_pairwise_distance <- function(data, 
                                        protein_id_column,
                                        filename_column,
                                        path) {
  
  
  protein_id_vector <- data %>% 
    pull({{protein_id_column}}) %>% 
    unique()
  
  print(protein_id_vector)
  
  output_data <- data.frame("protein_id" = character(0), 
                            "chain" = character(0), 
                            "resno" = numeric(0), 
                            "mean_distance_squared" = numeric(0))
  
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
      mutate(resno = pdb_2$atom$resno) %>% 
      left_join(reference_paes, by = c("chain")) %>% 
      mutate(distance_squared = (x_1 - x_2) ^ 2 + (y_1 - y_2) ^ 2 + (z_1 - z_2) ^ 2) %>% 
      group_by(chain, resno) %>% 
      mutate(mean_distance_squared = mean(distance_squared)) %>% 
      mutate(mean_b_1 = mean(b_1)) %>% 
      mutate(mean_b_2 = mean(b_2)) %>% 
      ungroup() %>%       
      mutate(protein_id = protein_id) %>% 
      distinct(protein_id, chain, resno, mean_distance_squared, mean_b_1, mean_b_2, pae)
    
    output_data <- output_data %>% rbind(aligned_coordinates)
    
    
  }
  
  return(output_data)
  
}