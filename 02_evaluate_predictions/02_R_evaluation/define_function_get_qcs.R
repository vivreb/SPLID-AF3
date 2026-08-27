# This function extracts QC data from the PDB files of the predicted structures.

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