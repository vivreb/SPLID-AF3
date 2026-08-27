# This function calculates changes in SASA at the ligand binding site

calculate_dSASA <- function(data,
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