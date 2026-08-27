# This script takes the domain annotations from the file "all_kinases_domain_annotations_for_input_json_files.csv"
# and generates the .json and config.yaml files in the correct folder structure to predict using jurgjn/batch-infer
# on the ETH Euler cluster. 



fasta_annotation_domain <- read.csv("all_kinases_domain_annotations_for_input_json_files.csv")

index_vector_df <- fasta_annotation_domain %>%
  group_by(uniprot_id) %>%
  mutate(domain_number = max(domain_number)) %>%
  ungroup() %>%
  distinct(uniprot_id, domain_number)

index_list <- list()
index_list[index_vector_df$uniprot_id] <- index_vector_df$domain_number


for(protein_id in index_vector_df$uniprot_id) {
  name <- paste(tolower(protein_id), "_no_ligand", sep = "")
  
  sequence_block <- c()
  
  max_domain <- as.numeric(index_list[protein_id])
  
  for (domain in 1:as.numeric(index_list[protein_id])) {
    sequence_ids <-
      fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_id)
    sequences <-
      fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_sequence)
    
    sequence_block <-
      paste(
        sequence_block,
        "\n    {\n      \"protein\": {\n        \"id\": [\"",
        paste(sequence_ids, collapse = '", "'),
        "\"],\n        \"sequence\": \"",
        paste(sequences, collapse = '", "'),
        "\"\n      }\n    }",
        sep = ""
      )
    
    if (domain < max_domain) {
      sequence_block <- paste(sequence_block, ",", sep = "")
    }
    
  }
  
  test <-
    paste(
      "{\n  \"name\": \"",
      name,
      "\",\n  \"sequences\": [",
      sequence_block,
      "\n  ],\n  \"modelSeeds\": [1],\n  \"dialect\": \"alphafold3\",\n  \"version\": 1\n}",
      sep = ""
    )
  
  dir.create(
    paste(
      "Z:/Viviane/Computational/json_files_kinases/",
      tolower(protein_id),
      "_no_ligand",
      sep = ""
    )
  )
  dir.create(
    paste(
      "Z:/Viviane/Computational/json_files_kinases/",
      tolower(protein_id),
      "_no_ligand",
      "/alphafold3_jsons",
      sep = ""
    )
  )
  
  write(
    test,
    paste(
      "Z:/Viviane/Computational/json_files_kinases/",
      tolower(protein_id),
      "_no_ligand",
      "/alphafold3_jsons/",
      tolower(protein_id),
      "_no_ligand.json",
      sep = ""
    )
  )
  
  filename <-
    file(
      paste(
        "Z:/Viviane/Computational/json_files_kinases/",
        tolower(protein_id),
        "_no_ligand/config.yaml",
        sep = ""
      )
    )
  
  text <-
    c(
      "alphafold3_databases: /cluster/project/alphafold/alphafold3\nalphafold3_models: /cluster/home/reberv/alphafold3_params\nalphafold3_docker: /cluster/project/beltrao/alphafold3-44e1fd5-lowspec.sif"
    )
  
  writeLines(text, filename)
  
  close(filename)
  
  
}

for(protein_id in index_vector_df$uniprot_id){
  
  ligand = "STU"
  
  name <- paste(tolower(protein_id), "_stu", sep = "")
  
  ligand_ids <-
    fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% pull(ligand_id) %>% unique()
  ligands <- "STU"
  
  sequence_block <- c()
  
  for (domain in 1:as.numeric(index_list[protein_id])) {
    sequence_ids <-
      fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_id)
    sequences <-
      fasta_annotation_domain %>% filter(uniprot_id == protein_id) %>% filter(domain_number == domain) %>% pull(domain_sequence)
    
    sequence_block <-
      paste(
        sequence_block,
        "\n    {\n      \"protein\": {\n        \"id\": [\"",
        paste(sequence_ids, collapse = '", "'),
        "\"],\n        \"sequence\": \"",
        paste(sequences, collapse = '", "'),
        "\"\n      }\n    },",
        sep = ""
      )
    
  }
  
  
  test <-
    paste(
      "{ \n  \"name\": \"",
      name,
      "\", \n  \"sequences\": [",
      sequence_block,
      "\n    { \n      \"ligand\": { \n        \"id\": [\"",
      ligand_ids,
      "\"], \n        \"ccdCodes\": [\"",
      ligands,
      "\"] \n      } \n    } \n  ], \n  \"modelSeeds\": [1], \n  \"dialect\": \"alphafold3\", \n  \"version\": 1 \n} ",
      sep = ""
    )
  
  dir.create(paste(
    "Z:/Viviane/Computational/json_files_kinases/",
    tolower(protein_id),
    "_stu",
    sep = ""
  ))
  dir.create(
    paste(
      "Z:/Viviane/Computational/json_files_kinases/",
      tolower(protein_id),
      "_stu",
      "/alphafold3_jsons",
      sep = ""
    )
  )
  
  
  write(
    test,
    paste(
      "Z:/Viviane/Computational/json_files_kinases/",
      tolower(protein_id),
      "_stu",
      "/alphafold3_jsons/",
      tolower(protein_id),
      "_stu",
      ".json",
      sep = ""
    )
  )
  
  
  filename <-
    file(
      paste(
        "Z:/Viviane/Computational/json_files_kinases/",
        tolower(protein_id),
        "_stu",
        "/config.yaml",
        sep = ""
      )
    )
  
  text <-
    c(
      "alphafold3_databases: /cluster/project/alphafold/alphafold3\nalphafold3_models: /cluster/home/reberv/alphafold3_params\nalphafold3_docker: /cluster/project/beltrao/shared/alphafold3/images/alphafold3-v3.0.1.sif"
    )
  
  writeLines(text, filename)
  
  close(filename)
  
}