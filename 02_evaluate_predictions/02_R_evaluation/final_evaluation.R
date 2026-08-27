library(protti)
library(bio3d)
library(cry)
library(dplyr)
library(msa)
library(stringr)


setwd("Z:/Viviane/Computational/AF3/")

# Find an index completed predictions, rerun predictions if needed

files <- list.files(path = "Z:/Viviane/Computational/AF3/all_top_models/")

ranked_0_files <- data.frame(filename = files) %>% 
  filter(grepl("*\\.pdb", filename)) %>% 
  mutate(uniprot_id = str_extract(filename, "[a-z0-9]*")) %>% 
  mutate(ligand = str_extract(filename, "\\_[a-z]*\\_(ligand\\_)?"))

complete_files <- data.frame(uniprot_id = c(rep(tolower(all_kinases_domain_annotations_for_input_json_files %>% pull(uniprot_id) %>% unique()), 2)),
                             ligand = c(rep("_stu_", 386), rep("_no_ligand_", 386))) %>% # 386 because that is how many models were generated
  left_join(ranked_0_files, by = c("uniprot_id", "ligand"))

usable_kinases_stu <- complete_files %>% 
  filter(!is.na(filename)) %>%
  filter(ligand %in% c("_stu_", "_no_ligand_")) %>% 
  group_by(uniprot_id) %>% 
  filter(n() == 2) %>% 
  ungroup()


# Calculate QC parameters, filter out any proteins if needed

qcs_stu <- get_qcs(data = usable_kinases_stu, 
                   protein_id_column = uniprot_id,
                   filename_column = filename,
                   path = "Z:/Viviane/Computational/AF3")

# Calculate RMSSD

prediction_error_stu <- calculate_RMSSD(data = usable_kinases_stu,
                                        protein_id_column = uniprot_id,
                                        filename_column = filename,
                                        path = "Z:/Viviane/Computational/AF3/all_top_models")

# Calculate dSASA

change_in_sasa_stu <- calculate_dSASA(data = usable_kinases_stu,
                                      protein_id_column = uniprot_id,
                                      filename_column = filename,
                                      all_domain_annotations = kinases_to_predict_first,
                                      path = "Z:/Viviane/Computational/AF3/all_top_models")

# Final output

rmsd_and_surface_stu <- prediction_error_stu %>%
  mutate(kinase_domain_filename = str_replace(kinase_domain_filename, "\\_model.pdb", "")) %>% 
  left_join((change_in_sasa_stu %>% 
               mutate(kinase_domain_file = str_replace(kinase_domain_file, "\\_model.csv", ""))), 
            by = c("uniprot_id", "kinase_domain_filename" = "kinase_domain_file"))