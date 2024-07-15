rm(list = ls())
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999)
pacman::p_load(tidyverse, depmap, ExperimentHub)

#
query_cpd <- "CX-5461"
query_genes <- c("MYC", "USP7")

# metadata
metadata_file <- "../data/metadata.Rdata"
if (!file.exists(metadata_file)) {
  df_cellline <-
    read_csv("sample_info.csv.gz") %>%
    select(DepMap_ID, cell_line_name, primary_disease, Subtype) %>%
    select(-primary_disease, -Subtype)

  df_oncotree <-
    read_csv("public_22q2_new_annotations.csv.gz") %>%
    select(DepMap_ID, OncotreePrimaryDisease, OncotreeSubtype) %>%
    mutate(
      full_disease_name = paste0(OncotreePrimaryDisease, " - ", OncotreeSubtype)
    )

  df_metadata <-
    left_join(df_cellline, df_oncotree, by = "DepMap_ID") %>%
    rename(depmap_id = DepMap_ID)

  save(df_cellline, df_oncotree, df_metadata, file = metadata_file)
}

# depmap
depmap_file <- "../data/depmap.Rdata"
if (!file.exists(depmap_file)) {
  # Load gene expression data
  df_exp <-
    read_csv("Expression_Public_23Q2.csv.gz") %>%
    rename("depmap_id" = "...1") %>%
    select(c("depmap_id", all_of(query_genes)))

  # Load copy number data
  df_copynumber <-
    read_csv("Copy_Number_Public_23Q2.csv.gz") %>%
    rename("depmap_id" = "...1") %>%
    select(c("depmap_id", all_of(query_genes)))

  # Load mutation data
  df_mut <-
    read_csv("Damaging_Mutations.csv.gz") %>%
    rename("depmap_id" = "...1") %>%
    select(c("depmap_id", all_of(query_genes)))

  # Drug sensitivity: PRISM
  df_drug_sen_prism <-
    read_csv(
      "Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4.csv.gz",
      show_col_types = FALSE
    ) %>%
    rename("depmap_id" = "...1")
  colnames(df_drug_sen_prism) <-
    sapply(strsplit(colnames(df_drug_sen_prism), " "), head, 1)

  # Drug sensitivity: PRISM AUC
  df_drug_sen_prism_auc <-
    read_csv(
      "Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4.csv.gz",
      show_col_types = FALSE
    ) %>%
    rename("depmap_id" = "...1")
  colnames(df_drug_sen_prism_auc) <-
    sapply(strsplit(colnames(df_drug_sen_prism_auc), " "), head, 1)

  # Drug sensitivity: GDSC AUC
  df_drug_sen_gdsc_auc <-
    read_csv("Drug_sensitivity_AUC_(Sanger_GDSC1).csv.gz",
                    show_col_types = FALSE) %>%
    rename("depmap_id" = "...1")
  colnames(df_drug_sen_gdsc_auc) <-
    sapply(strsplit(colnames(df_drug_sen_gdsc_auc), " "), head, 1)

  save(
    df_exp,
    df_copynumber,
    df_mut,
    df_drug_sen_prism,
    df_drug_sen_prism_auc,
    df_drug_sen_gdsc_auc,
    file = depmap_file
  )
}
