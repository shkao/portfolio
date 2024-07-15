rm(list = ls())
library(tidyverse)
Sys.setenv("VROOM_CONNECTION_SIZE" = 999999)

#
df_info <-
  read_csv("public_22q2_new_annotations.csv.gz", show_col_types = FALSE)
head(df_info)

#
df_prism <-
  read_csv(
    "Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4.csv.gz",
    show_col_types = FALSE
  ) %>% rename("DepMap_ID" = "...1")
df_prism_count <- data.frame(
  Dataset = "PRISM Primary",
  Primary = df_info %>%
    filter(DepMap_ID %in% df_prism$DepMap_ID) %>%
    pull(`OncotreePrimaryDisease`) %>%
    unique() %>%
    length(),
  Subtype = df_info %>%
    filter(DepMap_ID %in% df_prism$DepMap_ID) %>%
    pull(OncotreeSubtype) %>%
    unique() %>%
    length()
)
head(df_prism_count)

#
df_prism_auc <-
  read_csv(
    "Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4.csv.gz",
    show_col_types = FALSE
  ) %>% rename("DepMap_ID" = "...1")
df_prism_auc_count <- data.frame(
  Dataset = "PRISM AUC",
  Primary = df_info %>%
    filter(DepMap_ID %in% df_prism_auc$DepMap_ID) %>%
    pull(OncotreePrimaryDisease) %>%
    unique() %>%
    length(),
  Subtype = df_info %>%
    filter(DepMap_ID %in% df_prism_auc$DepMap_ID) %>%
    pull(OncotreeSubtype) %>%
    unique() %>%
    length()
)
head(df_prism_auc_count)

#
df_gdsc1 <-
  read_csv("Drug_sensitivity_AUC_(Sanger_GDSC1).csv.gz",
                  show_col_types = FALSE) %>%
  rename("DepMap_ID" = "...1")
df_gdsc1_count <- data.frame(
  Dataset = "GDSC1",
  Primary = df_info %>%
    filter(DepMap_ID %in% df_gdsc1$DepMap_ID) %>%
    pull(OncotreePrimaryDisease) %>%
    unique() %>%
    length(),
  Subtype = df_info %>%
    filter(DepMap_ID %in% df_gdsc1$DepMap_ID) %>%
    pull(OncotreeSubtype) %>%
    unique() %>%
    length()
)

#
df_all <- rbind(df_prism_count,
                df_prism_auc_count,
                df_gdsc1_count)

write_csv(df_all, "disease_cov.csv")
