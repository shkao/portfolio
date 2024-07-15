rm(list = ls())
pacman::p_load(tidyverse, stringr, GEOmetadb, pbapply, glue)

# Config for GEOmetadb
geo_sqlfile <- file.path("~/db/GEOmetadb.sqlite")
file.info(geo_sqlfile)
geo_con <- dbConnect(SQLite(), geo_sqlfile)

#
term_1 <- "CX-5461"
term_2 <- "CX5461"
sql_str <-
  glue(
    "SELECT gse.title, gse.summary, gse.gse, gsm.gsm, gsm.gpl, gsm.organism_ch1,
    gse.type, gse.submission_date, gse.pubmed_id, gsm.characteristics_ch1
    FROM gse JOIN gsm ON gse.gse=gsm.series_id
    WHERE gsm.organism_ch1 = 'Homo sapiens'
    AND
      (gse.type = 'Expression profiling by array' OR
      gse.type = 'Expression profiling by high throughput sequencing')
    AND
      (UPPER(gse.title) LIKE UPPER('% {term_1} %') OR
      UPPER(gse.summary) LIKE UPPER('% {term_1} %') OR
      UPPER(gse.title) LIKE UPPER('% {term_2} %') OR
      UPPER(gse.summary) LIKE UPPER('% {term_2} %'))"
  )

#
res <- dbGetQuery(geo_con, sql_str)
res %>% writexl::write_xlsx("20221215_GEO_search_full.xlsx")
res %>%
  dplyr::select(-gsm, -gpl, -organism_ch1, -characteristics_ch1) %>%
  dplyr::distinct() %>%
  writexl::write_xlsx("20221215_GEO_search.xlsx")

# Save the results and close db connection
dbDisconnect(geo_con)
