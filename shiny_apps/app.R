pacman::p_load(
  Cairo,
  shiny,
  shinydashboard,
  tidyverse,
  ggplot2,
  ggpubr,
  ggthemes,
  ggvenn,
  shinycssloaders
)

load("data/metadata.Rdata")
load("data/depmap.Rdata")

#
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Analysis for Senhwa Biosciences"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(fluidRow(box(
    width = 12,
    column(
      4,
      selectizeInput("query_cpd",
                     label = "Compound",
                     choices = NULL)
    ),
    column(4,
           selectizeInput(
             "query_gene",
             label = "Gene",
             choices = NULL
           )),
    column(
      4,
      selectizeInput(
        "query_gene_feature",
        label = "Gene feature",
        choices = c("Expression", "Copy number")
      )
    )
  )),
  fluidRow(column(
    6,
    box(
      width = "100%",
      height = "auto",
      DT::DTOutput("tbl_cor", height = "100%") %>% withSpinner()
    )
  ),
  column(
    6,
    tabBox(
      width = "100%",
      height = "auto",
      tabPanel(
        "Scatter plot",
        NULL,
        plotOutput("fig_scatter_out") %>% withSpinner(),
        downloadLink("download_fig", "Download figure"),
        plotOutput("fig_venn_out") %>% withSpinner()
      ),
      tabPanel(
        "Data",
        NULL,
        DT::DTOutput("tbl_scatter", height = "100%") %>% withSpinner(),
        downloadButton("download_tab", "Download table")
      )
    )
  )))
)

server <- function(input, output, session) {
  updateSelectizeInput(
    session,
    "query_cpd",
    choices = levels(factor(colnames(
      df_drug_sen_prism_auc[, -1]
    ))),
    selected = "CX-5461",
    server = TRUE
  )

  updateSelectizeInput(
    session,
    "query_gene",
    choices = levels(factor(colnames(df_exp[, -1]))),
    selected = "MYC",
    server = TRUE
  )

  df_drug_sen_prism_queried <- reactive({
    req(input$query_cpd)
    df <- df_drug_sen_prism %>%
      dplyr::select(depmap_id, contains(input$query_cpd))

    if (ncol(df) == 2) {
      df %>% dplyr::rename("PRISM" = input$query_cpd)
    } else {
      df$`PRISM` <- NA
      df
    }
  })

  df_drug_sen_prism_auc_queried <- reactive({
    req(input$query_cpd)
    df <- df_drug_sen_prism_auc %>%
      dplyr::select(depmap_id, contains(input$query_cpd))

    if (ncol(df) == 2) {
      df %>% dplyr::rename("PRISM AUC" = input$query_cpd)
    } else {
      df$`PRISM AUC` <- NA
      df
    }
  })

  df_drug_sen_gdsc_auc_queried <- reactive({
    req(input$query_cpd)
    df <- df_drug_sen_gdsc_auc %>%
      dplyr::select(depmap_id, contains(input$query_cpd))

    if (ncol(df) == 2) {
      df %>% dplyr::rename("GDSC AUC" = input$query_cpd)
    } else {
      df$`GDSC AUC` <- NA
      df
    }
  })

  df_exp_queried <- reactive({
    req(input$query_gene)
    df_exp %>%
      dplyr::select(depmap_id, input$query_gene) %>%
      dplyr::rename("Expression" = input$query_gene)
  })

  df_copynumber_queried <- reactive({
    req(input$query_gene)
    df_copynumber %>%
      dplyr::select(depmap_id, input$query_gene) %>%
      dplyr::rename("Copy number" = input$query_gene)
  })

  df_mut_queried <- reactive({
    req(input$query_gene)
    df_mut %>%
      dplyr::select(depmap_id, input$query_gene) %>%
      dplyr::rename("Mutation" = input$query_gene)
  })

  df_joined <- reactive({
    df <- df_metadata %>%
      dplyr::left_join(df_drug_sen_prism_queried()) %>%
      dplyr::left_join(df_drug_sen_prism_auc_queried()) %>%
      dplyr::left_join(df_drug_sen_gdsc_auc_queried())

    if (input$query_gene_feature == "Expression") {
      df <- df %>% dplyr::left_join(df_exp_queried())
    } else if (input$query_gene_feature == "Copy number") {
      df <- df %>% dplyr::left_join(df_copynumber_queried())
    } else if (input$query_gene_feature == "Mutation") {
      df <- df %>% dplyr::left_join(df_mut_queried())
    }

    df
  })

  df_joined_pivot <- reactive({
    df_joined() %>%
      dplyr::select(-depmap_id, -cell_line_name) %>%
      tidyr::pivot_longer(contains(c("PRISM", "PRISM AUC", "GDSC AUC")),
                          names_to = "Sensitivity_type",
                          values_to = "Sensitivity_value")
  })

  correlation <- reactive({
    nrow_prism <-
      df_joined() %>% dplyr::filter(!is.na(PRISM)) %>% nrow()
    nrow_prism_auc <-
      df_joined() %>% dplyr::filter(!is.na(`PRISM AUC`)) %>% nrow()
    nrow_gdsc_auc <-
      df_joined() %>% dplyr::filter(!is.na(`GDSC AUC`)) %>% nrow()

    df <- df_joined() %>% dplyr::select(full_disease_name) %>% distinct()

    if (nrow_prism > 0) {
      df_prism <-
        df_joined() %>%
        dplyr::select(full_disease_name,
                      PRISM,
                      input$query_gene_feature) %>%
        dplyr::rename("Value" = input$query_gene_feature) %>%
        dplyr::group_by(full_disease_name) %>%
        dplyr::mutate(n = n()) %>%
        dplyr::filter(n > 2) %>%
        dplyr::summarise(PRISM_cor = cor(`PRISM`,
                                         `Value`,
                                         method = "pearson",
                                         use = "pairwise.complete.obs")) %>%
        dplyr::mutate(PRISM_cor = round(PRISM_cor, 3))

      df <- dplyr::left_join(df, df_prism, by = "full_disease_name")
    }

    if (nrow_prism_auc > 0) {
      df_prism_auc <-
        df_joined() %>%
        dplyr::select(full_disease_name,
                      `PRISM AUC`,
                      input$query_gene_feature) %>%
        dplyr::rename("Value" = input$query_gene_feature) %>%
        dplyr::group_by(full_disease_name) %>%
        dplyr::mutate(n = n()) %>%
        dplyr::filter(n > 2) %>%
        dplyr::summarise(
          PRISM_AUC_cor = cor(`PRISM AUC`,
                              `Value`,
                              method = "pearson",
                              use = "pairwise.complete.obs")
        ) %>%
        dplyr::mutate(PRISM_AUC_cor = round(PRISM_AUC_cor, 3))
      
      df <-
        dplyr::left_join(df, df_prism_auc, by = "full_disease_name")
    }

    if (nrow_gdsc_auc > 0) {
      df_gdsc_auc <-
        df_joined() %>%
        dplyr::select(full_disease_name,
                      `GDSC AUC`,
                      input$query_gene_feature) %>%
        dplyr::rename("Value" = input$query_gene_feature) %>%
        dplyr::group_by(full_disease_name) %>%
        dplyr::mutate(n = n()) %>%
        dplyr::filter(n > 2) %>%
        dplyr::summarise(
          GDSC_AUC_cor = cor(`GDSC AUC`,
                              `Value`,
                              method = "pearson",
                              use = "pairwise.complete.obs")
        ) %>%
        dplyr::mutate(GDSC_AUC_cor = round(GDSC_AUC_cor, 3))

      df <-
        dplyr::left_join(df, df_gdsc_auc, by = "full_disease_name")
    }

    df
  })

  correlation_to_show <- reactive({
    correlation() %>%
      dplyr::filter((PRISM_cor < 0 & PRISM_AUC_cor < 0) |
                      (PRISM_cor < 0 & GDSC_AUC_cor < 0) |
                      (PRISM_AUC_cor < 0 & GDSC_AUC_cor < 0)
      ) %>%
      dplyr::filter(abs(PRISM_cor) != 1,
                    abs(PRISM_AUC_cor) != 1,
                    abs(GDSC_AUC_cor) != 1) %>%
      dplyr::arrange(PRISM_cor)
  })

  output$tbl_cor <- DT::renderDT(server = FALSE, {
    DT::datatable(
      correlation_to_show() %>%
        dplyr::rename("Primary disease - subtype" = "full_disease_name"),
      extensions = "Buttons",
      rownames = FALSE,
      selection = "single",
      options = list(
        dom = "Bfrtlip",
        lengthMenu = c(5, 25, 50),
        pageLength = 25,
        buttons = list(
          list(
            extend = "excel",
            text = "Export current results",
            filename = paste0(input$query_cpd, "_current_data"),
            exportOptions = list(modifier = list(page = "current"))
          ),
          list(
            extend = "excel",
            text = "Export full results",
            filename = paste0(input$query_cpd, "_full_data"),
            exportOptions = list(modifier = list(page = "all"))
          )
        )
      )
    )
  })

  data_scatter <- reactive({
    req(input$query_cpd)
    req(input$query_gene)

    if (length(input$tbl_cor_cell_clicked$row) > 0) {
      disease_clicked <-
        correlation_to_show()[input$tbl_cor_cell_clicked$row, ] %>%
        dplyr::pull(full_disease_name)
    } else {
      disease_clicked <-
        correlation_to_show()[1, ] %>%
        dplyr::pull(full_disease_name)
    }

    df <- df_joined() %>%
      tidyr::pivot_longer(contains(c("PRISM", "PRISM AUC", "GDSC AUC")),
                          names_to = "Sensitivity_type",
                          values_to = "Sensitivity_value") %>%
      dplyr::filter(full_disease_name == disease_clicked) %>%
      dplyr::select(-depmap_id, -full_disease_name) %>%
      dplyr::rename("Feature_value" = input$query_gene_feature) %>%
      dplyr::filter(!is.na(Sensitivity_value), !is.na(Feature_value)) %>%
      dplyr::group_by(Sensitivity_type) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::filter(n > 2)

    df
  })

  output$tbl_scatter <- DT::renderDT({
    DT::datatable(
      data_scatter() %>%
        dplyr::mutate(Sensitivity_value = round(Sensitivity_value, 3)) %>%
        dplyr::mutate(Feature_value = round(Feature_value, 3)) %>%
        dplyr::select(
          cell_line_name,
          Sensitivity_type,
          Sensitivity_value,
          Feature_value
        ) %>%
        dplyr::arrange(Sensitivity_value) %>%
        dplyr::arrange(desc(Feature_value)),
      rownames = FALSE,
      selection = "none",
      options = list(dom = "tlip")
    )
  })

  fig_scatter <- reactive({
    if (length(input$tbl_cor_cell_clicked$row) > 0) {
      disease_clicked <-
        correlation_to_show()[input$tbl_cor_cell_clicked$row, ] %>%
        dplyr::pull(full_disease_name)
    } else {
      disease_clicked <-
        correlation_to_show()[1, ] %>%
        dplyr::pull(full_disease_name)
    }

    df <- data_scatter()
    df$Sensitivity_type <-
      factor(df$Sensitivity_type,
             levels = c("PRISM", "PRISM AUC", "GDSC AUC"))

    ggscatter(
      df,
      x = "Sensitivity_value",
      y = "Feature_value",
      title = disease_clicked,
      add = "reg.line",
      color = "Sensitivity_type",
      facet.by = "Sensitivity_type",
      scales = "free_x",
      palette = "lancet",
      label = "cell_line_name",
      repel = TRUE
    ) +
      ylim(c(min(df$Feature_value), max(df$Feature_value) * 1.1)) +
      xlab(paste("Sensitivity of", input$query_cpd)) +
      ylab(paste(input$query_gene_feature, "of", input$query_gene)) +
      stat_cor(
        aes(color = Sensitivity_type),
        method = "pearson",
        label.y = max(df$Feature_value) * 1.1
      ) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  })

  output$fig_scatter_out <- renderPlot({
    fig_scatter()
  })

  fig_venn <- reactive({
    df <- data_scatter()
    df$Sensitivity_type <-
      factor(df$Sensitivity_type,
             levels = c("PRISM", "PRISM AUC", "GDSC AUC"))

    sets <- list(
      `PRISM` = df %>%
        dplyr::filter(Sensitivity_type == "PRISM") %>%
        dplyr::pull(cell_line_name),
      `PRISM AUC` = df %>%
        dplyr::filter(Sensitivity_type == "PRISM AUC") %>%
        dplyr::pull(cell_line_name),
      `GDSC AUC` = df %>%
        dplyr::filter(Sensitivity_type == "GDSC AUC") %>%
        dplyr::pull(cell_line_name)
    )

    ggvenn(
      sets,
      c("PRISM", "PRISM AUC", "GDSC AUC"),
      fill_color = c("#005096", "#EF0000", "#4BBD49"),
      fill_alpha = 0.45,
      digits = 0
    ) + ggtitle("Comparison of cell-lines used") +
      theme(plot.title = element_text(hjust = 0.5))
  })

  output$fig_venn_out <- renderPlot({
    fig_venn()
  })

  output$download_fig <- downloadHandler(
    filename = function() {
      paste0(input$query_cpd,
             " vs ",
             input$query_gene,
             " ",
             input$query_gene_feature,
             ".png")
    },
    content = function(file) {
      png(
        file = file,
        width = 10,
        height = 5,
        units = "in",
        res = 300
      )
      plot(fig_scatter())
      dev.off()
    }
  )
  
  output$download_tab <- downloadHandler(
    filename = function() {
      paste0(input$query_cpd,
             " vs ",
             input$query_gene,
             " ",
             input$query_gene_feature,
             ".csv")
    },
    content = function(file) {
      write.csv(
        data_scatter() %>%
          dplyr::mutate(Sensitivity_value = round(Sensitivity_value, 3)) %>%
          dplyr::mutate(Feature_value = round(Feature_value, 3)) %>%
          dplyr::select(
            cell_line_name,
            Sensitivity_type,
            Sensitivity_value,
            Feature_value
          ),
        file,
        row.names = FALSE
      )
    }
  )
}

shinyApp(ui, server)
