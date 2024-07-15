pacman::p_load(
  Cairo,
  shiny,
  shinydashboard,
  tidyverse,
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
      selectizeInput("compound_input",
                     label = "Compound",
                     choices = NULL)
    ),
    column(4,
           selectizeInput(
             "gene_input",
             label = "Gene",
             choices = NULL
           )),
    column(
      4,
      selectizeInput(
        "gene_feature_input",
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
      DT::DTOutput("correlation_table", height = "100%") %>% withSpinner()
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
        plotOutput("scatter_plot_output") %>% withSpinner(),
        downloadLink("download_scatter_plot", "Download figure"),
        plotOutput("venn_plot_output") %>% withSpinner()
      ),
      tabPanel(
        "Data",
        NULL,
        DT::DTOutput("scatter_data_table", height = "100%") %>% withSpinner(),
        downloadButton("download_data_table", "Download table")
      )
    )
  )))
)

server <- function(input, output, session) {
  # Update the compound selectize input with choices from df_drug_sen_prism_auc
  updateSelectizeInput(
    session,
    "compound_input",
    choices = levels(factor(colnames(
      df_drug_sen_prism_auc[, -1]
    ))),
    selected = "CX-5461",
    server = TRUE
  )

  # Update the gene selectize input with choices from df_exp
  updateSelectizeInput(
    session,
    "gene_input",
    choices = levels(factor(colnames(df_exp[, -1]))),
    selected = "MYC",
    server = TRUE
  )

  # Reactive function to query df_drug_sen_prism based on selected compound
  queried_prism <- reactive({
    req(input$compound_input)
    df <- df_drug_sen_prism %>%
      select(depmap_id, contains(input$compound_input))

    if (ncol(df) == 2) {
      df %>% rename("PRISM" = input$compound_input)
    } else {
      df$`PRISM` <- NA
      df
    }
  })

  # Reactive function to query df_drug_sen_prism_auc based on selected compound
  queried_prism_auc <- reactive({
    req(input$compound_input)
    df <- df_drug_sen_prism_auc %>%
      select(depmap_id, contains(input$compound_input))

    if (ncol(df) == 2) {
      df %>% rename("PRISM AUC" = input$compound_input)
    } else {
      df$`PRISM AUC` <- NA
      df
    }
  })

  # Reactive function to query df_drug_sen_gdsc_auc based on selected compound
  queried_gdsc_auc <- reactive({
    req(input$compound_input)
    df <- df_drug_sen_gdsc_auc %>%
      select(depmap_id, contains(input$compound_input))

    if (ncol(df) == 2) {
      df %>% rename("GDSC AUC" = input$compound_input)
    } else {
      df$`GDSC AUC` <- NA
      df
    }
  })

  # Reactive function to query df_exp based on selected gene
  queried_exp <- reactive({
    req(input$gene_input)
    df_exp %>%
      select(depmap_id, input$gene_input) %>%
      rename("Expression" = input$gene_input)
  })
  
  # Reactive function to query df_copynumber based on selected gene
  queried_copy_number <- reactive({
    req(input$gene_input)
    df_copynumber %>%
      select(depmap_id, input$gene_input) %>%
      rename("Copy number" = input$gene_input)
  })
  
  # Reactive function to query df_mut based on selected gene
  queried_mutation <- reactive({
    req(input$gene_input)
    df_mut %>%
      select(depmap_id, input$gene_input) %>%
      rename("Mutation" = input$gene_input)
  })
  
  # Reactive function to join multiple dataframes based on selected gene feature
  joined_data <- reactive({
    df <- df_metadata %>%
      left_join(queried_prism()) %>%
      left_join(queried_prism_auc()) %>%
      left_join(queried_gdsc_auc())
  
    if (input$gene_feature_input == "Expression") {
      df <- df %>% left_join(queried_exp())
    } else if (input$gene_feature_input == "Copy number") {
      df <- df %>% left_join(queried_copy_number())
    } else if (input$gene_feature_input == "Mutation") {
      df <- df %>% left_join(queried_mutation())
    }
  
    df
  })
  
  # Reactive function to pivot the joined dataframe for sensitivity analysis
  pivoted_data <- reactive({
    joined_data() %>%
      select(-depmap_id, -cell_line_name) %>%
      pivot_longer(contains(c("PRISM", "PRISM AUC", "GDSC AUC")),
                          names_to = "Sensitivity_type",
                          values_to = "Sensitivity_value")
  })
  
  # Reactive function to calculate correlation between sensitivity and gene feature
  correlation_calculated <- reactive({
    nrow_prism <-
      joined_data() %>% filter(!is.na(PRISM)) %>% nrow()
    nrow_prism_auc <-
      joined_data() %>% filter(!is.na(`PRISM AUC`)) %>% nrow()
    nrow_gdsc_auc <-
      joined_data() %>% filter(!is.na(`GDSC AUC`)) %>% nrow()
  
    df <- joined_data() %>% select(full_disease_name) %>% distinct()
  
    if (nrow_prism > 0) {
      df_prism <-
        joined_data() %>%
        select(full_disease_name,
                      PRISM,
                      input$gene_feature_input) %>%
        rename("Value" = input$gene_feature_input) %>%
        group_by(full_disease_name) %>%
        mutate(n = n()) %>%
        filter(n > 2) %>%
        summarise(PRISM_cor = cor(`PRISM`,
                                         `Value`,
                                         method = "pearson",
                                         use = "pairwise.complete.obs")) %>%
        mutate(PRISM_cor = round(PRISM_cor, 3))
  
      df <- left_join(df, df_prism, by = "full_disease_name")
    }
  
    if (nrow_prism_auc > 0) {
      df_prism_auc <-
        joined_data() %>%
        select(full_disease_name,
                      `PRISM AUC`,
                      input$gene_feature_input) %>%
        rename("Value" = input$gene_feature_input) %>%
        group_by(full_disease_name) %>%
        mutate(n = n()) %>%
        filter(n > 2) %>%
        summarise(
          PRISM_AUC_cor = cor(`PRISM AUC`,
                              `Value`,
                              method = "pearson",
                              use = "pairwise.complete.obs")
        ) %>%
        mutate(PRISM_AUC_cor = round(PRISM_AUC_cor, 3))
  
      df <-
        left_join(df, df_prism_auc, by = "full_disease_name")
    }
  
    if (nrow_gdsc_auc > 0) {
      df_gdsc_auc <-
        joined_data() %>%
        select(full_disease_name,
                      `GDSC AUC`,
                      input$gene_feature_input) %>%
        rename("Value" = input$gene_feature_input) %>%
        group_by(full_disease_name) %>%
        mutate(n = n()) %>%
        filter(n > 2) %>%
        summarise(
          GDSC_AUC_cor = cor(`GDSC AUC`,
                              `Value`,
                              method = "pearson",
                              use = "pairwise.complete.obs")
        ) %>%
        mutate(GDSC_AUC_cor = round(GDSC_AUC_cor, 3))
  
      df <-
        left_join(df, df_gdsc_auc, by = "full_disease_name")
    }
  
    df
  })

  correlation_filtered <- reactive({
    correlation_calculated() %>%
      filter((PRISM_cor < 0 & PRISM_AUC_cor < 0) |
                      (PRISM_cor < 0 & GDSC_AUC_cor < 0) |
                      (PRISM_AUC_cor < 0 & GDSC_AUC_cor < 0)
      ) %>%
      filter(abs(PRISM_cor) != 1,
                    abs(PRISM_AUC_cor) != 1,
                    abs(GDSC_AUC_cor) != 1) %>%
      arrange(PRISM_cor)
  })

  # Render the correlation table
  output$correlation_table <- DT::renderDT(server = FALSE, {
    DT::datatable(
      correlation_filtered() %>%
        rename("Primary disease - subtype" = "full_disease_name"),
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
            filename = paste0(input$compound_input, "_current_data"),
            exportOptions = list(modifier = list(page = "current"))
          ),
          list(
            extend = "excel",
            text = "Export full results",
            filename = paste0(input$compound_input, "_full_data"),
            exportOptions = list(modifier = list(page = "all"))
          )
        )
      )
    )
  })

  scatter_data <- reactive({
    req(input$compound_input)
    req(input$gene_input)

    # Determine the disease clicked in the correlation table
    if (length(input$correlation_table_cell_clicked$row) > 0) {
      disease_clicked <-
        correlation_filtered()[input$correlation_table_cell_clicked$row, ] %>%
        pull(full_disease_name)
    } else {
      disease_clicked <-
        correlation_filtered()[1, ] %>%
        pull(full_disease_name)
    }

    # Prepare the data for scatter plot
    df <- joined_data() %>%
      pivot_longer(contains(c("PRISM", "PRISM AUC", "GDSC AUC")),
                          names_to = "Sensitivity_type",
                          values_to = "Sensitivity_value") %>%
      filter(full_disease_name == disease_clicked) %>%
      select(-depmap_id, -full_disease_name) %>%
      rename("Feature_value" = input$gene_feature_input) %>%
      filter(!is.na(Sensitivity_value), !is.na(Feature_value)) %>%
      group_by(Sensitivity_type) %>%
      mutate(n = n()) %>%
      filter(n > 2)

    df
  })

  # Render the scatter plot table
  output$scatter_data_table <- DT::renderDT({
    DT::datatable(
      scatter_data() %>%
        mutate(Sensitivity_value = round(Sensitivity_value, 3)) %>%
        mutate(Feature_value = round(Feature_value, 3)) %>%
        select(
          cell_line_name,
          Sensitivity_type,
          Sensitivity_value,
          Feature_value
        ) %>%
        arrange(Sensitivity_value) %>%
        arrange(desc(Feature_value)),
      rownames = FALSE,
      selection = "none",
      options = list(dom = "tlip")
    )
  })

  scatter_plot <- reactive({
    # Determine the disease clicked in the correlation table
    if (length(input$correlation_table_cell_clicked$row) > 0) {
      disease_clicked <-
        correlation_filtered()[input$correlation_table_cell_clicked$row, ] %>%
        pull(full_disease_name)
    } else {
      disease_clicked <-
        correlation_filtered()[1, ] %>%
        pull(full_disease_name)
    }

    df <- scatter_data()
    df$Sensitivity_type <-
      factor(df$Sensitivity_type,
             levels = c("PRISM", "PRISM AUC", "GDSC AUC"))

    # Generate the scatter plot
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
      xlab(paste("Sensitivity of", input$compound_input)) +
      ylab(paste(input$gene_feature_input, "of", input$gene_input)) +
      stat_cor(
        aes(color = Sensitivity_type),
        method = "pearson",
        label.y = max(df$Feature_value) * 1.1
      ) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  })

  # Render the scatter plot
  output$scatter_plot_output <- renderPlot({
    scatter_plot()
  })

  venn_plot <- reactive({
    df <- scatter_data()
    df$Sensitivity_type <-
      factor(df$Sensitivity_type,
             levels = c("PRISM", "PRISM AUC", "GDSC AUC"))

    # Create sets for Venn diagram
    sets <- list(
      `PRISM` = df %>%
        filter(Sensitivity_type == "PRISM") %>%
        pull(cell_line_name),
      `PRISM AUC` = df %>%
        filter(Sensitivity_type == "PRISM AUC") %>%
        pull(cell_line_name),
      `GDSC AUC` = df %>%
        filter(Sensitivity_type == "GDSC AUC") %>%
        pull(cell_line_name)
    )

    # Generate the Venn diagram
    ggvenn(
      sets,
      c("PRISM", "PRISM AUC", "GDSC AUC"),
      fill_color = c("#005096", "#EF0000", "#4BBD49"),
      fill_alpha = 0.45,
      digits = 0
    ) + ggtitle("Comparison of cell-lines used") +
      theme(plot.title = element_text(hjust = 0.5))
  })

  # Render the Venn diagram
  output$venn_plot_output <- renderPlot({
    venn_plot()
  })

  # Download handler for scatter plot
  output$download_scatter_plot <- downloadHandler(
    filename = function() {
      paste0(input$compound_input,
             " vs ",
             input$gene_input,
             " ",
             input$gene_feature_input,
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
      plot(scatter_plot())
      dev.off()
    }
  )

  # Download handler for data table
  output$download_data_table <- downloadHandler(
    filename = function() {
      paste0(input$compound_input,
             " vs ",
             input$gene_input,
             " ",
             input$gene_feature_input,
             ".csv")
    },
    content = function(file) {
      write.csv(
        scatter_data() %>%
          mutate(Sensitivity_value = round(Sensitivity_value, 3)) %>%
          mutate(Feature_value = round(Feature_value, 3)) %>%
          select(
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

# Create and run the Shiny app
shinyApp(ui, server)
