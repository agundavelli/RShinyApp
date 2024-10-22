library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(viridis)
library(colourpicker)
library(ggbeeswarm)

options(shiny.maxRequestSize = 100 * 1024^2) # Set max upload size to 100 MB

# Function to convert metadata
metadata_converted <- function(metadata) {
  # Extract and convert data
  diagnosis <- substr(metadata[, 2], 12, nchar(metadata[, 2]))
  pmi_numeric <- as.numeric(gsub("[^0-9.]", "", substr(metadata[, 3], 6, nchar(metadata[, 3]))))
  age_of_death_numeric <- as.numeric(gsub("\\D", "", substr(metadata[, 4], 15, nchar(metadata[, 4]))))
  rin_numeric <- as.numeric(gsub("[^0-9.]", "", substr(metadata[, 5], 6, nchar(metadata[, 5]))))
  mrna_seq_reads_numeric <- as.numeric(gsub("\\D", "", substr(metadata[, 6], 17, nchar(metadata[, 6]))))
  
  converted_data <- data.frame(
    "Sample" = metadata[, 1],
    "Diagnosis" = diagnosis,
    "PMI" = pmi_numeric,
    "Age of Death" = age_of_death_numeric,
    "Rin" = rin_numeric,
    "mRNA-Seq Reads" = mrna_seq_reads_numeric
  )
  converted_data <- na.omit(converted_data)
  
  return(converted_data)
}

# File Input Validation Function
validate_file <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("File not found.")
  }
  
  # Check if the file extension is either CSV or TSV
  file_extension <- tools::file_ext(file_path)
  if (!(file_extension %in% c("csv", "tsv"))) {
    stop("Invalid file format. Only CSV or TSV files are allowed.")
  }
  
  # Additional checks for well-formatted file can be added if needed
}

# Define UI for application
ui <- dashboardPage(
  dashboardHeader(title = "BF591 Final Project"),
  dashboardSidebar(
    fileInput("metadata_file", "Choose Metadata File (CSV)", accept = c(".csv")),
    actionButton("submit_btn", "Submit")
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML(".subheader {font-size: 18px; font-weight: bold;}"))
    ),
    tags$div(class = "subheader", "Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"),
    tabsetPanel(
      tabPanel("Samples", 
               tabsetPanel(
                 tabPanel("Summary", DTOutput("summary_table")),  # Use DTOutput for the summary table
                 tabPanel("Table", DTOutput("metadata_table")),
                 tabPanel("Density Plots", 
                          sidebarPanel(
                            selectInput("density_plot_column", "Choose Numeric Column to Plot", 
                                        choices = c("PMI", "Age of Death", "Rin", "mRNA-Seq Reads")),
                            actionButton("submit_plot_btn", "Submit")  # Changed action button ID
                          ),
                          mainPanel(
                            plotlyOutput("density_plot")
                          )
                 )
               )
      ),
      tabPanel("Counts",
               dashboardSidebar(
                 width = 300,
                 fileInput("file", "Choose Counts File (CSV)", accept = c(".csv")),
                 sliderInput("variance_threshold", "Variance Threshold",
                             min = 0, max = 100, value = 0,
                             step = 1, post = "%"),
                 sliderInput("nonzero_samples_threshold", "Non-Zero Samples Threshold",
                             min = 0, max = 69, value = 0,
                             step = 1),
                 actionButton("submit_button", "Submit")
               ),
               tabsetPanel(
                 tabPanel("Filter Summary", DTOutput("filter_summary")),
                 tabPanel("Diagnostic Plots",
                          plotOutput("median_vs_variance"),
                          plotOutput("median_vs_zeros")),
                 tabPanel("Clustered Heatmap", plotOutput("clustered_heatmap", width = "100%", height = "600px")),
                 tabPanel("PCA Scatter Plot",
                          sidebarPanel(
                            # Add input fields to allow selection of PCs for x and y axes
                            selectInput("x_pc", "Choose X-axis Principal Component:", 
                                        choices = paste0("PC", 1:10), selected = "PC1"),
                            selectInput("y_pc", "Choose Y-axis Principal Component:", 
                                        choices = paste0("PC", 1:10), selected = "PC2")
                          ),
                          mainPanel(
                            # Plot output for PCA scatter plot
                            plotOutput("pca_scatter_plot")
                          )
                 )
               )
               
      ),
      tabPanel("Differential Expression",
               dashboardSidebar(
                 fileInput("file_expr", "Choose Differential Expression Data (CSV)", accept = c(".csv")),
                 actionButton("submit_btn_expr", "Submit")
               ),
               tabsetPanel(
                 tabPanel("Differential Expression Results",
                          # Add a container for the table with horizontal scrolling enabled
                          div(style = 'overflow-x: scroll', DTOutput("expressionTable"))
                 ),
                 tabPanel("Plot Results",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("xAxis", "Choose the column for the x-axis: ", choices = NULL),
                              selectInput("yAxis", "Choose the column for the y-axis", choices = NULL),
                              colourInput("color1", "Base point color", value = "#22577A"),
                              colourInput("color2", "Highlight point color", value = "#FFCF56"),
                              sliderInput("magnitudeSlider", "Select the magnitude of the p-adjusted coloring:",
                                          min = -35, max = 0, value = -15),
                              actionButton("submit_plot_btn", "Submit"),
                              # Set a fixed height for the sidebar panel
                              style = "height: 600px; overflow-y: auto;"
                            ),
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Volcano Plot", 
                                         # Set the height of the volcano plot to match the sidebar
                                         plotOutput("volcano", height = "557px")
                                ),
                                tabPanel("P-adjusted Filtered Table", 
                                         # Add horizontal scrolling to the table, set height to match
                                         div(style = 'overflow-x: scroll; height: 580px;', 
                                             DTOutput("table")
                                         )
                                )
                              )
                            )
                          )
                 )
               )
      ),
      # Gene Expression Visualization tab
      tabPanel("Gene Expression Visualization",
               dashboardSidebar(
                 fileInput("counts_file", "Choose Gene Counts CSV File", accept = c(".csv")),
                 fileInput("info_file", "Choose Sample Information CSV File", accept = c(".csv")),
                 selectInput("categorical_var", "Select Categorical Variable", "Diagnosis", choices = "Diagnosis"),
                 selectInput("gene_var", "Select Gene", ""),
                 selectInput("plot_type", "Select Plot Type",
                             choices = c("Bar Plot", "Boxplot", "Violin Plot", "Beeswarm Plot")),
                 actionButton("submit_button", "Submit")
               ),
               plotOutput("gene_plot")
      )
      
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive function to read metadata when a file is uploaded
  metadata <- reactive({
    req(input$metadata_file)
    read.csv(input$metadata_file$datapath)
  })
  
  # Reactive function to perform summary analysis
  summary_df <- reactive({
    req(metadata())
    
    metadata_converted_data <- metadata_converted(metadata())
    
    # Create summary dataframe
    summary_df <- data.frame(
      "Column Name" = names(metadata_converted_data)[-1],
      "Type" = sapply(metadata_converted_data[-1], function(x) {
        if (is.numeric(x) && all(x %% 1 == 0)) {
          return("integer")
        } 
        else if (is.numeric(x)) {
          return("double")
        } 
        else if (is.factor(x)) {
          return(paste(levels(x), collapse = "/"))
        } 
        else {
          return(class(x))
        }
      }),
      "Mean (sd) or Distinct Values" = sapply(metadata_converted_data[-1], function(x) {
        if (is.numeric(x)) {
          return(paste0(round(mean(x), 2), " (+/- ", round(sd(x), 2), ")"))
        } else {
          return(paste0(paste(levels(factor(x)), collapse = ", ")))
        }
      })
    )
    rownames(summary_df) <- NULL  # Remove row names
    
    return(summary_df)
  })
  
  # Render summary table
  output$summary_table <- renderDT({
    req(summary_df())
    if (!input$submit_btn) {
      return(NULL)
    }
    datatable(summary_df(), options = list(scrollX = TRUE))
  })
  
  # Render metadata table
  output$metadata_table <- renderDT({
    req(metadata())
    datatable(metadata_converted(metadata()), options = list(scrollX = TRUE))
  })
  
  observe({
    data_df <- metadata_converted(metadata())
    numeric_columns <- sapply(data_df, is.numeric)
    updateSelectInput(session, "density_plot_column", choices = names(data_df)[numeric_columns])
  })
  
  # Render density plot
  output$density_plot <- renderPlotly({
    req(metadata(), input$submit_plot_btn)
    
    # Allow the user to choose columns
    column_to_plot <- input$density_plot_column
    
    # Check if the selected column is valid
    if (!(column_to_plot %in% names(metadata_converted(metadata())))) {
      return(NULL)
    }
    
    # Create a density plot
    plot_data <- ggplot(metadata_converted(metadata()), aes(x = .data[[column_to_plot]])) +
      geom_density(fill = "skyblue", color = "black", alpha = 0.7) +  # Adjusted fill color and added transparency
      labs(title = paste("Density Plot of", column_to_plot),
           x = column_to_plot,
           y = "Density") +
      theme_minimal() +  # Used a minimal theme for better clarity
      theme(plot.title = element_text(hjust = 0.5))  # Centered the plot title
    
    ggplotly(plot_data)
  })
  
  # Read the counts matrix
  counts_data <- reactive({
    req(input$file)
    validate_file(input$file$datapath)
    read.csv(input$file$datapath, stringsAsFactors = FALSE, row.names = 1)  # Assuming gene names are in the first column
  })
  
  # Filter genes based on variance and non-zero samples thresholds
  filtered_counts_data <- reactive({
    req(counts_data())
    
    variance_threshold <- quantile(apply(counts_data(), 1, var), input$variance_threshold / 100)
    non_zero_samples_threshold <- input$nonzero_samples_threshold
    
    filtered_genes <- rownames(counts_data())[apply(counts_data(), 1, function(x) var(x) >= variance_threshold & sum(x > 0) >= non_zero_samples_threshold)]
    
    counts_data()[filtered_genes, , drop = FALSE]
  })
  
  # Display filter summary table
  output$filter_summary <- renderDT({
    total_genes <- nrow(counts_data())
    total_samples <- ncol(counts_data())
    filtered_genes <- nrow(filtered_counts_data())
    remaining_genes <- total_genes - filtered_genes
    
    percentage_passed <- paste0(round(filtered_genes / total_genes * 100, 1), "%")
    percentage_not_passed <- paste0(round(remaining_genes / total_genes * 100, 1), "%")
    
    filter_summary_data <- data.frame(
      Metric = c("Total Samples", "Total Genes", "Number of Genes Passing Current Filter", "% of Genes Passing Current Filter", 
                 "Number of Genes Not Passing Current Filter", "% of Genes Not Passing Current Filter"),
      Value = c(total_samples, total_genes, filtered_genes, percentage_passed,remaining_genes, percentage_not_passed)
    )
    
    datatable(filter_summary_data)
  })
  
  # Diagnostic scatter plots
  output$median_vs_variance <- renderPlot({
    # Plot median count vs variance
    plot(log(apply(counts_data(), 1, median) + 1), log(apply(counts_data(), 1, var) + 1),
         col = ifelse(rownames(counts_data()) %in% rownames(filtered_counts_data()), "blue", "red"),
         pch = 16, cex = 0.7,
         xlab = "Log(Median Count)", ylab = "Log(Variance)",
         main = "Log(Median Count) vs Log(Variance)")
    
    legend("topright", legend = c("Passed Filter", "Did Not Pass Filter"), col = c("blue", "red"), pch = 16)
  })
  
  output$median_vs_zeros <- renderPlot({
    # Plot median count vs number of zeros
    plot(log(apply(counts_data(), 1, median) + 1), log(apply(counts_data() == 0, 1, sum) + 1),
         col = ifelse(rownames(counts_data()) %in% rownames(filtered_counts_data()), "blue", "red"),
         pch = 16, cex = 0.7,
         xlab = "Log(Median Count)", ylab = "Log(Number of Zeros)",
         main = "Log(Median Count) vs Log(Number of Zeros)")
    
    legend("topright", legend = c("Passed Filter", "Did Not Pass Filter"), col = c("blue", "red"), pch = 16)
  })
  
  # Clustered heatmap
  output$clustered_heatmap <- renderPlot({
    # Plot clustered heatmap of counts matrix
    heatmap(log1p(as.matrix(filtered_counts_data())), scale = "row", col = viridis::viridis(256), main = "Clustered Heatmap",
            xlab = "Samples", ylab = "Genes")
  })
  
  # PCA scatter plot
  output$pca_scatter_plot <- renderPlot({
    # Perform PCA on the filtered counts matrix
    pca_result <- prcomp(t(log1p(filtered_counts_data())))
    
    # Get the selected PCs from the user input
    x_pc <- as.numeric(sub("PC", "", input$x_pc))  # Extract the number from "PCx"
    y_pc <- as.numeric(sub("PC", "", input$y_pc)) 
    
    # Plot PCA scatter plot with the selected PCs
    plot(pca_result$x[, x_pc], pca_result$x[, y_pc], col = "blue", pch = 16, cex = 0.7,
         xlab = paste("PC", x_pc), ylab = paste("PC", y_pc),
         main = paste("PCA Scatter Plot: PC", x_pc, "vs PC", y_pc))
  })
  
  # Reactive function to read data when a file is uploaded in the Differential Expression tab
  data_expr <- reactive({
    req(input$file_expr)
    read.csv(input$file_expr$datapath)
  })
  
  
  # Render the expression table
  output$expressionTable <- renderDT({
    if (!input$submit_btn_expr) {
      return(NULL)
    }
    datatable(
      data_expr(),  # Fix: Change from data() to data_expr()
      options = list(
        pageLength = 10,
        searchHighlight = TRUE
      ),
      filter = list(position = 'top', clear = FALSE),
      class = 'cell-border stripe',
      rownames = FALSE
    )
  })
  
  observe({
    data_df <- data_expr()
    updateSelectInput(session, "xAxis", choices = names(data_df))
    updateSelectInput(session, "yAxis", choices = names(data_df))
  })
  
  
  # Reactive function to load data
  load_data <- reactive({
    req(input$file$datapath)
    read.csv(input$file$datapath)
  })
  
  # Observe block to update choices for x and y axes
  observe({
    data_df <- data_expr()
    # Remove columns 1, 2, 4, and 5 from the choices
    choices <- setdiff(names(data_df), c("Gene", "symbol", "HD.mean", "Control.mean"))
    updateSelectInput(session, "xAxis", choices = choices, selected = "log2FoldChange")
    updateSelectInput(session, "yAxis", choices = choices, selected = "padj")
  })
  
  # Volcano plot function
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    plot <- ggplot(dataf, aes_string(x = x_name, y = sprintf("-log10(%s) + 1e-10", y_name), color = sprintf("factor(padj < 10^%d)", slider))) +
      geom_point(aes(), shape = 19) +
      scale_color_manual(values = c("FALSE" = color1, "TRUE" = color2, "NA" = "gray"), name = sprintf("padj < 1 × 10^%d", slider)) +
      labs(x = x_name, y = sprintf("-log10(%s)", y_name)) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    return(plot)
  }
  
  # Draw and filter table function
  draw_table_expr <- function(dataf, slider) {
    # remove rows with NAs
    dataf <- na.omit(dataf)
    # filter the data frame based on the slider magnitude
    filtered_data <- dataf[dataf$padj < 10^slider, ]
    # format p-value columns to display more digits
    filtered_data$pvalue <- formatC(filtered_data$pvalue, format = "e", digits = 3)
    filtered_data$padj <- formatC(filtered_data$padj, format = "e", digits = 3)
    
    return(filtered_data)
  }
  
  # Reactive function to track if the "Submit" button in the "Plot Results" tab is clicked
  plot_submit_clicked_expr <- reactiveVal(FALSE)
  
  # Render the volcano plot
  output$volcano <- renderPlot({
    req(input$submit_btn_expr)
    
    data <- data_expr()
    
    # Call the volcano_plot function and store the plot
    plot <- volcano_plot(data, input$xAxis, input$yAxis, input$magnitudeSlider, input$color1, input$color2)
    
    # Return the plot
    return(plot)
  })
  
  # Render the filtered table for differential expression results
  output$table <- renderDT({
    req(input$submit_btn_expr)
    
    data <- data_expr()
    filtered_table <- draw_table_expr(data, input$magnitudeSlider)
    
    # Use datatable function with options for sorting
    datatable(
      filtered_table,
      options = list(order = list(list(4, 'asc'))),  # Sort by the fourth column (padj) in ascending order
      class = 'cell-border stripe',
      rownames = FALSE
    )
  })
  
  # Observers to track the "Submit" button clicks in the "Differential Expression Results" and "Plot Results" tabs
  observeEvent(input$submit_btn_expr, {
    plot_submit_clicked_expr(FALSE)  # Reset the flag when clicking Differential Expression Results tab
  })
  
  observeEvent(input$submit_plot_btn, {
    plot_submit_clicked_expr(TRUE)  # Set the flag when clicking Plot Results tab
  })
  
  # Read the gene counts matrix
  counts_data2 <- reactive({
    req(input$counts_file)
    validate_file(input$counts_file$datapath)
    
    # Read the counts file and transpose it
    counts_matrix <- read.csv(input$counts_file$datapath, stringsAsFactors = FALSE, row.names = 1)
    counts_matrix_transposed <- t(counts_matrix)
    
    # Convert row names to a new column named "Sample"
    counts_df <- data.frame(Sample = row.names(counts_matrix_transposed), counts_matrix_transposed)
    
    # Reset row names
    row.names(counts_df) <- NULL
    
    counts_df
  })
  
  # Read the sample information matrix
  info_data <- reactive({
    req(input$info_file)
    validate_file(input$info_file$datapath)
    read.csv(input$info_file$datapath, stringsAsFactors = FALSE)
  })
  
  # Update the categorical variable selection dropdown based on the sample information data
  observe({
    req(info_data())
    variable_choices <- names(info_data())
    updateSelectInput(session, "categorical_var", choices = variable_choices[-1])  # Exclude the first column
  })
  
  # Update the gene selection dropdown based on the gene counts data
  observe({
    req(counts_data2())
    gene_choices <- names(counts_data2())[2:length(names(counts_data2()))]  # Exclude the first column (Sample)
    updateSelectInput(session, "gene_var", choices = gene_choices)
  })
  
  # Generate and display the gene expression plot
  output$gene_plot <- renderPlot({
    req(counts_data2(), info_data(), input$gene_var, input$categorical_var, input$plot_type)
    
    # Assuming the selected column for gene counts
    selected_gene <- input$gene_var
    selected_category <- input$categorical_var
    plot_type <- input$plot_type
    
    # Prepare data for plotting
    data_to_plot <- data.frame(
      Sample = counts_data2()$Sample,
      Expression = counts_data2()[[selected_gene]],
      Category = info_data()[[selected_category]]
    )
    
    # Generate plot based on selected plot type
    if (plot_type == "Bar Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_bar(stat = "identity", fill = "steelblue2") +  # Fill with light blue, border with black
        labs(title = paste("Bar Plot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Boxplot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_boxplot(fill = "olivedrab4") +
        labs(title = paste("Boxplot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Violin Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_violin(fill = "darkgoldenrod1") +
        labs(title = paste("Violin Plot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Beeswarm Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_beeswarm(color = "lightcoral", size = 3, cex = 0.7) +  # Proper beeswarm configuration
        labs(title = paste("Beeswarm Plot of", selected_gene, "by", selected_category),
             x = selected_category, y = "Expression")
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)