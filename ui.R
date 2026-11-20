#################################################
# Differential Analysis of Proteomics Data
# User Interface (UI) Components
#################################################

shinyUI(fluidPage(

  # Application header
  headerPanel(
    title = "Differential analysis of label-free proteomic experiments"
  ),

  # Sidebar panel with inputs
  wellPanel(

    # File type selection
    checkboxInput("spectronaut", "Spectronaut input file"),
    helpText("Check if you're using Spectronaut format, uncheck for MaxQuant"),

    # File upload
    fileInput('file1', 'Choose your input file ("proteinGroups.txt" from MaxQuant, "...Report.csv" from Spectronaut):',
              accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.tsv')),

    # Intensity type selection
    radioButtons('nameKind', 'Kind of intensities used for the analysis:',
                 choices = c(Raw = 'Raw', LFQ = 'LFQ'),
                 selected = 'Raw',
                 inline = TRUE),

    # Columns flagging the proteins to delete
    htmlOutput("protnameDel"),

    HTML("<p style= 'font-size : 14px;font-weight: bold;'>Important: if you want to use SAINT, the first condition has to be the control samples and the second one the bait samples. Also, the names of baits have to be insert in the column names by renaming columns Peptides, Unique peptides, Intensity, iBAQ,LFQ intensity as : \"Column_name bait_name sample_name\" (e.g. LFQ intensity CD59 CD59_1GFP).</p>"),

    checkboxInput('saint', HTML("<p style= 'font-size : 14px;font-weight: bold;'> Check if you want to use SAINT <p>")),

    textInput("nameC1", "Name of the first condition:", "Condition 1"),
    # Variable selection condition 1 :
    htmlOutput("varselectC1"),

    textInput("nameC2", "Name of the second condition:", "Condition 2"),
    # Variable selection condition 2 :
    htmlOutput("varselectC2"),

    radioButtons('nameIBAQ', 'IBAQ values:', c(None = 'No', Yes = "Yes"), 'Yes'),

    numericInput("N_Unique_Identified", "Minimal number of unique identified peptides for a quantified value:", 1),

    numericInput("numFilter", "Minimal number of observed values in at least one condition:", 2),

    numericInput("numFilter2", "Minimal percentage of observed values in at least one condition (between 0 and 100):", 0),

    radioButtons('nameNorm', 'Normalization method:', c(None = 'None', Median_Centering = 'Median'), 'Median'),

    radioButtons('nameImpute', 'Imputation method:', c(None = 'None', MLE = "MLE", SLSA = "SLSA"), 'MLE'),

    numericInput("FCT", "Select a minimal ratio (fold-change) between conditions (2 by default, i.e. log2(fold-change) threshold is fixed at log2(2)=1):", 2),

    radioButtons('nameTest', 'Statistical test:', c(LIMMA = 'LIMMA', paired.LIMMA = 'paired.LIMMA'), 'LIMMA'),

    radioButtons('nameProp', 'Estimator of the proportion of differentially abundant proteins:', c(Pounds = 'pounds', Slim = 'slim'), 'pounds'),

    checkboxInput('condP', 'Make the p-values conditional on |log2(fold-change)| > threshold before the FDR adjustment (each p-value is divided by the p-value the protein would have with |log2(fold-change)| = threshold)', value = FALSE),

    numericInput("FDRT", "FDR level:", 0.01),

    numericInput('avgp', 'AvgP threshold (only when SAINT is used):', value = 0.5),

    checkboxInput('modif', 'PPTX with modifiable plots'),

    # Sample overlap (Jaccard / overlap indices) and clustering robustness (pvclust AU/BP)
    checkboxInput('overlap', 'Add Jaccard/overlap indices and pvclust AU/BP p-values (sample clustering)', value = FALSE),
    conditionalPanel(
      condition = "input.overlap == true",
      textInput('revPrefix', 'Prefix of the output files:', 'AllProteins'),
      textInput('groupsManual',
                'Group of each sample, comma-separated, in the order of the selected samples (leave empty to use the sample name without its trailing number, e.g. Ctrl3 -> Ctrl):',
                ''),
      numericInput('nboot', 'Number of bootstrap replicates for pvclust (10000 recommended):', 10000, min = 100),
      numericInput('nperm', 'Number of permutations for Jaccard tests:', 9999, min = 99),
      numericInput('ncores', 'Number of cores for pvclust (pvclust results depend on this number):', 5, min = 1)
    ),

    numericInput('seed', 'Random seed (only used by the SLSA imputation):', 1234),

    actionButton('run', 'Run the analysis', class = "btn-primary"),
    uiOutput("download_ui")
  ),

  mainPanel(
    fluidRow(
      column(width = 6, verbatimTextOutput("stdout_output")),
      column(width = 6, verbatimTextOutput("error_output"))
    )
  )
))
