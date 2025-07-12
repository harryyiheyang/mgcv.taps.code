library(shiny)
library(shinydashboard)
library(DT)

# Define data
NAM_Phenotype <- c("ALB","ALP",'ALT',"APOA","APOB","AST","TBL","BMI","BUN","CA","CRP","CysC","DBP", 
                   "FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1", 
                   "PLA","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","UA","VTD","WBC","CAD",
                   "Stroke","HF","AF","PAD","T2D","SMK","DRNK")

phenotype_fullnames <- c(
  ALB = "Albumin", 
  ALP = "Alkaline Phosphatase", 
  ALT = "Alanine Aminotransferase", 
  APOA = "Apolipoprotein A", 
  APOB = "Apolipoprotein B", 
  AST = "Aspartate Aminotransferase", 
  TBL = "Total Bilirubin", 
  BMI = "Body Mass Index", 
  BUN = "Blood Urea Nitrogen", 
  CA = "Calcium", 
  CRP = "C-reactive Protein", 
  CysC = "Cystatin C", 
  DBP = "Diastolic Blood Pressure", 
  FPG = "Fasting Plasma Glucose", 
  GGT = "Gamma-Glutamyl Transferase", 
  HB = "Hemoglobin", 
  HT = "Hematocrit", 
  HBA1C = "Hemoglobin A1c", 
  HDL = "HDL Cholesterol", 
  HEG = "Height", 
  LDL = "LDL Cholesterol", 
  LYM = "Lymphocyte Count", 
  Mono = "Monocyte Count", 
  Neutro = "Neutrophil Count", 
  IGF1 = "Insulin-like Growth Factor 1", 
  PLA = "Platelet Count", 
  RBC = "Erythrocyte Count", 
  SBP = "Systolic Blood Pressure", 
  SHBG = "Sex Hormone-Binding Globulin", 
  RET = "Reticulocyte Count", 
  TCh = "Total Cholesterol", 
  TG = "Triglycerides", 
  TT = "Total Testosterone", 
  TP = "Total Protein", 
  UA = "Uric Acid", 
  VTD = "Vitamin D", 
  WBC = "Leukocytes Count", 
  CAD = "Coronary Artery Disease", 
  Stroke = "Ischemic Stroke", 
  HF = "Heart Failure", 
  AF = "Atrial Fibrillation", 
  PAD = "Peripheral Artery Disease", 
  T2D = "Type 2 Diabetes", 
  SMK = "Smoking Initiation", 
  DRNK = "Drinks per Week"
)

# Create selection options
trait_choices <- setNames(NAM_Phenotype, phenotype_fullnames)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Visualization"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Visualization", tabName = "analysis", icon = icon("chart-bar"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "analysis",
              fluidRow(
                box(
                  title = "Visualization of Results of TAPS", 
                  status = "primary", 
                  solidHeader = TRUE,
                  width = 12,
                  
                  fluidRow(
                    column(6,
                           h4("Analysis Type"),
                           selectInput("analysis_type", 
                                       label = NULL,
                                       choices = c("PRS Linearity Test" = "PRS",
                                                   "Age-PRS Linear Interaction Test" = "Age"),
                                       selected = "PRS")
                    ),
                    
                    column(6,
                           h4("Trait"),
                           selectizeInput("trait_select",
                                          label = NULL,
                                          choices = trait_choices,
                                          selected = "ALB",
                                          options = list(
                                            placeholder = "Search trait...",
                                            onInitialize = I('function() { this.setValue(""); }')
                                          ))
                    )
                  ),
                  
                  hr(),
                  
                  conditionalPanel(
                    condition = "input.trait_select != ''",
                    div(
                      style = "text-align: center;",
                      h4(textOutput("plot_title")),
                      br(),
                      imageOutput("analysis_plot", height = "540px"),
                      br(),
                      p(em("Image path: "), textOutput("image_path", inline = TRUE))
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "input.trait_select == ''",
                    div(
                      style = "text-align: center; margin-top: 50px;",
                      h4("Please select a trait to view analysis results", style = "color: #999;")
                    )
                  )
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Update trait selection with search functionality
  observe({
    updateSelectizeInput(session, "trait_select",
                         choices = trait_choices,
                         server = TRUE)
  })
  
  # Generate image path
  image_path <- reactive({
    if (input$trait_select != "") {
      folder <- input$analysis_type
      trait_code <- input$trait_select
      paste0(folder, "/", trait_code, "_plot.png")
    } else {
      NULL
    }
  })
  
  # Output plot title
  output$plot_title <- renderText({
    if (input$trait_select != "") {
      trait_name <- phenotype_fullnames[input$trait_select]
      analysis_name <- ifelse(input$analysis_type == "PRS", 
                              "PRS Linearity Test", 
                              "Age-PRS Linear Interaction Test")
      paste(analysis_name, "-", trait_name)
    }
  })
  
  # Output image path
  output$image_path <- renderText({
    image_path()
  })
  
  # Display image
  output$analysis_plot <- renderImage({
    if (input$trait_select != "") {
      path <- image_path()
      
      # Check if file exists
      if (file.exists(path)) {
        list(src = path,
             alt = "Analysis Plot",
             style = "max-width: 100%; height: 100%;")
      } else {
        # Return placeholder image if file doesn't exist
        list(src = "placeholder.png",
             alt = "Image not found",
             style = "max-width: 100%; height: 100%;")
      }
    }
  }, deleteFile = FALSE)
}

# Run the app
shinyApp(ui = ui, server = server)