library(shiny)
library(shinydashboard)
library(DT)

# Define data
NAM_Phenotype <- c("ALB","ALP",'ALT',"APOA","APOB","AST","TBL","BMI","BUN","CA","CRP","CysC","DBP", 
                   "FPG","GGT","HB","HT","HBA1C","HDL","HEG","LDL","LYM","Mono","Neutro","IGF1", 
                   "PLA","PP","RBC","SBP","SHBG","RET","TCh","TG","TT","TP","UA","VTD","WBC","CAD",
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
  PP = "Pulse Pressure",
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
                                       choices = c("Age Interaction Test" = "Age",
                                                   "Retirement RDD/RKD Test" = "Retirement"),
                                       selected = "Age")
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
                      uiOutput("plot_area")
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
  
  # 添加静态文件资源路径
  addResourcePath("images", ".")
  
  # Update trait selection with search functionality
  observe({
    updateSelectizeInput(session, "trait_select",
                         choices = trait_choices,
                         server = TRUE)
  })
  
  # Generate web-accessible image path(s)
  web_image_path <- reactive({
    if (input$trait_select != "") {
      trait_code <- input$trait_select
      
      if (input$analysis_type == "Age") {
        # Age 目录下二选一：优先 *_plot.png，其次 *_plot_quan.png
        cand <- c(
          paste0("Age/", trait_code, "_plot.png"),
          paste0("Age/", trait_code, "_plot_quan.png")
        )
        
        for (path in cand) {
          if (file.exists(path)) {
            return(paste0("images/", path))
          }
        }
        return(NA_character_)
        
      } else {
        # Retirement：并排两张
        p1 <- paste0("Retirement/", trait_code, "_RDD_quan.png")
        p2 <- paste0("Retirement/", trait_code, "_RKD_quan.png")
        
        web_p1 <- if (file.exists(p1)) paste0("images/", p1) else NA_character_
        web_p2 <- if (file.exists(p2)) paste0("images/", p2) else NA_character_
        
        return(c(web_p1, web_p2))
      }
    } else {
      return(NULL)
    }
  })
  
  # Output plot title
  output$plot_title <- renderText({
    if (input$trait_select != "") {
      trait_name <- phenotype_fullnames[input$trait_select]
      analysis_name <- ifelse(input$analysis_type == "Age", 
                              "Age Interaction Test", 
                              "Retirement RDD/RKD Test")
      paste(analysis_name, "-", trait_name)
    }
  })
  
  # Display image(s): Age 单图；Retirement 并排两图
  output$plot_area <- renderUI({
    if (input$trait_select == "") return(NULL)
    paths <- web_image_path()
    
    # 占位图策略
    placeholder <- NULL
    if (file.exists("www/placeholder.png")) {
      placeholder <- "placeholder.png"  # www目录下的文件可以直接访问
    } else if (file.exists("placeholder.png")) {
      placeholder <- "images/placeholder.png"
    }
    
    if (input$analysis_type == "Age") {
      src <- if (!is.null(paths) && length(paths) >= 1 && !is.na(paths[1])) {
        paths[1]
      } else {
        placeholder
      }
      
      if (is.null(src)) {
        return(tags$div(style="color:#999; padding: 100px 0; text-align: center;",
                        h4("Image not found")))
      }
      
      tags$img(src = src, alt = "Analysis Plot",
               style = "max-width: 100%; height: 540px; object-fit: contain;")
      
    } else {
      # Retirement：两列（RDD | RKD）
      p1 <- if (length(paths) >= 1 && !is.na(paths[1])) paths[1] else placeholder
      p2 <- if (length(paths) >= 2 && !is.na(paths[2])) paths[2] else placeholder
      
      fluidRow(
        column(
          6,
          tags$div(
            style = "text-align: center;",
            if (is.null(p1)) {
              tags$div(style="color:#999; padding: 50px 0;", h5("RDD image not found"))
            } else {
              tags$img(src = p1, alt = "RDD Plot",
                       style = "width: 100%; height: 540px; object-fit: contain;")
            }
          )
        ),
        column(
          6,
          tags$div(
            style = "text-align: center;",
            if (is.null(p2)) {
              tags$div(style="color:#999; padding: 50px 0;", h5("RKD image not found"))
            } else {
              tags$img(src = p2, alt = "RKD Plot",
                       style = "width: 100%; height: 540px; object-fit: contain;")
            }
          )
        )
      )
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)