# Global variables can go here
library(stringr)
.fix <- function(x){
  x <- tolower(x)
  x <- str_replace_all(x, "[[:punct:]]", "_")
  x <- str_replace_all(x, " ", "_")
  return(x)
}


# Define the UI
ui <- fluidPage(
  # Application title
  titlePanel("Create project name"),

  sidebarLayout(
    # Sidebar with a slider and selection inputs
    sidebarPanel(
      textInput('hbc', 'hbc-code (no letters)', value = "00000"),
      textInput('pi', 'What is PI last name:', value = "lastname"),
      textInput('scientist', 'What is the scientist last name:', value = "scientist"),
      textInput('tech', 'What is the technology:', value = "rnaseq"),
      textInput('tissue', 'What is the tissue:', value = "mix|cells|heart"),
      textInput('org', 'What is the organism:', value = "mix|human"),
      textInput('proj', 'What is the project name:', value = "this_analysis_is_cool"),

    ),

    # Show Word Cloud
    mainPanel(
      br("Suggested project name:"),
      br(),
      verbatimTextOutput('project')
    )
  )
)


# Define the server code
server <- function(input, output, session) {
  output$project <- renderText({
    hbc_code <- paste0("hbc", input$hbc)
    project_full <- paste(input$tech, .fix(input$pi), .fix(input$scientist),
                          .fix(input$proj),
                          input$tissue, input$org, hbc_code, sep="_")
    project_full
  })
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)
