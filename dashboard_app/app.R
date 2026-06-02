library(shiny)
ui <- fluidPage(
  h1("BRCA Navigator - Test Mode"),
  p("If you see this text, Shiny Server + Docker works."),
  textOutput("status")
)
server <- function(input, output, session) {
  output$status <- renderText("Server session is running. Time: Sys.time()")
}
shinyApp(ui, server)
