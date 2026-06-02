port <- as.integer(Sys.getenv("PORT", "7860"))
shiny::runApp("/app", host = "0.0.0.0", port = port, launch.browser = FALSE)
