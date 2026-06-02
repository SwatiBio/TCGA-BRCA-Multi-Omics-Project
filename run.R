port <- as.integer(Sys.getenv("PORT", "7860"))
shiny::runApp("/srv/shiny-server", host = "0.0.0.0", port = port, launch.browser = FALSE)
