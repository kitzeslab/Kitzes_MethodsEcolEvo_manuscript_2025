
#
##
### LOAD LIBRARIES
##
#

library(shiny)
library(unmarked)
library(ggplot2)
set.seed(123)

##############################################################################

#
##
### SHINY APP
##
#

#DEFINE UI 
ui <- fluidPage(
  
  #APP TITLE
  titlePanel("Simulation of effects of error rates throughout sensors + AI pipeline"),
  
  p("This app supports the exploration of how different stages of the the sensors + AI
  pipeline affect bias and precision in single-season static occupancy models. Select parameters using the sliders,
  then click the button to simulate 10 data sets and fitted models. The parameter ranges below correspond roughly
  to those for a typical temperate bird survey in which recorders are deployed at many points
  for multiple days, although the same simulations appy to other vocalizing taxa."),
  
  p("NOTES: The models may fail to estimate upper and lower confidence interval 
  bounds for psi (returning NA) under certain combinations of parameters. This
  can occur when the simulated detection history consists either of all 0's or
  all 1's for every site (i.e., there is no variability in detection of presences
  within any site). Try adjusting parameters to lower detectability and/or the
  false positive rate to address this issue."),
  
  #SET UP SIDEBAR FIELDS
  fluidRow(
    column(6,
          # Header
          h3("Study Design"),
          # Sliders for model parameters
           sliderInput("n_sites","Number of sites surveyed", 
                       min = 10, max = 200, value = 50, step = 10),
           sliderInput("n_surveys", "Number of days of recording", 
                       min = 2, max = 10, value = 4, step = 1), 
          
          # Header
          h3("Ecological Context"),
          sliderInput("psi","Probability a site is occupied (ψ)", 
                      min = 0.0, max = 1.00, value = 0.25, step = 0.05),
           sliderInput("availability","Probability at least one bird is available for detection during a visit",
                       min = 0.1, max = 1.0, value = 1, step = 0.1),
           sliderInput("sing_rate","Probability at least one bird sings in a 5 second window, given site occupancy and daily availability", 
                       min = 0.0, max = 0.25, value = 0.05, step = 0.01)),
    
    column(6,
          h3("Stages 1-3: Hardware, Deployment, Data"),
          sliderInput("recording_length","Number of minutes recorded per day",
                       min = 5, max = 60, value = 10, step = 5),
           sliderInput("songs_captured","Probability that a song made in a 5 second window is recorded successfully", 
                       min = 0.0, max = 1.0, value = 1, step = 0.1),

          h3("Stage 4: AI"),                       
           sliderInput("p_false_pos","Probability a true negative file is predicted positive by the classifier (false positive)",
                       min = 0.0, max = 0.01, value = 0.001, step = 0.001),
           sliderInput("p_false_neg","Probability a true positive file is predicted negative by the classifier (false negative)",
                       min = 0.0, max = 0.5, value = 0.25, step = 0.05),

          h3("Stage 5: Statistics"),                       
          sliderInput("min_num_pred","Minimum number of positives needed to declare presence for visit",
                       min = 1, max = 10, value = 1, step = 1), 
           checkboxInput("review","Check this box if recordings are reviewed by a human
           (human reviewed clips are labeled with no errors)", 
                         value = FALSE),
           sliderInput("n_rev","If human review, maximum number of predicted positive clips reviewed per visit",
                       min = 1, max = 10, value = 1, step = 1),

           
           #BUTTON TO KICK OFF THE MODEL
          h3("Stage 6: Insight"),  
           actionButton("run", "Simulate & Fit Model")
    ), 
    
    mainPanel(
      br(), br(),
      verbatimTextOutput("model_output"),
      plotOutput("plot_output")
    )
  )
)

##############################################################################

#
##
### DEFINE SERVER LOGIC
##
#

server <- function(input, output) {
  
  #OCCUPANCY MODEL RESULTS
  model_result <- eventReactive(input$run, {
    results <- data.frame(
      Run = integer(),
      Estimate = numeric(),
      Lower = numeric(),
      Upper = numeric()
    )
    
    for (run in 1:10) {
      #SET MODEL INPUTS FROM USER SELECTIONS
      n_sites <- input$n_sites
      psi <- input$psi
      sing_rate <- input$sing_rate
      songs_captured <- input$songs_captured
      availability <- input$availability
      p_false_pos <- input$p_false_pos
      p_false_neg <- input$p_false_neg
      n_surveys <- input$n_surveys
      recording_length <- input$recording_length
      review <- input$review
      min_num_pred <- input$min_num_pred
      n_rev <- input$n_rev
      
      # Set up detection history matrix
      n_occup <- round(n_sites * psi)
      z <- c(rep(1, n_occup), rep(0, n_sites - n_occup))
      y <- matrix(0, nrow = n_sites, ncol = n_surveys)
      
      # Loop through each site and survey to make detection history
      for (i in 1:n_sites) {
        for (j in 1:n_surveys) {

          # If site is occupied, find number of songs and true positives
          if (z[i] == 1 && rbinom(1, 1, availability) == 1) {
              n_songs <- rbinom(1, recording_length * 12, sing_rate * songs_captured)
              n_true_pos_pred <- rbinom(1, n_songs, (1 - p_false_neg))

          # If site is unoccupied, assign zero songs and true positives
          } else {
              n_songs <- 0
              n_true_pos_pred <- 0
          }

          # Find number of false positives
          n_false_pos_pred <- rbinom(1, recording_length * 12 - n_songs, p_false_pos)

          # Find number of predicted positives (true and false)
          n_pos_pred <- n_true_pos_pred + n_false_pos_pred
          
          # If there's no review, assign site-visit present if n_pos_pred greater than minimum
          if (!review) {
            y[i,j] <- ifelse(n_pos_pred >= min_num_pred, 1, 0)
            
          # If there's review, assign site-visit present if at least one confirmed positive
          } else {
            
            # But first, if min needed positives is more than number of clips to review,
            # throw an error as you can never declare a site present
            #print(min_num_pred)
            if (min_num_pred > n_rev) {
              stop("The number of clips to review must be greater than the
                   minimum number of positives needed to declare presence.")
            }
            
            n_rev_realized <- min(n_pos_pred, n_rev)
            n_true_pos_confirmed <- sum(sample(
              c(rep(1, n_true_pos_pred), rep(0, n_false_pos_pred)),
              n_rev_realized, replace = FALSE
            ))
            y[i,j] <- ifelse(n_true_pos_confirmed >= min_num_pred, 1, 0)
          }
        }
      }
      
      # Fit occupancy model using unmarked
      umf <- unmarkedFrameOccu(y = y)
      model <- occu(~1 ~1, data = umf)
      estimate <- as.numeric(model@estimates@estimates$state@estimates)
      conf <- confint(model, type = 'state')
      
      results <- rbind(results, data.frame(
        Run = run,
        Naive = sum(rowSums(y) > 0),
        Estimate = plogis(estimate),
        Lower = plogis(conf[1]),
        Upper = plogis(conf[2])
      ))
    }
    
    list(
      results = results,
      psi = psi)
    
  })
  
  output$model_output <- renderPrint({
    print(model_result()$results)
  })
  
  output$plot_output <- renderPlot({
    model_out <- model_result()
    df <- model_out$results
    psi <- model_out$psi
    
    ggplot(df, aes(x = factor(Run), y = Estimate)) +
      geom_point(size = 3, color = "steelblue") +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray30") +
      geom_hline(yintercept = psi, linetype = "dashed", color = "red") +
      ylim(0, 1) +
      theme_minimal() +
      labs(
        x = "Run",
        y = "Occupancy Estimate",
        title = "Occupancy Estimates Across 10 Simulations"
      )
  })
  
}

##############################################################################

#
##
### RUN THE APP
##
#

shinyApp(ui = ui, server = server)
