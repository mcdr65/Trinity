library(shiny)
library(ggplot2)
# Define the log-likelihood function
log_likelihood <- function(MLE,theta, sigma) {
  -0.5 * ((MLE - theta)^2 / sigma^2)
}
# Parameters
theta_MLE <- 0
sigma <- 1
theta_seq <- seq(-3, 3, length.out = 100)
# Calculate log-likelihood values
ll_values <- log_likelihood(theta_seq, theta_MLE, sigma)
ll_MLE <- log_likelihood(theta_MLE, theta_MLE, sigma)



ui <- fluidPage(
    titlePanel("One-parameter model: Trinity of asymptotic tests"),
    verticalLayout(sliderInput(inputId="theta_0",label="Null Hypothesis (Constraint)",min=-2,max=2,value=0,step=0.1),
                    plotOutput(outputId="LL",width="50%")))


server <- function(input, output) {
    output$LL<-renderPlot( {
        ll_null <- log_likelihood(input$theta_0, theta_MLE, sigma)
        slope <- (theta_MLE - input$theta_0) / sigma^2
        tangent_range <- seq(input$theta_0 - 0.5, input$theta_0 + 0.5, length.out = 100)
        tangent_line <- ll_null + slope * (tangent_range - input$theta_0)
        data <- data.frame(theta = theta_seq,ll = ll_values)
        tangent_data <- data.frame(theta = tangent_range,ll = tangent_line)
        ggplot(data, aes(x = theta, y = ll)) +
            geom_line() +
            geom_vline(xintercept = theta_MLE, linetype = "dashed", color = "blue") +
            geom_vline(xintercept = input$theta_0, linetype = "dashed", color = "red") +
            geom_hline(yintercept = 0, linetype = "dotted",color="black") +
            geom_hline(yintercept = ll_null, linetype = "dotted",color="black") +
            geom_segment(aes(x = theta_MLE, y = -3, xend = input$theta_0, yend = -3),
                         arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm")), color = "blue",lwd=1.5) +
            geom_segment(data = tangent_data, aes(x = theta, y = ll, xend = theta + 0.02, yend = ll + 0.01 * slope),  lwd=1.5,color = "red") +
            geom_segment(aes(x = input$theta_0+.5, y = ll_null, xend = input$theta_0+.5, yend = ll_MLE),
                         arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm")), lwd=1.5,color = "black") +
            labs(x = expression(theta),y = "Log-Likelihood") +
            scale_x_continuous(breaks = c(input$theta_0, theta_MLE), labels = c(expression(theta[0]), expression(hat(theta)[MLE]))) +
            annotate("text", x = theta_MLE+1, y = -3, label = expression(Wald~Test), color = "blue", vjust = -1, hjust = 1.5) +
            annotate("text", x = input$theta_0+0.1, y = ll_null-0.5, label = expression(Score~Test), color = "red", vjust = -1, hjust = -0.5) +
            annotate("text", x = input$theta_0+.5, y = (ll_null + ll_MLE) / 2, label = expression(LR~Test), color = "black", vjust = -1, hjust = -0.5)
    })
}    

shinyApp(ui = ui, server = server)


