library(shiny)
library(ggplot2)

ll<-function(theta,x,n) {
    dbinom(x=x,size=n,prob=theta,log=TRUE)-dbinom(x=x,size=n,prob=x/n,log=TRUE)
}

qa<-function(theta,mle,se) {
    -0.5 * ((mle - theta)^2 / se^2)
}

wald_test <- function(x, n, theta0) {
  p_hat <- x / n
  se <- sqrt(p_hat * (1 - p_hat) / n)
  z <- (p_hat - theta0) / se
  p_value <- 2 * (1 - pnorm(abs(z)))
  return(p_value)
}

score_test <- function(x, n, theta0) {
  se <- sqrt(theta0 * (1 - theta0) / n)
  z <- (x - n * theta0) / (n * se)
  p_value <- 2 * (1 - pnorm(abs(z)))
  return(p_value)
}

lr_test <- function(x, n, theta0) {
  p_hat <- x / n
  log_likelihood_null <- x * log(theta0) + (n - x) * log(1 - theta0)
  log_likelihood_alt <- x * log(p_hat) + (n - x) * log(1 - p_hat)
  lr_stat <- -2 * (log_likelihood_null - log_likelihood_alt)
  p_value <- 1 - pchisq(lr_stat, df = 1)
  return(p_value)
}

ui <- fluidPage(
  withMathJax(),
    titlePanel("Binomial likelihood: Trinity of asymptotic tests"),
    sidebarLayout(          
        sidebarPanel(
            titlePanel("Data and Hypothesis"),
            sliderInput(inputId="x",label="Data: Number of successes",min=50,max=350,value=100,step=5),
            sliderInput(inputId="n",label="Data: Sample size",min=0,max=400,step=5,value=300),
            sliderInput(inputId="theta_0",label="Constraint (Null hypothesis)",min=0.1,max=0.9,value=0.5,step=0.01)
        ),
        mainPanel(
            titlePanel("Wald, Score and LR-Test"),
            plotOutput(outputId="LL",width="50%"),
            )))




server <- function(input, output,session) {

 observeEvent(input$x, {
    updateSliderInput(inputId = "n", min = input$x)
 })

    
    
    output$LL<-renderPlot({
        n<-input$n
        mle<-input$x/input$n
        se<-sqrt(mle*(1-mle)/input$n)
        theta_seq<-seq(0.0,1,.01)
        # Calculate log-likelihood values
        ll_values <- ll(theta_seq,input$x,input$n)
        ll_MLE <- ll(mle,input$x,input$n)
        qa_values <- qa(theta_seq,mle=mle,se=se)
        qa_MLE <- qa(mle,mle=mle,se=se)
        ll_null <- ll(input$theta_0,input$x,input$n)
        qa_null <- qa(input$theta_0,mle=mle,se=se)
        slope <- input$x/input$theta_0-(input$n-input$x)/(1-input$theta_0)
        tangent_range <- seq(input$theta_0 - 0.1, input$theta_0 + 0.1, length.out = 100)
        tangent_line <- ll_null + slope * (tangent_range - input$theta_0)
        data   <- data.frame(theta = theta_seq,ll = ll_values)
        dataqa <- data.frame(theta = theta_seq,qa = qa_values)
        tangent_data <- data.frame(theta = tangent_range,ll = tangent_line)
        ggplot(data, aes(x = theta, y = ll))+ ylim(-150,10)+xlim(0,1)+
            geom_line() +
            geom_vline(xintercept = mle, linetype = "dashed", color = "blue") +
            geom_vline(xintercept = input$theta_0, linetype = "dashed", color = "red") +
            geom_hline(yintercept = 0, linetype = "dotted",color="brown") +
            geom_hline(yintercept = ll_null, linetype = "dotted",color="brown") +
            geom_segment(aes(x = mle, y = qa_null, xend = input$theta_0, yend = qa_null),
                         arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm")), linetype="dashed",color = "blue",lwd=1.5) +
            geom_segment(data = tangent_data, aes(x = theta[1], y = ll[1], xend = theta[100], yend = ll[100]),  lwd=1.5,color = "red") +
            geom_segment(aes(x = input$theta_0, y = ll_null, xend = input$theta_0, yend = ll_MLE),
                         arrow = arrow(type = "closed", ends = "both", length = unit(0.2, "cm")), lwd=1.5,color = "brown") +
            labs(x = expression(theta),y = "Log-Likelihood (black), quadratic approximation (dashed)") +
            scale_x_continuous() +
             ## scale_x_continuous(breaks = c(input$theta_0, mle), labels = c(expression(theta[0]), expression(hat(theta)[MLE])))+ 
            annotate("text", x = mle, y = qa_null, label = expression(Wald),color = "blue", vjust = -1, hjust = 1.5) +
            annotate("text", x = input$theta_0, y = ll_null, label = expression(Score), color = "red", vjust = -1, hjust = -0.5)+
            annotate("text", x = input$theta_0, y = (ll_null + ll_MLE) / 2, label = expression(LR), color = "brown", vjust = -1, hjust = -0.5)+
            geom_line(data=dataqa,aes(x=theta,y=qa),linetype="dashed")
    })
}    

shinyApp(ui = ui, server = server)


