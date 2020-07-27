## Piers MacNaughton - 2020
#VRE Colonization - Six Compartment Deterministic Model
#Based on paper by Wolkewitz et al. 2008
#With additional term for daylight decontamination


library(shiny)
library(gridExtra)
library(tidyverse)
library(deSolve)
library(EpiModel)
library(ggplot2)

## Configuration
theme_set(theme_minimal(base_size = 18))


hai_fun <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {
        ##main ODEs
        
        #Population sizes
        n.p <- c.p + u.p
        n.s <- c.s + u.s
        n.e <- c.e + u.e
        
        #Contaminated Patients
        dc.p <- (lam.u*n.p+(lam.c-lam.u)*c.p)*phi + (beta.sp*c.s/n.s+beta.ep*c.e/n.e)*(n.p-c.p) - lam.c*c.p
        
        #Uncontaminated Patient
        du.p <- (lam.u*n.p+(lam.c-lam.u)*c.p)*(1-phi) - lam.u*(n.p-c.p) - (beta.sp*c.s/n.s+beta.ep*c.e/n.e)*(n.p-c.p)
        
        #Contaminated HCW
        dc.s <- (beta.ps*c.p/n.p+beta.es*c.e/n.e)*(n.s-c.s) - mu*c.s
        
        #Uncontaminated HCW
        du.s <- mu*c.s - (beta.ps*c.p/n.p+beta.es*c.e/n.e)*(n.s-c.s)
        
        #Contaminated Environment
        dc.e <- (beta.se*c.s/n.s+beta.pe*c.p/n.p)*(n.e-c.e) - kappa*c.e - alpha*c.e
        
        #Uncontamintated Environment
        du.e <- kappa*c.e + alpha*c.e - (beta.se*c.s/n.s+beta.pe*c.p/n.p)*(n.e-c.e)
        
        ##output
        
        list(c(dc.p, du.p, dc.s, du.s, dc.e, du.e), n.p=n.p, Prevalence=(c.p/n.p))
        
    })
}

solve_ode <- function(alpha_one, alpha_two){
    
    param <- param.dcm(phi=0.1, lam.u=0.1/24, lam.c=0.05/24, mu=24/24, kappa=1/24, alpha=c(0,alpha_one/24,alpha_two/24),
                       beta.sp=0.3/24, beta.se=2/24, beta.ps=2/24, beta.pe=2/24, beta.es=2/24, beta.ep=0.3/24)
    
    init <- init.dcm(c.p=1, u.p=19, c.s=0, u.s=5, c.e=0, u.e=100)
    
    control <- control.dcm(nsteps=120*24, new.mod=hai_fun)
    
    mod_hai <- dcm(param, init, control)
    
    mod_output <- as.data.frame(mod_hai)
    mod_output$Day <- (mod_output$time+24)/24
    
    mod_output
}


#Plot of Prevalence of VRE Colonization
plot_result <- function(mod_output, alpha_one, alpha_two, max_time) {    
    solve_ode(alpha_one, alpha_two)
    
    pp <- ggplot(data=mod_output, aes(x=Day, y=Prevalence)) + geom_point(aes(color=as.factor(run))) + xlim(c(0,max_time)) +
          scale_color_manual(name="Decontamination Rate", labels=c("None", "Low", "High"), values=c("#000000", "#EBB172", "#EB9234")) +
          theme(legend.title = element_text(size=16))
    
    print(pp)
}


ui <- fluidPage(
    
    titlePanel("VRE Transmission in Hospital - Deterministic Compartment Model"),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("x_max", 
                        "Max days for the model", 
                        min   = 3, 
                        max   = 120, 
                        value = 30,
                        step  =  1),
            hr(),
            br(),
            sliderInput("alpha_one",
                        "Daylight Decontamination Rate (/day) - Low",
                        min   = 0,
                        max   = 1.0,
                        value = 0.5, 
                        step  = 0.1), 
            sliderInput("alpha_two",
                        "Daylight Decontamination Rate (/day) - High",
                        min   = 1,
                        max   = 10,
                        value = 5, 
                        step  = 1), 
            br(),
            hr(),         
            tags$div("Piers MacNaughton, 2020"), 
        ),
        
        mainPanel(
            # p("Note: This tool is not intended to create a prediction."),
            plotOutput("chart", height = "500px"), 
        )
    )
)


server <- function(input, output) {
    
    res <- reactive({
        solve_ode(alpha_one=input$alpha_one, alpha_two=input$alpha_two)
    })
    
    
    output$chart <- renderPlot({
        plot_result(res(), input$alpha_one, input$alpha_two, input$x_max)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
