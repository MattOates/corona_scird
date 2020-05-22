#
# Shiny app to explore SCIRD models for looking at Lockdown
#

library(shiny)
library(deSolve)

# Function to return derivatives of SCIRD model
seir_ode <- function(t,state,parameters) {
    # With our state and parameter names
    with(as.list(c(state, parameters)), {

    # ODE system definition
    
    # Our susceptible pop
    dSdt <- lambda + q * C - p * S - (beta/N) * I * S
    
    # Our confined pop
    # ToDo: Add in an exponential? ramp up for q given t 
    #       to capture lockdown fatigue. Could make it depend
    #       on deaths/infected too for informed pop
    dCdt <- p * S - q * C
    
    # Our infected pop
    dIdt <- (beta/N) * S * I - (gamma + mu) * I
    
    # The dead
    dDdt <- mu * I
    
    # The recovered
    dRdt <- gamma * I

    # List of equations to return
    list(c(dSdt, dCdt, dIdt, dDdt, dRdt))
    
    })
}

# Set parameter values
beta<-520/365;
# Births per day (2018 averaged)
lambda<-731213/365
# Total population, default UK
N<-65640000;
# Rate of recovery from infection
gamma<-1/30;
# Death rate
#mu<-774835/(N*365) # UK birth and population figures 2016
mu<-0.01
# Proportion locking down per day
q<-0.8
# Proportion breaking lock down per day
p<-0.2
# Time range, assume a year
t<-seq(0,365)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Interactive SCIRD Model"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("N",
                        "N (Total population):",
                        min = 1,
                        max = 7594000000,
                        value = 65640000),
            sliderInput("beta",
                        "Beta β (Contact rate per person, per unit time):",
                        min = 1,
                        max = 1000,
                        value = 520/365),
            sliderInput("gamma",
                        "Gamma (Rate of recovery from infection):",
                        min = 0.0,
                        max = 1.0,
                        value = 1/30),
            sliderInput("mu",
                        "Mu (Death rate from infection):",
                        min = 0.0,
                        max = 1.0,
                        value = 774835/(65640000*365)),
            sliderInput("q",
                        "q (Proportion locking down):",
                        min = 0.0,
                        max = 1.0,
                        value = 0.9),
            sliderInput("p",
                        "p (Proportion coming out of lock down for work):",
                        min = 0.0,
                        max = 1.0,
                        value = 0.2)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        
        # Bundle up params
        # ToDo: Get these from Shiny
        parameters<-c(beta=beta,gamma=gamma,mu=mu,q=q,p=p,N=N,lambda=lambda)
        # Initial state of one infected and total pop susceptible
        initial_state<-c(S=N,C=0,I=1,R=0,D=0)
        
        # Solve system using lsoda
        sol<-lsoda(initial_state,t,seir_ode,parameters)
        
        # Plot solution
        plot(t,sol[,2],type="l",col="blue",ylim=c(0,N),ylab="Population")
        lines(t,sol[,3],col="orange")
        lines(t,sol[,4],col="red")  
        lines(t,1-rowSums(sol[,2:4]),col="green")
        legend(
            300,0.7*N,
            legend=c("S","C","I","R","D"),
            col=c("blue","orange","red","green","purple"), 
            lty=1, 
            cex=0.8
        )
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
