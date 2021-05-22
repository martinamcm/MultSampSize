
library('shinydashboard')
library('shiny')
library('DT')
library('ggplot2')
library('shinycssloaders')
library('data.table')
library('dplyr')
library('formattable')
library('tidyr')
library('ggpubr')
library('caTools')
library('knitr')
library('stats')
library('MASS')
library('shinyjs')
library('xml2')
library('rvest')
library('rmarkdown')
library('mvtnorm')
library('tableschema.r')

#xml2::write_html(rvest::html_node(xml2::read_html("OneCtsOneBinSamp.html"), "body"), file = "OneCtsOneBinSamp.html")

source('functions.R',local=TRUE)

ui <- dashboardPage(
  dashboardHeader(title = "MultSampSize"),
  
  dashboardSidebar(
    sidebarMenu(id="tab",
                menuItem("Home", tabName = "home", icon = icon("home")),
                menuItem("Sample size", tabName = "sampsize", icon = icon("calculator"), 
                   menuSubItem("Co-primary", tabName="coprim", icon=icon("sitemap")),
                   menuSubItem("Multiple Primary", tabName = "multprim", icon = icon("arrows-alt")),
                   menuSubItem("Composite", tabName="comp", icon=icon("project-diagram"))
                ),
                menuItem("Source code", icon = icon("file-code"),href="https://github.com/martinamcm/MultSampSize")
    )
  ),
  
  dashboardBody(
   useShinyjs(),
    tabItems(
      tabItem(tabName = "home",
              fluidRow(
                includeHTML("LandPageSamp.html")
              )),
      
      tabItem(tabName = "multprim",
              fluidRow(
                box(title="Endpoint",width=3,solidHeader = TRUE,status="primary",    
                    numericInput("Ctsnomult", "Number of continuous outcomes", value=1, min=1, max=2),
                    numericInput("Binnomult", "Number of binary outcomes", value=1, min=0, max=1),
                    actionButton("GetModelMult","Generate Model")
                ),
                
                box(title="Model Summary",width=8,solidHeader = TRUE,status="primary",
                    uiOutput("markdownMult")
                ),
                
                box(title="Model Parameter Inputs",width=11,solidHeader = TRUE,status="primary",
                    column(width=3,
                           numericInput("MeanY1mult", label=HTML("&delta; <sub>1</sub>"),value=1,min=0, max=100),
                           
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==0",
                                            numericInput("MeanY2mult", label=HTML("&delta; <sub>2</sub>"),value=1,min=0, max=100)
                           ),
                           conditionalPanel("input.Ctsnomult==1 && input.Binnomult==1",
                                            numericInput("piT2multa", label=HTML("&pi; <sub>T2</sub><sub style='position: relative; left: -.5em;'>2</sub>"),value=0.6,min=0, max=1, step=0.05),
                                            numericInput("piC2multa", label=HTML("&pi; <sub>C2</sub><sub style='position: relative; left: -.5em;'>2</sub>"),value=0.5,min=0, max=1, step=0.05)
                                            
                           ),
                           
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==1",
                                            numericInput("MeanY2multb", label=HTML("&delta; <sub>2</sub>"),value=1,min=0, max=100),
                                            numericInput("piT3mult", label=HTML("&pi; <sub>T3</sub><sub style='position: relative; left: -.5em;'>3</sub>"),value=0.6,min=0, max=1, step=0.05),
                                            numericInput("piC3mult", label=HTML("&pi; <sub>C3</sub><sub style='position: relative; left: -.5em;'>3</sub>"),value=0.5,min=0, max=1, step=0.05)
                           )
                           
                    ),
                    
                    column(width=3, 
                           numericInput("SigmaY1mult", label=HTML("&sigma; <sub>1</sub>"), value=1,min=0, max=100),
                           
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==0",
                                            numericInput("SigmaY2mult", label=HTML("&sigma; <sub>2</sub>"), value=1,min=0, max=100)                 
                           ),
                           
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==1",
                                            numericInput("SigmaY2multa", label=HTML("&sigma; <sub>2</sub>"), value=1,min=0, max=100)               
                           )
                    ),
                    
                    column(width=4,
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==0",
                                            sliderInput("rho12mult",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1)                
                           ),
                           
                           conditionalPanel("input.Ctsnomult==1 && input.Binnomult==1",
                                            sliderInput("rho12multa",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1)             
                           ),
                           
                           conditionalPanel("input.Ctsnomult==2 && input.Binnomult==1",
                                            sliderInput("rho12multb",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1),
                                            sliderInput("rho13mult",label=HTML("&rho; <sub>13</sub>"),value=0,min=-1,max=1,step=0.1),
                                            sliderInput("rho23mult",label=HTML("&rho; <sub>23</sub>"),value=0,min=-1,max=1,step=0.1)              
                           )
                           
                    )),  
                
                box(title="Sample Size Estimation",width=11,solidHeader = TRUE,status="primary",    
                    column(width=4,    
                           selectInput("alphamult", "One-Sided Significance Level",c("Alpha = 0.01", "Alpha = 0.05", "Alpha = 0.10"),selected="Alpha = 0.05"),
                           sliderInput("targetmult", "Power Target", value = .8, min = 0, max = 1),
                           numericInput("maxnmult", "Maximum Number of Subjects", value = 400, min = 0, max = 100000),
                           infoBoxOutput("samplesizemult",width=11)
                    ),
                    
                    column(width=7,
                           uiOutput("plotmult")
                    )
                )
              )
         ),
      
      tabItem(tabName="coprim",
              fluidRow(
                box(title="Endpoint",width=3,solidHeader = TRUE,status="primary",    
                           numericInput("Ctsnoco", "Number of continuous outcomes", value=1, min=1, max=2),
                           numericInput("Binnoco", "Number of binary outcomes", value=1, min=0, max=1),
                           actionButton("GetModel","Generate Model")
                    ),
            
                box(title="Model Summary",width=8,solidHeader = TRUE,status="primary",
                    uiOutput("markdown")
                ),
              
              box(title="Model Parameter Inputs",width=11,solidHeader = TRUE,status="primary",
                  column(width=3,
                         numericInput("MeanY1", label=HTML("&delta; <sub>1</sub>"),value=1,min=0, max=100),
                        
                          conditionalPanel("input.Ctsnoco==2 && input.Binnoco==0",
                                          numericInput("MeanY2", label=HTML("&delta; <sub>2</sub>"),value=1,min=0, max=100)
                         ),
                          conditionalPanel("input.Ctsnoco==1 && input.Binnoco==1",
                                           numericInput("piT2a", label=HTML("&pi;<sub>T2</sub><sub style='position: relative; left: -.5em;'>2</sub>"),value=0.6,min=0, max=1, step=0.05),
                                           numericInput("piC2a", label=HTML("&pi;<sub>C2</sub><sub style='position: relative; left: -.5em;'>2</sub>"),value=0.5,min=0, max=1, step=0.05)
                                           ),
                          
                          conditionalPanel("input.Ctsnoco==2 && input.Binnoco==1",
                                           numericInput("MeanY2b", label=HTML("&delta; <sub>2</sub>"),value=1,min=0, max=100),
                                           numericInput("piT3", label=HTML("&pi; <sub>T3</sub><sub style='position: relative; left: -.5em;'>3</sub>"),value=0.6,min=0, max=1, step=0.05),
                                           numericInput("piC3", label=HTML("&pi; <sub>C3</sub><sub style='position: relative; left: -.5em;'>3</sub>"),value=0.5,min=0, max=1, step=0.05)
                          )
                          
                  ),
                  
                  column(width=3, 
                         numericInput("SigmaY1", label=HTML("&sigma; <sub>1</sub>"), value=1,min=0, max=100),
                         
                         conditionalPanel("input.Ctsnoco==2 && input.Binnoco==0",
                                          numericInput("SigmaY2", label=HTML("&sigma; <sub>2</sub>"), value=1,min=0, max=100)                 
                         ),
                         
                         conditionalPanel("input.Ctsnoco==2 && input.Binnoco==1",
                                    numericInput("SigmaY2a", label=HTML("&sigma; <sub>2</sub>"), value=1,min=0, max=100)               
                         )
                         ),
                  
                  column(width=4,
                         conditionalPanel("input.Ctsnoco==2 && input.Binnoco==0",
                                    sliderInput("rho12",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1)                
                         ),
                         
                         conditionalPanel("input.Ctsnoco==1 && input.Binnoco==1",
                                          sliderInput("rho12a",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1)             
                         ),
                         
                         conditionalPanel("input.Ctsnoco==2 && input.Binnoco==1",
                                          sliderInput("rho12b",label=HTML("&rho; <sub>12</sub>"),value=0,min=-1,max=1,step=0.1),
                                          sliderInput("rho13",label=HTML("&rho; <sub>13</sub>"),value=0,min=-1,max=1,step=0.1),
                                          sliderInput("rho23",label=HTML("&rho; <sub>23</sub>"),value=0,min=-1,max=1,step=0.1)              
                          )
                         
                         )),  
                
              box(title="Sample Size Estimation",width=11,solidHeader = TRUE,status="primary",    
                    column(width=4,    
                        selectInput("alphacop", "One-Sided Significance Level",c("Alpha = 0.01", "Alpha = 0.05", "Alpha = 0.10"),selected="Alpha = 0.05"),
                         sliderInput("targetcop", "Power Target", value = .8, min = 0, max = 1),
                         numericInput("maxncop", "Maximum Number of Subjects", value = 400, min = 0, max = 100000),
                        infoBoxOutput("samplesizecopr",width=11)

                  ),
                  
                  column(width=7,
                         uiOutput("plotcoprim")
                         )
                  
                )
              
             )),
      
     
               tabItem(tabName="comp",
              fluidRow(
                box(title="Endpoint",width=3,solidHeader = TRUE,status="primary",    

                           numericInput("Ctsno", "Number of continuous components", value=1, min=1, max=2),
                           numericInput("Binno", "Number of binary components", value=1, min=0, max=1),
                           numericInput("dichY1",label=HTML("Y<sub>1</sub> responder threshold"),value=0),
                           conditionalPanel("input.Ctsno==2",
                                            numericInput("dichY2", label=HTML("Y<sub>2</sub> responder threshold"), value=0)),
                           actionButton("GetModelComp", "Get Model")
                    ),
                
                box(title="Model Summary",width=8,solidHeader = TRUE,status="primary",
                    uiOutput("markdownComp")
                ),
                
                box(title="Parameter Estimates",width=11,solidHeader = TRUE,status="primary",
                   
                    column(width=4,
                     # Input: Select a file ----
                    fileInput("file1", "Choose CSV File",
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    # Horizontal line ----
                    tags$hr(),
                    
                    # Input: Checkbox if file has header ----
                    checkboxInput("header", "Header", TRUE),
                    
                    # Input: Select separator ----
                    radioButtons("sep", "Separator",
                                 choices = c(Comma = ",",
                                             Semicolon = ";",
                                             Tab = "\t"),
                                 selected = ","),

                     actionButton("GetParams","Obtain estimates")),
                
                column(width=7,

                   uiOutput("ResultsTable"),
                   uiOutput("LatEstTable")
                    
                   )),
                
                box(title="Sample Size Estimation",width=11,solidHeader = TRUE,status="primary",    
                   
                    conditionalPanel("output.Analysis",
                     column(width=4,
                           selectInput("alpha", "One-Sided Significance Level",c("Alpha = 0.01", "Alpha = 0.05", "Alpha = 0.10"),selected="Alpha = 0.05"),
                           #numericInput("sigmaval","Sigma Value",value=0.01,min=0,max=1000),
                           #numericInput("deltaval","Delta Value",value=0.08,min=-2,max=2),
                           sliderInput("target", "Power Target", value = .8, min = 0, max = 1),
                           numericInput("maxn", "Maximum Number of Subjects", value = 400, min = 0, max = 100000),
                           infoBoxOutput("samplesizecomp",width=11)
                     )),
                    
                    conditionalPanel("output.Analysis",
                    column(width=7,
                           plotOutput("powercomp")
                          ))
                )
                    )
            
              
                    )
              ))
)





# Define server logic required to draw a histogram
server <- function(input, output) {
  
  InputData <- eventReactive(input$file1,{
    read.csv(input$file1$datapath)
  })
  
  DataInf <- reactive(
    rawData <- InputData()
  )

  observe({
    if(input$piC2a>input$piT2a)
      showModal(modalDialog(
        title = "Warning",
        "A positive treatment effect is required in all outcomes for co-primary endpoints",
        easyClose = TRUE,
        footer = NULL
      ))
  })
  
  observe({
    if(input$piC3>input$piT3)
      showModal(modalDialog(
        title = "Warning",
        "A positive treatment effect is required in all outcomes for co-primary endpoints",
        easyClose = TRUE,
        footer = NULL
      ))
  })
  
  
  GenAnalysis <- eventReactive(input$GetParams,{
    
    if(input$Ctsno==2 && input$Binno==1){
      
      source('LatVarEst_21.R', local=TRUE)
      
      Analysis<-LatVarfunc(DataInf(),c(input$dichY1,input$dichY2))
      
    }
    
    else if(input$Ctsno==2 && input$Binno==0){
      source('LatVarEst_20.R', local=TRUE)
      Analysis<-LatVarfunc(DataInf(),c(input$dichY1,input$dichY2))
    }
    
    else if(input$Ctsno==1 && input$Binno==1){
      source('LatVarEst_11.R', local=TRUE)
      
      Analysis<-LatVarfunc(DataInf(),input$dichY1)
    }
    
    else{
      source('LatVarEst_10.R', local=TRUE)
      
      Analysis<-LatVarfunc(DataInf(),input$dichY1)
    }
  })
  
  
  
  output$ResultsTable <- renderUI({
    formattableOutput("Analysis")%>% withSpinner(color="#0dc5c1")
  })
  
  output$Analysis <- renderFormattable({
    
    Parameters <- c(HTML("&mu; <sub>T</sub>"),HTML("&mu;<sub>C</sub>"),HTML("&delta;"),HTML("&sigma;"))
    LatVar <- c(GenAnalysis()[c(3:4,1)],0.5*dim(DataInf())[1]*GenAnalysis()[2])
    Bin <- c(GenAnalysis()[c(7:8,5)],0.5*dim(DataInf())[1]*GenAnalysis()[6])
    
    dataresultstable <- data.frame(Parameters,LatVar,Bin)
    formattable(dataresultstable,col.names=(c("Parameters","Estimates","Binary")),digits=3)
    
  })
  
  output$LatEstTable <- renderUI({
    formattableOutput("parameterest")
  })
  
  output$parameterest<- renderFormattable({
    
    Parameters <- c(HTML("&delta;<sub>1</sub>"),HTML("&sigma;<sub>1</sub>"))
    LatVarests <- c(GenAnalysis()[c(9,10)])
    
    dataresultstable <- data.frame(Parameters,LatVarests)
    formattable(dataresultstable,col.names=(c("Parameters","LatVarests")),digits=3)
  })
  
  
  powercalc <- reactive({
     mean <- GenAnalysis()[1]
     var <- (ceiling(0.5*dim(DataInf())[1])*GenAnalysis()[2])
     maxn <- input$maxn
     alpha <- switch(input$alpha,
                     "Alpha = 0.01" = 0.01, 
                     "Alpha = 0.05" = 0.05, 
                     "Alpha = 0.10" = 0.10)
     Ns<-as.matrix(seq(0,maxn,1))
     betas <- apply(X=Ns,MARGIN = 2,FUN = powerfunc, 
                          mean=mean, var=var,alpha=alpha)
     return(betas)
     
  })
  
  powercalcbin <- reactive({
    mean <- GenAnalysis()[1]
    var <- (ceiling(0.5*dim(DataInf())[1])*GenAnalysis()[6])
    maxn <- input$maxn
    alpha <- switch(input$alpha,
                    "Alpha = 0.01" = 0.01, 
                    "Alpha = 0.05" = 0.05, 
                    "Alpha = 0.10" = 0.10)
    Ns<-as.matrix(seq(0,maxn,1))
    betas <- apply(X=Ns,MARGIN = 2,FUN = powerfunc, 
                   mean=mean, var=var,alpha=alpha)
    return(betas)
    
  })
  

  coprimpower <- reactive({
    
    maxn <- input$maxncop
    Ns<-as.matrix(seq(0,maxn,1))
    
    alphaco <- switch(input$alphacop,
                      "Alpha = 0.01" = 0.01, 
                      "Alpha = 0.05" = 0.05, 
                      "Alpha = 0.10" = 0.10)
    
    zalph = qnorm(1-alphaco)
    
    if(input$Ctsnoco==1 && input$Binnoco==0){
      Sigma<-input$SigmaY1
      z <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1,sigma=input$SigmaY1)
      betacoprim <- apply(X=as.matrix(z),MARGIN=2,pnorm,mean=0,sd=Sigma)
    }
    
    if(input$Ctsnoco==2 && input$Binnoco==0){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1,sigma=input$SigmaY1)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY2,sigma=input$SigmaY2)
      Sigma <- matrix(c(input$SigmaY1^2, input$rho12*input$SigmaY1*input$SigmaY2,
                      input$rho12*input$SigmaY1*input$SigmaY2, input$SigmaY2^2) ,nrow=2,ncol=2)
      betacoprim <- apply(X=as.matrix(cbind(z1,z2)),MARGIN=1,pmvnorm,lower=c(-Inf,-Inf),mean=c(0,0),sigma=Sigma)
    }
    
    if(input$Ctsnoco==1 && input$Binnoco==1){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1,sigma=input$SigmaY1)
      MeanY2a<-qnorm(input$piT2a)-qnorm(input$piC2a)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=MeanY2a,sigma=1)
      Sigma <-  matrix(c(input$SigmaY1^2, input$rho12a*input$SigmaY1, input$rho12a*input$SigmaY1, 
                       1),nrow=2,ncol=2)
      betacoprim <- apply(X=as.matrix(cbind(z1,z2)),MARGIN=1,pmvnorm,lower=c(-Inf,-Inf),mean=c(0,0),sigma=Sigma)
    }
    
    if(input$Ctsnoco==2 && input$Binnoco==1){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1,sigma=input$SigmaY1)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY2b,sigma=input$SigmaY2a)
      MeanY3<-qnorm(input$piT3)-qnorm(input$piC3)
      z3 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=MeanY3,sigma=1)
      Sigma <-  matrix(c(input$SigmaY1^2, input$rho12b*input$SigmaY1*input$SigmaY2a, input$rho13*input$SigmaY1,
                       input$rho12b*input$SigmaY1*input$SigmaY2a, input$SigmaY2a^2, input$rho23*input$SigmaY2a,
                       input$rho13*input$SigmaY1, input$rho23*input$SigmaY2a, 1), nrow=3,ncol=3)
      betacoprim <- apply(X=as.matrix(cbind(z1,z2,z3)),MARGIN=1,pmvnorm,lower=c(-Inf,-Inf,-Inf),mean=c(0,0,0),sigma=Sigma)
      
    }
    
    return(betacoprim)
    })
  
  multpower <- reactive({
    
    maxn <- input$maxnmult
    Ns<-as.matrix(seq(0,maxn,1))
    
    alphaco <- switch(input$alphamult,
                      "Alpha = 0.01" = 0.01, 
                      "Alpha = 0.05" = 0.05, 
                      "Alpha = 0.10" = 0.10)
    
     zalph = qnorm(1-alphaco)
    
    if(input$Ctsnomult==1 && input$Binnomult==0){
      Sigma<-input$SigmaY1mult
      z <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1mult,sigma=input$SigmaY1mult)
      betamult <- apply(X=as.matrix(z),MARGIN=2,pnorm,mean=0,sd=Sigma)
    }
    
    if(input$Ctsnomult==2 && input$Binnomult==0){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1mult,sigma=input$SigmaY1mult)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY2mult,sigma=input$SigmaY2mult)
      Sigma <- matrix(c(input$SigmaY1mult^2, input$rho12mult*input$SigmaY1mult*input$SigmaY2mult,
                        input$rho12mult*input$SigmaY1mult*input$SigmaY2mult, input$SigmaY2mult^2) ,nrow=2,ncol=2)
      betamult <- apply(X=as.matrix(cbind(z1,z2)),MARGIN=1,multfunc2,sigmat=Sigma)
    }
    
    if(input$Ctsnomult==1 && input$Binnomult==1){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1mult,sigma=input$SigmaY1mult)
      MeanY2multa<-qnorm(input$piT2multa)-qnorm(input$piC2multa)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=MeanY2multa,sigma=1)
      Sigma <-  matrix(c(input$SigmaY1mult^2, input$rho12multa*input$SigmaY1mult, input$rho12multa*input$SigmaY1mult, 
                         1),nrow=2,ncol=2)
      betamult <- apply(X=as.matrix(cbind(z1,z2)),MARGIN=1,multfunc2,sigmat=Sigma)
    }
    
    if(input$Ctsnomult==2 && input$Binnomult==1){
      z1 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY1mult,sigma=input$SigmaY1mult)
      z2 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=input$MeanY2multb,sigma=input$SigmaY2multa)
      MeanY3mult<-qnorm(input$piT3mult)-qnorm(input$piC3mult)
      z3 <- apply(X=Ns,MARGIN=2,FUN=zfunc,zalph=zalph,mean=MeanY3mult,sigma=1)
      Sigma <-  matrix(c(input$SigmaY1mult^2, input$rho12multb*input$SigmaY1mult*input$SigmaY2multa, input$rho13mult*input$SigmaY1mult,
                         input$rho12multb*input$SigmaY1mult*input$SigmaY2multa, input$SigmaY2multa^2, input$rho23mult*input$SigmaY2multa,
                         input$rho13mult*input$SigmaY1mult, input$rho23mult*input$SigmaY2multa, 1), nrow=3,ncol=3)
      betamult <- apply(X=as.matrix(cbind(z1,z2,z3)),MARGIN=1,multfunc3,sigmat=Sigma)
      
    }
    
    return(betamult)
  })
  
  output$powercomp <- renderPlot({
   
      Ns <- seq(0,input$maxn,1)
      Ntot<-c(Ns,Ns)
      betares <- as.numeric(powercalc())
      betabin <- powercalcbin()
      beta <- c(betares,betabin)
      Method<-c(rep("Latent",length(Ns)),rep("Binary",length(Ns)))
      
      datagg<-data.frame(Ntot,beta,Method)
      #datanew <- data.frame(Ns,betares)
      ggplot(datagg,aes(x=Ntot,y=beta,group=Method))+geom_line(aes(color=Method),size=1)+theme_bw(base_size=16)+
        geom_hline(yintercept=input$target, linetype="dashed", color = "black")+
        xlab("Number of Subjects")+ylab("Power")+theme(legend.position = "right")+ scale_color_manual(values=c("slategray3","steelblue4"))
    
      #ggplot(datanew,aes(x=Ns,y=betares))+geom_line()+theme_bw(base_size=16)+
      #    geom_hline(yintercept=input$target, linetype="dashed", color = "black")+
      #    xlab("Number of Subjects")+ylab("Power")
      
      })
  
  
    output$powercop <- renderPlot({
    
    Ns <- seq(0,input$maxncop,1)
    betares <- coprimpower()
    
    datagg<-data.frame(Ns,betares)
    
     ggplot(datagg,aes(x=Ns,y=betares))+geom_line(size=1)+theme_bw(base_size=16)+
      geom_hline(yintercept=input$targetcop, linetype="dashed", color = "black")+
      xlab("Number of Subjects")+ylab("Power")+ scale_color_manual(values="steelblue4")
    
    })
    
    
    output$powermult <- renderPlot({
      
      Ns <- seq(0,input$maxnmult,1)
      betares <- multpower()
      
      datagg<-data.frame(Ns,betares)
      
      ggplot(datagg,aes(x=Ns,y=betares))+geom_line(size=1)+theme_bw(base_size=16)+
        geom_hline(yintercept=input$targetmult, linetype="dashed", color = "black")+
        xlab("Number of Subjects")+ylab("Power")+ scale_color_manual(values="steelblue4")
      
    })
    
  
  
  output$plotcomp <- renderUI(
    plotOutput("powercomp",width="110%",height="350px")
  )
  
  output$plotcoprim <- renderUI(
    plotOutput("powercop",width="100%",height="350px")
  )
  
  output$plotmult <- renderUI(
    plotOutput("powermult",width="100%",height="350px")
  )
  
  
  output$samplesizecomp <- renderInfoBox({
    
    infoBox(
     "Individuals per arm", nrequired(), icon = icon("users", lib = "font-awesome"),
      color = "light-blue",fill=TRUE
    )
  }) 
  
 nrequired<- reactive({
   
    Ns <- seq(0,input$maxn,1)
    nrequired <- Ns[which.max(powercalc()>=input$target)]
  
    return(nrequired)
  })
 
 
 output$samplesizecopr <- renderInfoBox({
   infoBox(
     "Individuals per arm", nrequiredcopr(), icon = icon("users", lib = "font-awesome"),
     color = "light-blue", fill =TRUE
   )
 }) 
   
   nrequiredcopr <- reactive({
     Ns <- seq(0,input$maxncop,1)
     nrequiredcopr <- Ns[which.max(coprimpower()>=input$targetcop)]
     
     return(nrequiredcopr)
   })
  
   output$samplesizemult <- renderInfoBox({
     infoBox(
       "Individuals per arm", nrequiredmult(), icon = icon("users", lib = "font-awesome"),
       color = "light-blue", fill =TRUE
     )
   }) 
     
     nrequiredmult <- reactive({
     Ns <- seq(0,input$maxnmult,1)
     nrequiredmult <- Ns[which.max(multpower()>=input$targetmult)]
     
     return(nrequiredmult)
   })
     
   getPage <- function() {
     
     if(input$Ctsnoco==2 && input$Binnoco==1){
       getPage <- withMathJax(includeHTML("twoctsonebin.html")) }
     else if(input$Ctsnoco==1 && input$Binnoco==1){
       getPage <- withMathJax(includeHTML("onectsonebin.html"))}
     else if(input$Ctsnoco==2 && input$Binnoco==0){
       getPage <- withMathJax(includeHTML("twocts.html"))}
     else { getPage <- withMathJax(includeHTML("onects.html"))}
     
     return(getPage)
   }
   
   ModelPage <- eventReactive(input$GetModel,getPage())
   
   output$markdown <- renderUI({
     ModelPage()
   })
   
   getPageMult <- function() {
     
     if(input$Ctsnomult==2 && input$Binnomult==1){
       getPageMult <- withMathJax(includeHTML("TwoCtsOneBinMult.html")) }
     else if(input$Ctsnomult==1 && input$Binnomult==1){
       getPageMult <- withMathJax(includeHTML("OneCtsOneBinMult.html"))}
     else if(input$Ctsnomult==2 && input$Binnomult==0){
       getPageMult <- withMathJax(includeHTML("TwoCtsMult.html"))}
     else { getPageMult <- withMathJax(includeHTML("OneCtsMult.html"))}
     
     return(getPageMult)
   }
   
   ModelPageMult <- eventReactive(input$GetModelMult,getPageMult())
   
   output$markdownMult <- renderUI({
     ModelPageMult()
   })
   
   getPageComp <- function() {
     
     if(input$Ctsno==2 && input$Binno==1){
       getPageComp <- withMathJax(includeHTML("TwoCtsOneBinSamp.html")) }
     else if(input$Ctsno==1 && input$Binno==1){
       getPageComp <- withMathJax(includeHTML("OneCtsOneBinSamp.html"))}
     else if(input$Ctsno==2 && input$Binno==0){
       getPageComp <- withMathJax(includeHTML("TwoCtsSamp.html"))}
     else { getPageComp <- withMathJax(includeHTML("OneCtsSamp.html"))}
     
     return(getPageComp)
   }
   
   ModelPageComp <- eventReactive(input$GetModelComp,getPageComp())
   
   output$markdownComp <- renderUI({
     ModelPageComp()
   })
   
  
  }


# Run the application 
shinyApp(ui = ui, server = server)

