# Script Author: JaeJin Choi, OCT-2020
# Description: Shiny app that visualize FFP distance matrices to determine optimal point(s)


require(shiny)

# required for actual work
require(ggplot2)

#require(reshape)
require(dplyr)
require(plotrix) #visualize two y-axis plot

#source(knitr::purl("visual_func_src.R", quiet=TRUE))
source("./visual_func_src.R") #load custome function

# R widget gallery
# https://shiny.rstudio.com/gallery/widget-gallery.html

ui <- fluidPage(

    # Application title
    titlePanel("Feature Frequency Profile optimum determine"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # sliderInput(inputId = "box_ratio",
            #             label = "Threshold for choosing boxes (bag1 versus bag2)",
            #             min = 0,
            #             max = 1,
            #             value = 0.5
            #             ),
            # 
            # sliderInput(inputId = "n_repeats",
            #             label = "Simulation repeats",
            #             min = 5,
            #             max = 1000,
            #             value = 100
            #             ),
            # 
            # numericInput(inputId = "rnd_seed", 
            #              label = "Random seed",
            #              value = 123,
            #              step = 1
            #              ),
            # 
            # radioButtons(inputId = "ball_color"
            #              , label = ("Choose ball color")
            #              , choices = list("Red" = "red", "Blue" = "blue", "white" = "white")
            #              , selected = "blue" # select value (not the label)
            #              ),

            selectInput("select_project", label = ("Select project"), 
                        choices = list("Select project title"="none"
                                       , "NCBI CoV 66 virus sample" = "1"
                                       , "Shakespear" = "/home/jjc/Desktop/project-work/shakespeare_project/FFP_alphanumeric-play/bionj"
                                       , "Tree of Life" = "2"
                                       , "test matrices" = "/DATA/workspace/berkeley_schedule/semester_progress/fall_2020/stat_130/project/01_optimal_range_plot/matrix"), #path = title
                        #selected = "none"
                        selected = "/DATA/workspace/berkeley_schedule/semester_progress/fall_2020/stat_130/project/01_optimal_range_plot/matrix"
                        #selected = "/home/jjc/Desktop/project-work/shakespeare_project/FFP_alphanumeric-play/bionj"
                        ),
            
            fluidRow(column(3, paste0("Project data path", "\n"), verbatimTextOutput("select_project"))),
            
            textInput("name_regex", label = ("File selection by regular expression"), value = "sym.mat"),
            
            #actionButton("do_run", "Plot"),

            checkboxGroupInput("other_indicators", label = ("Check other indicators (only if applicable; affected by l-mer range)"), 
                               choices = list("JSD_upperlimit (gray)" = 1, "Maximum JSD_sd (red)" = 2, "Minimum RF (blue)" = 3),
                               selected = NA #unselected
                               ),

            #initial slider
            sliderInput("lmer_range", label = ("Limit l-mer range")
                        , min = 1
                        , max = 10
                        , value = c(1, 10)
            ),
            
            hr() #,
            #tableOutput(outputId = "proportions")

        ),
        
        # Show plot
        mainPanel(
           #tableOutput(outputId = "proportions"), 
           plotOutput("main_plot")
           
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    output$select_project <- renderPrint({input$select_project})
    
    m_combined_df <- reactive({
        if (input$select_project!="none")
        {
            #read distance matrix file, either symmetric or triangular matrix
            file_list <- list.files(path=input$select_project
                                    , full.names=T #full file path, F will get only file name
                                    , pattern = input$name_regex
            ) 
            
            #print(paste0("what project: ", input$select_project))
            #ffp_visual_optimum(load_path=input$select_project, name_regex="sym.mat.")
            m_df <- data.frame
                
            if (length(file_list)>0)
            {
                m_df <- ffp_visual_optimum(file_list=file_list, name_regex=input$name_regex)
            }

            m_df
                
        } else
        {
            NULL
        }
        
    })

    #https://shiny.rstudio.com/reference/shiny/0.14/updateSliderInput.html
    observeEvent(!is.null(m_combined_df()), {
        feature_range <- m_combined_df()$feature_length
        #print(feature_range)
        #update slider parameters
        updateSliderInput(session, inputId = "lmer_range"
                          , label = ("Limit l-mer range-updated")
                          , min = min(feature_range)
                          , max = max(feature_range)
                          , value = c(min(feature_range), max(feature_range))
        )
    })
    
    
    #outputOptions(output, "lmer_range", suspendWhenHidden = FALSE)
    
    output$main_plot <- renderPlot({
        if (!is.null(m_combined_df())) #render if not null
            {
            pm_df <- m_combined_df()
            #pm_df <- na.omit(pm_df) #not necessary
            
            # print(pm_df$feature_length)
            # lmer_from <- which(pm_df$feature_length == input$lmer_range[1], arr.ind=TRUE)
            # lmer_to <- which(pm_df$feature_length == input$lmer_range[2], arr.ind=TRUE)
 
            

            #print(pm_df[lmer_from:lmer_to, ])
            
            #pm_df <- pm_df[lmer_from:lmer_to, ]
            
            twoord.plot(lx = pm_df$feature_length
                        , ly = pm_df$tree_dist
                        , rx = pm_df$feature_length
                        , ry = pm_df$jsd_sd
                        , xlab = "Feature length (k-mer)"
                        , ylab = "Tree distance (k-mer versus k-mer+1)"
                        , rylab = "JSD-sd"
                        , main = "RF variation and JSD_sd trend"
                        , lcol = "black" #black
                        , rcol = "red" #red
                        ) #+
    
            
            ### display other_indicators (check box)
            if (1 %in% input$other_indicators & sum(pm_df$jsd_upperlimit)!=0) #show a point of JSD-upperlimit if exists
            {
                #print(pm_df$feature_length[pm_df$jsd_upperlimit!=0][1])
                abline(v = pm_df$feature_length[pm_df$jsd_upperlimit!=0][1] #the first feature length showing JSD_upperlimit
                       , col = "gray") #color red==2
            }
    
            if (2 %in% input$other_indicators) #show a point of maximum JSD_sd, in red vertical lines
            {
                #jsd_sd_max <- max(pm_df$jsd_sd)
                #print(pm_df$feature_length[pm_df$jsd_sd==jsd_sd_max])
                jsd_sd_max_lmer <- pm_df$feature_length[pm_df$jsd_sd==max(pm_df$jsd_sd)]
                abline(v = jsd_sd_max_lmer #the first feature length showing JSD_upperlimit
                        , col = "red" #color blue==3
                        , lty = 2) #long dashed line
                
                ## overlapping using + not working, "non-numeric argument to binary operator"
                # ty_plot <- ty_plot + abline(v = pm_df$feature_length[pm_df$jsd_sd==jsd_sd_max] #the first feature length showing JSD_upperlimit
                #        , col = "blue" #color blue==3
                #        , lty = 3) #long dashed line
            }
    
            if (3 %in% input$other_indicators) #show a point of minimum RF, in blue vertical lines
            {
                #jsd_sd_max <- max(pm_df$jsd_sd)
                #print(pm_df$feature_length[pm_df$jsd_sd==jsd_sd_max])
                rf_min_lmer <- pm_df$feature_length[pm_df$tree_dist==min(pm_df$tree_dist, na.rm=T)]
                print(rf_min_lmer)
                
                abline(v = rf_min_lmer #point(s) of minimal RF (e.g., RF==0)
                       , col = "blue" #color blue==3
                       , lty = 2) #fine dashed line
    
            }
            
        } else
        {
            NULL
        }
    })

    
}

# Run the application 
shinyApp(ui = ui, server = server)
