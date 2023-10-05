##Libraries
library(shiny)
library(shinycssloaders)
library(rtracklayer)
library(dplyr)
library(tidyverse)
library(gggenomes)
library(GenomicRanges)
library(cowplot)
library(paletteer)

## Functions
source("Functions.R")

##UserInterface
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
      @import url('https://fonts.googleapis.com/css2?family=Yusei+Magic&display=swap');
      h1, h2, h3, h4, h5, h6 {
      font-family: 'Yusei Magic', sans-serif;
      }
      ")
    )
  ),
  titlePanel(
    h1("WormSynteny", align = "center")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h3("Plot controls", align = "center"),
      
      helpText("Zoom in or out by 1000b on the plot (500b at both ends)", align = "center"),
      
      fluidRow(
        actionButton("ZoomIn","Zoom In",
                     icon = icon("magnifying-glass-plus", lib = "font-awesome"),
                     width = "200px",
                     style = "margin-left:60px;"
        ),
        tags$hr(style = "margin-top: 5px; margin-bottom: 5px;"),
        
      ),
      fluidRow(
        actionButton("ZoomOut","Zoom Out",
                     icon = icon("magnifying-glass-minus", lib = "font-awesome"),
                     width = "200px",
                     style = "margin-left:60px;"
        ),
      ),
      
      helpText("Select species sequences to switch them from bottom to top or vice versa", align = "center"),
      
      div(
        checkboxGroupInput("species", "Select Species:",
                           choices = c("bovis", "becei", "panamensis", "inopinata", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
        ),
        style = "margin-left: 80px;"
      ),
      
      actionButton("switch", "Switch", style = "margin-left:80px;"),
      
    ),
    mainPanel(
      h3("Coordinates selection", align = "center"),
      fluidRow(
        column(4,
               selectInput("Chr", "Choose a chromosome", choices = c("I", "II", "III", "IV", "V", "X"))
        ),
        column(4,
               numericInput("Start", "Enter start position", value = NULL)),
        column(4,
               numericInput("End", "Enter end position", value = NULL)),
      ),
      fluidRow(
        column(4,
               numericInput("Gap", "Enter a gap value", value = 20, min = 20, max = 1000)),
        column(4,
               numericInput("Filter", "Enter a filtering percentage", value = 5, min = 5, max = 90)),
        hr(),
        column(4,
               actionButton("Search","Search",icon = icon("magnifying-glass", lib = "font-awesome"), width = "200px", style = "margin-left:50px"))
      ),
      tabsetPanel(
        tabPanel(title = "Plot",
                 withSpinner(plotOutput("gggenomes_plot")))),
      fluidRow(
        column(12,
               verbatimTextOutput("explanation")))
    )
  )
)

##Server
server <- function(input, output, session){
  
  #To avoid the spinner to be triggered before Search button has been clicked 
  output$gggenomes_plot <- renderPlot(NULL)
  
  #Dependent of Search button : must be click to start plotting process
  observeEvent(input$Search, {
    
    #Use isolate() to avoid dependency on input$Start, input$End, input$Gap and input$Filter
    node <- isolate(
      data.frame(Chr = input$Chr, Start = input$Start, End = input$End, Gap = input$Gap, Filter = input$Filter)
    )
    
    #Create a plot that depends on the selected coordinates 
    output$gggenomes_plot <- renderPlot({
      #Save node for zoom_in and zoom_out
      write.csv(node,paste0(folder_ubuntu,"node.csv"))
      
      #Create the C.elegans region 
      region <- create_region(node$Chr, node$Start, node$End)
      
      #Run the pipeline 
      data <- general_function(region, node$Gap, node$Filter)
      
      #Save data for switching function
      genes <- data$genes
      write.csv(genes, paste0(folder_ubuntu,"genes.csv"))
      seqs <- data$seqs
      write.csv(seqs, paste0(folder_ubuntu,"seqs.csv"))
      links <- data$links
      write.csv(links, paste0(folder_ubuntu,"links.csv"))
      
      #Gggenomes plotting
      # If no genes and no links meaning there is no aligned regions : plot only sequences, otherwise plot everything
      if(!(is.null(data$genes)) | !(is.null(data$links))){
        
        if(any(data$genes$seq_id %in% data$links$seq_id2)){
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- node$End - node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          #Plot 
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          ## Create a second legend for C.elegans genes
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
          
        }else{
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- node$End - node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
          
        }
        
      }else{
        elegans_length <- node$End - node$Start
        
        p3 <- gggenomes(seqs = data$seqs) +
          geom_seq() +
          geom_bin_label(fontface = "italic",size = 5) +
          scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
          theme(axis.text.x=element_text(size=15))
        
        print(p3) # Return the plot
        
      }
      
    })
    
    
    output$explanation <- renderText({

      #First sentence always appearing
      string_normal <- paste0("Showing synteny for chromosome ", input$Chr, " from ", input$Start, " to ", input$End,".\n")
      
      ## Genes 
      #Load genes
      genes <- read.csv(paste0(folder_ubuntu,"genes.csv"), header = FALSE)
      l <- vector("character", length = 0)
      
      # Check if the file is empty
      if (nrow(genes) == 0) {
        string_no_genes <- "No genes found in this selected region of C.elegans.\n"
      } else {
        
        genes <- read.csv(paste0(folder_ubuntu,"genes.csv")) %>% select(-X)
        
        #Select only top genes
        genes_top <- genes %>% filter(!(str_detect(seq_id, "//(filtered)//")))
        
        #Check if all the 11 species are present
        spe <- c("bovis", "becei", "panamensis", "inopinata", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
        j <- 1
        
        for(species in spe){
          check <- genes_top %>% filter(str_detect(seq_id, species))
          if(nrow(check)> 0 ){
            l[j]<- species
            j <- j + 1
          }
        }
        
        species_with_no_genes <- setdiff(spe, l)
        if(length(species_with_no_genes)!=0){
          species_string <- paste(species_with_no_genes, collapse = ", ")
          string_no_genes <- paste0("Some of the species do not have gene in the selected region of C.elegans.\n", " Here is the list of species that do not contain genes on the plot :",species_string, ".")
        }else{
          string_no_genes <- "Every species contain at least one gene in the selected region of C.elegans."
        }
      }
      
      
      ## Alignments
      species_no_alignment <- read.csv(paste0(folder_ubuntu,"species_no_alignment.csv"), header = FALSE)
      if(nrow(species_no_alignment)!=0){
        species_no_alignment <- read.csv(paste0(folder_ubuntu,"species_no_alignment.csv"))
        no_ali <- paste(species_no_alignment$V1, collapse = ", ")
        string_no_alignment <- paste0("Some of the species do not have alignments in the considered region.\n", " Here is the list of the missing species that do not appear on the plot :",no_ali, ".")
        
      }else{
        
        string_no_alignment <- "All species have aligned regions in the selected region of C.elegans."
        
      }
      
      message <- paste0(string_normal," ", string_no_alignment, " ", string_no_genes)
      return(message)
      
    })
  })
  
  #Dependent of ZoomOut button : must be click to start plotting process
  observeEvent(input$ZoomOut, {
    
    #Load node previously saved in zoom_in/out or search sections
    node <- read.csv(paste0(folder_ubuntu,"node.csv"))
    
    #Zoom out : -500bp and +500bp 
    new_node <- node %>% mutate(Chr = node$Chr, Start = ifelse(node$Start>500,node$Start-500,node$Start), End = node$End + 500, width = (node$End-node$Start)+1)
    
    #Use isolate() to avoid dependency on input$Gap and input$Filter
    new_node <- isolate(new_node %>% mutate(Gap = input$Gap, Filter = input$Filter))
    
    #Create a plot that depends on the previous plot and update the settings (filter and gap values)
    output$gggenomes_plot <- renderPlot({
      
      #Create the C.elegans region 
      region <- create_region(new_node$Chr, new_node$Start, new_node$End)
      
      #Run the pipeline 
      data <- general_function(region, new_node$Gap, new_node$Filter)
      
      #Save data for switching function
      genes <- data$genes
      write.csv(genes,paste0(folder_ubuntu, "genes.csv"))
      seqs <- data$seqs
      write.csv(seqs, paste0(folder_ubuntu,"seqs.csv"))
      links <- data$links
      write.csv(links, paste0(folder_ubuntu,"links.csv"))
      
      #Gggenomes plotting
      # If no genes and no links meaning there is no aligned regions : plot only sequences, otherwise plot everything
      if(!(is.null(data$genes)) | !(is.null(data$links))){
        if(any(data$genes$seq_id %in% data$links$seq_id2)){
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- new_node$End - new_node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          #Plot 
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          ## Create a second legend for C.elegans genes
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
          
          
        }else{
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- new_node$End - new_node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
        }
        
      }else{
        elegans_length <- new_node$End - new_node$Start
        
        p3 <- gggenomes(seqs = data$seqs) +
          geom_seq() +
          geom_bin_label(fontface = "italic",size = 5) +
          scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
          theme(axis.text.x=element_text(size=15))
        
        print(p3) # Return the plot
        
      }
      
    })
    
    #Dependent of ZoomOut button : must be click to display the message 
    output$explanation <- renderText({
      
      zoom <- paste0("Start : ", node$Start, " became ",new_node$Start," and End : ",node$End, " became ",new_node$End)
      string <- paste0("Zooming out.\n",zoom)
      
      string
    })
    
    #Save node for zoom_in and zoom_out
    write.csv(new_node, paste0(folder_ubuntu,"node.csv"))
    
  })
  
  #Dependent of ZoomIn button : must be click to start plotting process
  observeEvent(input$ZoomIn, {
    
    #Load node previously saved in zoom_in/out or search sections
    node <- read.csv(paste0(folder_ubuntu,"node.csv"))
    
    #Zoom In : -500bp and +500bp 
    new_node <- node %>% mutate(Chr = node$Chr, Start2 = ifelse((node$End-500)-(node$Start+500)<=0, node$Start, node$Start+500), End =  ifelse((node$End-500)-(node$Start+500)<=0, node$End, node$End-500), width = (node$End-node$Start)+1) %>% select(-Start) %>% dplyr::rename(Start = "Start2")
    
    #Use isolate() to avoid dependency on input$Gap and input$Filter
    new_node <- isolate(new_node %>% mutate(Gap = input$Gap, Filter = input$Filter))
    
    #Create a plot that depends on the previous plot and update the settings (filter and gap values)
    output$gggenomes_plot <- renderPlot({
      
      
      #Create the C.elegans region 
      region <- create_region(new_node$Chr, new_node$Start, new_node$End)
      
      #Run the pipeline 
      data <- general_function(region, new_node$Gap, new_node$Filter)
      
      #Save data for switching function
      genes <- data$genes
      write.csv(genes, paste0(folder_ubuntu,"genes.csv"))
      seqs <- data$seqs
      write.csv(seqs, paste0(folder_ubuntu,"seqs.csv"))
      links <- data$links
      write.csv(links, paste0(folder_ubuntu,"links.csv"))
      
      #Gggenomes plotting
      # If no genes and no links meaning there is no aligned regions : plot only sequences, otherwise plot everything
      if(!(is.null(data$genes)) | !(is.null(data$links))){
        if(any(data$genes$seq_id %in% data$links$seq_id2)){
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- new_node$End - new_node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          #Plot 
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          ## Create a second legend for C.elegans genes
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
          
        }else{
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- new_node$End - new_node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
        }
        
      }else{
        elegans_length <- new_node$End - new_node$Start
        
        p3 <- gggenomes(seqs = data$seqs) +
          geom_seq() +
          geom_bin_label(fontface = "italic",size = 5) +
          scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(new_node$Start, new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (new_node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-new_node$Start)/5))) +
          theme(axis.text.x=element_text(size=15))
        
        print(p3) # Return the plot
        
      }
      
    })
    
    #Dependent of ZoomIn button : must be click to display the message
    output$explanation <- renderText({
      
      zoom <- paste0("Start : ", node$Start, " became " , new_node$Start," and End : ",node$End, " became ",new_node$End)
      string <- paste0("Zooming in.\n",zoom)
      
      string
    })
    #Save node for zoom_in and zoom_out
    write.csv(new_node, paste0(folder_ubuntu,"node.csv"))
  })
  
  #Dependent of switch button : must be click to start plotting process
  observeEvent(input$switch, {
    
    #Initialize a list to store the selected species (checked boxes)
    selected_species <- list()
    
    #Get the currently selected species
    selected_species <- isolate(input$species)
    
    #Create a plot that depends on the previous plot and update the settings (switch sequences)
    output$gggenomes_plot <- renderPlot({
      
      #Load the previously set of genes created by either search, zoom_in, or zoom_out buttons
      data <- list()
      data$genes <- read.csv(paste0(folder_ubuntu,"genes.csv"), header = FALSE)
      data$seqs <- read.csv(paste0(folder_ubuntu,"seqs.csv")) %>% select(-X)
      
      if(nrow(data$genes)==0){
        #Run the function switch 
        new_genes <- switch_filtered_genes(selected_species, data$seqs, data$genes)
        data$links <- read.csv(paste0(folder_ubuntu,"links.csv"), header = FALSE)
        
      }else{
        data$genes <- read.csv(paste0(folder_ubuntu,"genes.csv")) %>% select(-X)
        
        #Run the function switch 
        new_genes <- switch_filtered_genes(selected_species, data$seqs, data$genes)
        
        data$links <- read.csv(paste0(folder_ubuntu,"links.csv"))
      }
      
      
      #Run the function seqs 
      new_seqs <- update_seqs(selected_species, data$seqs, new_genes)
      
      #Store the up-to-date genes and seqs set (switched bottom/top sequences)
      write.csv(new_genes,paste0(folder_ubuntu, "genes.csv"))
      write.csv(new_seqs, paste0(folder_ubuntu,"seqs.csv"))
      
      #Update data list with new genes and new seqs sets 
      data$genes <- new_genes
      
      data$seqs <- new_seqs
      
      #Gggenomes plotting
      # If no genes and no links meaning there is no aligned regions : plot only sequences, otherwise plot everything
      if(!(nrow(data$genes)==0) | !(nrow(data$links)==0)){
        if(any(data$genes$seq_id %in% data$links$seq_id2)){
          
          node <- read.csv(paste0(folder_ubuntu,"./node.csv"))
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- node$End - node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          #Plot 
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          ## Create a second legend for C.elegans genes
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes, links = data$links) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            geom_link(offset = 0.25) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
          
          
        }else{
          node <- read.csv(paste0(folder_ubuntu,"./node.csv"))
          
          color_palette <- paletteer_d("ggthemes::Tableau_20")
          color_palette <- paletteer_d("ggthemes::Tableau_20")[-c(13,14)]
          data$genes$color <- color_palette[rank(data$genes$Orthogroup, na.last = "keep")]
          
          #Create the elegans gene_id data frame to link gene_id, orthogroup and color information
          elegans_length <- node$End - node$Start
          elegans_genes <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- subset_by_overlap(elegans, elegans_genes) %>% select(gene_id, Orthogroup)
          elegans_orthogroups <- data$genes %>% filter(str_detect(seq_id, "elegans"))
          elegans_genes <- merge(elegans_genes, elegans_orthogroups, by = "Orthogroup")
          
          p1 <- gggenomes(seqs = data$seqs, genes = data$genes) +
            geom_seq() +
            geom_gene(aes(fill = ifelse(is.na(Orthogroup) & nrow(data$genes) == 1, "NA", Orthogroup)), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            theme(axis.text.x=element_text(size=15))+
            labs(fill = "Orthogroups")+
            scale_fill_manual(values = setNames(data$genes$color, data$genes$Orthogroup))+
            theme(legend.text=element_text(size=15))+
            theme(legend.title = element_text(size=15))
          
          p2 <- gggenomes(seqs = data$seqs, genes = elegans_genes) +
            geom_seq() +
            geom_gene(aes(fill = gene_id), stroke = 0.5) +
            geom_bin_label(fontface = "italic", size = 5,expand_left =0.8) +
            scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
            labs(fill = "C.elegans genes")+
            scale_fill_manual( values =setNames(elegans_genes$color, elegans_genes$gene_id), label = elegans_genes$gene_id)+
            theme(legend.text=element_text(size=15))+ # adjust size of legend labels here
            theme(legend.title = element_text(size=15))
          
          #Remove the first legend 
          p3 <- p1 + theme(legend.position = "none")
          
          # Combine the two legends using cowplot
          legend1 <- get_legend(p1)
          legend2 <- get_legend(p2)
          combined_legend <- plot_grid(legend1, legend2, ncol = 1, align = "v", rel_heights = c(1, 1))
          
          # Add the combined legend to the plot
          plot_with_legend <- plot_grid(p3, combined_legend, ncol = 2, align = "h", axis = "tb", rel_widths = c(0.8, 0.2))
          
          # Display the plot
          print(plot_with_legend)
        }
        
      }else{
        node <- read.csv(paste0(folder_ubuntu,"./node.csv"))
        elegans_length <- node$End - node$Start
        
        p3 <- gggenomes(seqs = data$seqs) +
          geom_seq() +
          geom_bin_label(fontface = "italic",size = 5) +
          scale_x_continuous(breaks = seq(0, max(data$seqs$length), max(data$seqs$length)/5), labels = round(seq(node$Start, node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0)), (node$End+(ifelse(max(data$seqs$length)-elegans_length > 0,max(data$seqs$length)-elegans_length,0))-node$Start)/5))) +
          theme(axis.text.x=element_text(size=15))
        
        print(p3) # Return the plot
        
      }
    })
  })
}

##Run the application 
shinyApp(ui = ui, server = server)

