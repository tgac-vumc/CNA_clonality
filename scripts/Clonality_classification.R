#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clonality_classification.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Classify clonality using CNA and NGS mutations
#
# Authors: Barbara Andrade Barbosa
# Edited and compliled by Jurriaan Janssen (j.janssen4@amsterdamumc.nl) 
#
# TODO:
# 1) 
#
# History:
#  13-08-2024: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
if(!'ggsankey' %in% installed.packages()){remotes::install_github("davidsjoberg/ggsankey")}
library(ggplot2)
library(dplyr)
library(ggsankey)

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input <- snakemake@input[["Clonality_statistics"]]
    input_GMM <- snakemake@input[["GMM_model"]]
    input_SampleTable <- snakemake@input[["SampleTable"]]
    dataset <- snakemake@wildcards[["dataset"]]
    output <-  snakemake@output[["Clonality_classification"]]
}else{
    input_CNA <- 'output/CNA_Clonality_statistics_AUMC.txt'
    input_Mutation <- 'output/Mutations_Clonality_statistics_AUMC.txt'
    input_SampleTable <- 'reference/SampleTables/SampleTable_AUMC.csv'
    dataset <- 'AUMC'
    output <- 'output/Clonality_classification_AUMC.pdf'
    output_Sankey <- 'output/SankeyPlot_AUMC.pdf'
}


#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
CNA_Clonality_statistics <- read.delim(input_CNA, stringsAsFactors = F)
Mutation_Clonality_statistics <- read.delim(input_Mutation, stringsAsFactors = F)
SampleTable <- read.csv(input_SampleTable, stringsAsFactors = F)
#-------------------------------------------------------------------------------
# 2.1 Classify clonality
#-------------------------------------------------------------------------------
Clonality_classification <-
    CNA_Clonality_statistics %>%
    tidyr::separate(patient,sep = '-', into=c('Sample1','Sample2')) %>% 
    mutate(Sample1 = gsub('SU_T1.|SU_T2.|SU_','',Sample1),Sample2 = gsub('SU_T1.|SU_T2.|SU_','',Sample2)) %>%
    left_join(Mutation_Clonality_statistics) %>% 
    # Fetch two metric classifcation
    mutate(
        llr2.selected = as.numeric(gsub(',','.',llr2.selected)),
        cor = as.numeric(gsub(',','.',cor)),
        TwoMetric_classification = dplyr::case_when(
                                                 cor > 0.54 & llr2.selected > 0 ~ 'Clonal',
                                                 llr2.selected < -5 | cor < 0.45 ~  'Non-Clonal',
                                                 TRUE ~ 'Inconclusive'),
           Clonality_MC = factor(Clonality_MC,levels = c('Clonal','Non-Clonal','Probably Non-Clonal','Inconclusive')))

if(dataset == 'TRACERx'){
    # For TRACERx create subset to amount of intratumoral pairs (n=41)
    set.seed(123)
    ClonalSubsampled <- Clonality_classification %>%
        filter(True_clonality == 'Clonal') %>%
        mutate(patient = substr(Sample1,1,8)) %>%
        group_by(patient) %>%
        # Randomly slice one within group row
        slice_sample(n=1) %>%
        ungroup()

        # Create the same amount of intrapatient non-clonal pairs (n=41)

    set.seed(123)
    NonClonalSubsampled <-
        Clonality_classification %>%
        filter(True_clonality == 'Non-Clonal') %>%
        slice_sample(n=nrow(ClonalSubsampled),replace = T)
    Clonality_classification <- bind_rows(ClonalSubsampled,NonClonalSubsampled)
}


#-------------------------------------------------------------------------------
# 3.1 Plot and save data to files
#-------------------------------------------------------------------------------
if(dataset == 'TRACERx'){
    pdf(output, height = 5, width = 8)
    Clonality_classification %>%
        mutate(Clonality_MC = as.character(Clonality_MC),Clonality_MC = ifelse(Clonality_MC == 'Probably Non-Clonal','test',Clonality_MC)) %>%
        ggplot( aes(x=llr2.selected, y=cor)) +
        geom_rect(aes(xmin = -Inf,    xmax =Inf, ymin =-Inf , ymax =Inf), alpha=0.25,  fill = "lightyellow") +
        geom_rect(aes(xmin = -5, xmax = Inf,   ymin = 0.45,    ymax = Inf), alpha = 0.25, fill = "lightgrey")+
        geom_rect(aes(xmin = 0,    xmax = Inf, ymin = 0.54, ymax = Inf), alpha = 0.25, fill = "lightblue") +
        geom_point(
            aes(color=Clonality_MC,shape = True_clonality),
            alpha = 2,size= 4) +
        theme_bw(base_size = 18) +
        xlim(-20, 120) + ylim(-0.25,1) + 
        geom_hline(yintercept = 0.54, linetype = "dashed", alpha = 0.5)+
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values= c("Clonal" = "red", "Non-Clonal" = "forestgreen",'test' = '#7ba04dff' ,'Inconclusive' = 'grey38')) +
        scale_shape_manual(values=c(16,17)) +
        labs(y = 'Pearson Correlation', x = 'Log-likelihood ratio', shape = 'True clonality',color = 'NGS clonality')
    dev.off()

    pdf(output_Sankey, height = 5, width = 6)
    Clonality_classification %>%
                mutate(Clonality_MC = factor(Clonality_MC, levels = c('Non-Clonal','Probably Non-Clonal','Inconclusive','Clonal'))) %>% 
        make_long(Clonality_MC,True_clonality,TwoMetric_classification) %>%
        ggplot(aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node, levels = c('Non-Clonal','Probably Non-Clonal','Inconclusive','Clonal')),
                   label = node)) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) +
        theme_sankey(base_size = 16)+
        theme(legend.position="bottom")
    dev.off()

}else if(dataset == 'AUMC'){
    pdf(output, height = 5, width = 8)
    Clonality_classification %>%
        mutate(Clonality_MC = as.character(Clonality_MC),Clonality_MC = ifelse(Clonality_MC == 'Probably Non-Clonal','test',Clonality_MC)) %>%

        filter(!is.na(True_clonality)) %>%
        ggplot( aes(x=llr2.selected, y=cor)) +
        geom_rect(aes(xmin = -Inf,    xmax =Inf, ymin =-Inf , ymax =Inf), alpha=0.25,  fill = "lightyellow") +
        geom_rect(aes(xmin = -5, xmax = Inf,   ymin = 0.45,    ymax = Inf), alpha = 0.25, fill = "lightgrey")+
        geom_rect(aes(xmin = 0,    xmax = Inf, ymin = 0.54, ymax = Inf), alpha = 0.25, fill = "lightblue") +
        geom_point(
            aes(color=Clonality_MC,shape = True_clonality),
            alpha = 2,size= 4) +
        theme_bw(base_size = 18) +
        xlim(-20, 120) + ylim(-0.25,1) + 
        geom_hline(yintercept = 0.54, linetype = "dashed", alpha = 0.5)+
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values= c("Clonal" = "red", "Non-Clonal" = "forestgreen",'test' = '#7ba04dff' ,'Inconclusive' = 'grey38')) +
        scale_shape_manual(values=c(16,17)) +
        geom_point(data = Clonality_classification %>% filter(is.na(True_clonality)),aes(llr2.selected, y=cor),alpha=0.5,size=1,shape = 8,color =  '#F28E2B') +
        labs(y = 'Pearson Correlation', x = 'Log-likelihood ratio', shape = 'True clonality',color = 'NGS clonality')
    dev.off()

    pdf(output_Sankey, height = 5, width = 6)
    Clonality_classification %>%
        filter(!is.na(True_clonality)) %>%
        make_long(Clonality_MC,True_clonality,TwoMetric_classification) %>%
        ggplot(aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node),
                   label = node)) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) +
        theme_sankey(base_size = 16)+
        theme(legend.position="bottom")
    dev.off()
    
}







