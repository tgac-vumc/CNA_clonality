#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clonality_classification.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Classify clonality using GMM and 
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
library(mclust)
library(ggplot2)
library(dplyr)

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
    input <- 'output/Clonality_statistics_UMCG.txt'
    input_GMM <-'output/GMM_model.Rds'
    input_SampleTable <- 'reference/SampleTables/SampleTable_UMCG.csv'
    dataset <- 'UMCG'
    output <- 'output/Clonality_classification_UMCG.pdf'
}

#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
Clonality_statistics <- read.delim(input)
GMM <- readRDS(input_GMM)
SampleTable <- read.csv(input_SampleTable)
#-------------------------------------------------------------------------------
# 2.1 Classify clonality
#-------------------------------------------------------------------------------
# Match Subtype information
if(dataset == 'UMCG'){
    Clonality_statistics <-
        Clonality_statistics %>%
        mutate(patient = paste0(substring(patient,1,4),'-',substring(patient,nchar(patient)-4,nchar(patient)-1)))
    SampleTable <-  SampleTable[match(purrr::map_chr(Clonality_statistics$patient,~strsplit(.x,'-')[[1]][1]),SampleTable$sample),]
}else{
    SampleTable <- SampleTable[match(purrr::map_chr(Clonality_statistics$patient,~strsplit(.x,'_')[[1]][1]),SampleTable$sample),]
}
Subtypes <- SampleTable$subtype

# Fetch LLR for subtype
LLR <- c()
for(i in 1:nrow(Clonality_statistics)){
    if(Subtypes[i] == 'Other'){
        LLR <- c(LLR, Clonality_statistics$llr2.split[i])
    }else{
        LLR <- c(LLR, Clonality_statistics[i,paste0('llr2.',tolower(Subtypes[i]),'.split')])
    }
}

#-------------------------------------------------------------------------------
# Create GMM input
GMM_input <- data.frame(llr = LLR,cor = Clonality_statistics$cor)

# Fit data to GMM
GMM_classification <- predict.Mclust(GMM, newdata = GMM_input)
#-------------------------------------------------------------------------------
Clonality_classification <-
    Clonality_statistics %>%
    mutate(
        llr2 = LLR,
        # Fetch GMM classification
        GMM_classification = ifelse(GMM_classification$classification<3,'Clonal','Non-Clonal'),
        # Fetch two metric classifcation
        TwoMetric_classification = dplyr::case_when(
                                                             cor > 0.54 & llr2 > 0 ~ 'Clonal',
                                                             llr2 < -5 | cor < 0.45 ~  'Non-Clonal',
                                                             TRUE ~ 'Inconclusive'),
        # Fetch true clonality
        pat1 = purrr::map_chr(patient,~strsplit(strsplit(.x,'\\-')[[1]][1],'_')[[1]][1]),

        pat2 = purrr::map_chr(patient,~strsplit(strsplit(.x,'\\-')[[1]][2],'_')[[1]][1]),

        True_clonality = ifelse(pat1 == pat2, 'Clonal','Non-Clonal'))

Clonality_classification
if(dataset == 'AUMC'){
    Clonality_classification$True_clonality <- 'Unknown'
}

#-------------------------------------------------------------------------------
# 3.1 Plot data
#-------------------------------------------------------------------------------
pdf(output, height = 5, width = 8)
Clonality_classification %>%
    ggplot( aes(x=llr2, y=cor)) +
    geom_rect(aes(xmin = -Inf,    xmax =Inf, ymin =-Inf , ymax =Inf), alpha=0.25,  fill = "lightyellow") +
    geom_rect(aes(xmin = -5, xmax = Inf,   ymin = 0.45,    ymax = Inf), alpha = 0.25, fill = "lightgrey")+
    geom_rect(aes(xmin = 0,    xmax = Inf, ymin = 0.54, ymax = Inf), alpha = 0.25, fill = "lightblue") +
    geom_point(
        aes(color= GMM_classification,shape = True_clonality),
        alpha = 2,size= 4) +
    theme_bw(base_size = 18) +
    xlim(-20, 120) + ylim(-0.25,1) + 
    geom_hline(yintercept = 0.54, linetype = "dashed", alpha = 0.5)+
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values= c("Clonal" = "red", "Non-Clonal" = "forestgreen")) +
    scale_shape_manual(values=c(16,4,1)) +
    labs(y = 'Pearson Correlation', x = 'Log-likelihood ratio', shape = 'True clonality',color = 'GMM clonality', size = 'Jaccard index')
dev.off()
