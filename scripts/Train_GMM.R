#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Train_GMM.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Train Gaussian Mixture Model on TRACERx clonality statistics
#
# Authors: Barbara Andrade Barbosa
# Edited and compliled by Jurriaan Janssen (j.janssen4@amsterdamumc.nl) 
#
# TODO:
# 1) 
#
# History:
#  07-08-2024: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
library(mclust)
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input <- snakemake@input[["Clonality_statistics"]]
    input_SampleTable <- snakemake@input[["SampleTable"]]
    output <-  snakemake@output[["GMM_model"]]
    output_summary <-  snakemake@output[["GMM_summary"]]

}else{
    input <-'output/Clonality_statistics_TRACERx.txt'
    input_SampleTable <- 'reference/SampleTables/SampleTable_TRACERx.csv'
    output <-'output/GMM_model.Rds'
}

#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
Clonality_statistics <- read.delim(input)
SampleTable <- read.delim(input_SampleTable, sep = ',')

#-------------------------------------------------------------------------------
# 3.1 Train GMM
#-------------------------------------------------------------------------------
# Match Subtype information
SampleTable <- SampleTable[match(purrr::map_chr(Clonality_statistics$patient,~strsplit(.x,'_')[[1]][1]),SampleTable$sample),]
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
# Create training data
GMM_input <- data.frame(llr = LLR,cor = Clonality_statistics$cor)

# Fit model with three components
BIC <- mclustBIC(GMM_input, 3)
GMM <- Mclust(GMM_input, x=BIC)

# Plot summary
pdf(output_summary)
plot(GMM)
dev.off()
#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
saveRDS(GMM , output)
