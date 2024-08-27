#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate_Clonality_statistics.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Calculate Clonality statistics from pairs of segmented CNA profiles
#
# Authors: Barbara Andrade Barbosa / Tim Mocking
# Edited and compliled by Jurriaan Janssen (j.janssen4@amsterdamumc.nl) 
#
# TODO:
# 1) 
#
# History:
#  31-07-2024: File creation, compile clonality code
#  07-08-2024: Create dataset-specfic comparisons
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
if(!'QDNAseq.dev' %in% installed.packages()){devtools::install_github("tgac-vumc/QDNAseq.dev", ref="clonality")}
library(QDNAseq.dev)
library(Biobase)
library(DNAcopy)
library(Clonality)
source('scripts/clonality_analysis_PFREQ_MOD.R')
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------

if(exists("snakemake")){
    input_segments <- snakemake@input[["Segments"]]
    input_SampleTable <- snakemake@input[["SampleTable"]]
    dataset <- snakemake@wildcards[["dataset"]]
    output <-  snakemake@output[["Clonality_statistics"]]
}else{
    input_segments <-'data/SupplementaryData_TRACERx.Rds'
    input_SampleTable <- 'reference/SampleTables/SampleTable_TRACERx.csv'
    dataset <- 'TRACERx'
    output <-'output/Clonality_statistics_TRACERx.txt'
}

#-------------------------------------------------------------------------------
# 1.2 Specify functions
#-------------------------------------------------------------------------------
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
    # Taken from https://stackoverflow.com/a/17171754
    x <- unique(x)

    y <- unique(y)

    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])

        if(length(z)) cbind(x[i], z, deparse.level=0)
    }

    do.call(rbind, lapply(seq_along(x), g))
}


get_ref_data <- function(path, split) {
  refData <- read.table(path, skip=0, fill=TRUE, header=TRUE, sep="\t",quote="") 
  duplicated(refData[,3]) -> exc
  CNA(refData[!exc,c(-1,-2,-3,-4)], refData[!exc,2], refData[!exc,4] + refData[!exc,3] / 2) -> refCNA
  refCNA = refCNA[!grepl("X|Y|23", refCNA$chrom),]
  if (split == TRUE){
    refCNA$maploc = refCNA$maploc/1000
    refCNA$chrom = as.numeric(refCNA$chrom)
    arms1 = splitChromosomes(refCNA$chrom, refCNA$maploc)
    refCNA$chrom = arms1
    refCNA$maploc = refCNA$maploc*1000
  }
  return(refCNA)
}

get_sample_data <- function(data, split){
  CNA(assayDataElement(data, "copynumber"), chromosomes(data), bpstart(data)) -> dataCNA2
  dataCNA2 = dataCNA2[!grepl("X|Y", dataCNA2$chrom),]
  dataCNA = dataCNA2
  if (split == TRUE){
    dataCNA$maploc = dataCNA$maploc/1000
    dataCNA$chrom = as.numeric(dataCNA$chrom)
    arms = splitChromosomes(dataCNA$chrom, dataCNA$maploc)
    dataCNA$chrom = arms
    dataCNA$maploc = dataCNA$maploc*1000
  }
  return(dataCNA)
}

add_chr <- function(column){
  new_values <- c()
  for (i in column){
    if (i < 10){
      new_values <- append(new_values, paste('chr0', i, sep=''))
    }
    else {
      new_values <- append(new_values, paste('chr', i, sep=''))
    }
  }
  return(new_values)
}

# LLR references
LLR2_ref <- function(data,refData, split){
  colnames(refData) <- c("chrom", "pg", "pl", "pn")
  refData$chrom <- as.character(refData$chrom)
  if (split == TRUE){
    clonality.analysis.mod(sampleCNA_split, ptlist=rep(pjct, ncol(data)), pfreq=refData) -> llrData
  }
  if (split == FALSE){
    clonality.analysis.mod(sampleCNA_whole, ptlist=rep(pjct, ncol(data)), pfreq=refData) -> llrData
  }
  clonalityTest(data, sbjctLst=rep(pjct, ncol(data)), pfreq=refData, llrData=llrData) -> cln
  return(data.frame(cln$clonTab)$llr2)
}


#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
Segments <- readRDS(input_segments)
SampleTable <- read.csv(input_SampleTable)
#-------------------------------------------------------------------------------
# 2.1 Fetch sample/patient names
#-------------------------------------------------------------------------------
SampleNames <- colnames(Segments)
PatientNames <- as.vector(sapply(SampleNames,function(x)strsplit(x,'_')[[1]][1]))
nPatients <- length(unique(PatientNames))

if(dataset == 'UMCG'){
    PatientNames <- substring(PatientNames,1,nchar(PatientNames)-1)
}
#-------------------------------------------------------------------------------
# 2.1 Create Comparions
#-------------------------------------------------------------------------------
if(dataset %in% c('TRACERx','UMCG')){
    # Create empty data frames to store comparisons
    IntraPatient <- data.frame()
    InterPatient <- data.frame()
    for(Patient in unique(PatientNames)){
        # skip samples with only 1 tumor
        if(sum(Patient == PatientNames) == 1){
            Patient <- sample(PatientNames,1)
        }
        # Fetch patient indices
        PatientIx <- which(Patient == PatientNames)
        if(dataset == 'TRACERx'){
            # randomly match two tumors from same patient
            IntraPatient <- rbind(IntraPatient,matrix(sample(PatientIx,2),ncol=2))
            # randomly match two tumors from different patients
            # Select a sample with the same subtype
            Subtype <- SampleTable$subtype[SampleTable$sample == Patient]
            if(Subtype == 'Other'){
                Ix_to_consider <- 1:nPatients
            }else{
                Ix_to_consider <- which(SampleTable$subtype == Subtype)
            }
            InterPatient <- rbind(InterPatient,
                                  matrix(c(PatientIx[1],sample(subset(Ix_to_consider,Ix_to_consider != PatientIx),1)),ncol=2))

        }else if(dataset == 'UMCG'){
            # For UMCG match all intra patient tumors and generate an equal amount of interpatient comparisons
            IntraPatient <- rbind(IntraPatient,expand.grid.unique(PatientIx,PatientIx))
            InterPatient_Comparisons <- expand.grid(PatientIx,subset(seq_along(PatientNames),!seq_along(PatientNames) %in% PatientIx))
            InterPatient <- rbind(InterPatient, InterPatient_Comparisons[sample(1:nrow(InterPatient_Comparisons),nrow(expand.grid.unique(PatientIx,PatientIx))),])
        }
    }
    colnames(InterPatient) <- colnames(IntraPatient)
    Comparisons <- rbind(IntraPatient,InterPatient)
}else if(dataset == 'AUMC'){
    # For the AUMC dataset, we only make Intrapatient comparisons
    # Create empty data frame
    Comparisons <- data.frame()
    for(patient in unique(PatientNames)){
        # Append all pairwise combinations of one patient 
        Comparisons <- rbind(Comparisons, as.data.frame(expand.grid.unique(which(PatientNames == patient),which(PatientNames == patient))))
    }
}

#-------------------------------------------------------------------------------
# 3.1 Calculate clonality statistics of pairs
#-------------------------------------------------------------------------------
# Get reference data
refCNA_split <- get_ref_data('reference/Corrected_Imp_raw_output_Ref08.txt', TRUE)
refCNA_whole <- get_ref_data('reference/Corrected_Imp_raw_output_Ref08.txt', FALSE)
refData_LUAD <- read.table('reference/pfreq_luad_tcga_split.csv', skip=0, fill=TRUE, header=TRUE, sep=",")
refData_LUSC <- read.table('reference/pfreq_lusc_tcga_split.csv', skip=0, fill=TRUE, header=TRUE, sep=",", row.names = 'X')




# Intialize output dataframe
output_df <- data.frame()

# Iterate over comparisons
for(i in 1:nrow(Comparisons)){
    # Subset Segments
    data <- Segments[,c(Comparisons$V1[i],Comparisons$V2[i])]
    
    sampleCNA_split <- get_sample_data(data, TRUE)
    sampleCNA_whole <- get_sample_data(data, FALSE)
    sampleCNA_split <- sampleCNA_split[sampleCNA_split['chrom']!='chr21p',]
    sampleCNA_split <- sampleCNA_split[order(sampleCNA_split[,'chrom']),]
    pjct <- paste0(sampleNames(data),collapse = '-')
    # Get LLR2 for split arms        
    clonality.analysis(sampleCNA_split, ptlist=rep(pjct, ncol(data)), refdata=refCNA_split) -> llrData
    clonalityTest(data, sbjctLst=rep(pjct, ncol(data)), refdata=refCNA_split, llrData=llrData) -> cln
    regular_output <- data.frame(cln$clonTab)
    regular_output$llr2.split <- data.frame(cln$clonTab)$llr2
    regular_output$cor <- 1 - regular_output$cor

    # LUAD (split)
    regular_output$llr2.luad.split <- LLR2_ref(data,refData_LUAD, split=TRUE)
    # LUSC (split)
    regular_output$llr2.lusc.split <- LLR2_ref(data,refData_LUSC, split=TRUE)

    output_df <- rbind(output_df,regular_output)
}


#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
write.table(output_df , output, sep = '\t',quote=F,row.names = F)
