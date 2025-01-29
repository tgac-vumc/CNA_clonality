#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Mutation_clonality_TRACERx.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Perform clonality classification using mutations in TRACERx
# (based on https://github.com/tgac-vumc/Mutation_Clonality/)
#
# a) Filter TRACERx mutations by panel regions (https://github.com/tgac-vumc/Mutation_Clonality/scripts/Filter_mutations.R)
# b) Make same comparisons as CNAs
# c) Apply molecular classification (MC) algorithm (https://github.com/tgac-vumc/Mutation_Clonality/scripts/Call_Clonality_MCalgorithm.R)
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  13-11-2024: File creation, write code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))

#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_segments <- snakemake@input[["Segments"]]
    input_mutations <- snakemake@input[["Mutations"]]
    input_panel <- snakemake@input[["Panel"]]
    output <-  snakemake@output[["Clonality_MC"]]
}else{
    input_segments <- 'output/CNA_Clonality_statistics_TRACERx.txt'
    input_mutations <- 'data/Mutations/Mutations_TRACERx.txt'
    input_panel <- 'reference/InhouseLungPanel.bed'
    output <- 'output/Mutations_Clonality_statistics_TRACERx.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
mutations <- read.delim(input_mutations,stringsAsFactors = F)
segments <- read.delim(input_segments, stringsAsFactors = F)
panel <- read.delim(input_panel, header = F, stringsAsFactors = F)

#-------------------------------------------------------------------------------
# 1.1 Reformat mutations in long format
#-------------------------------------------------------------------------------
# Transform mutation table to indicate presence of mutation in recurrence
mutations <- mutations %>%
    # filter by subtype
    # split rows per recurrence
    tidyr::separate_rows(RegionSum, sep = ';') %>% 
    # recover number of alt reads and region number
    mutate(Region = purrr::map_chr(RegionSum,~strsplit(.x,':')[[1]][1]),
           ALT =  purrr::map_chr(RegionSum,~strsplit(strsplit(.x,':')[[1]][2],'/')[[1]][1]),
           Depth = purrr::map_chr(RegionSum,~strsplit(strsplit(.x,':')[[1]][2],'/')[[1]][2]), 
           VAF = as.numeric(ALT) / as.numeric(Depth),
           chr= paste0('chr',chr),
           # create boolean to indicate presence of mutation
           Mutation_present = ALT != 0,
           NucleotideChange = '',
           AAChange = '') %>%
    select(SampleID,Region,Hugo_Symbol,chr,start,stop,ref,var,func,MutationID,NucleotideChange, AAChange,VAF,Mutation_present)



#-------------------------------------------------------------------------------
# 1.2 Filter mutations by panel
#-------------------------------------------------------------------------------
colnames(panel)[1:3] <- c('chr','chromStart','chromEnd')
# Join panel regions and filter by panel
Filtered_mutations <-
    mutations  %>% 
    inner_join(panel) %>%
    filter(start >= chromStart & start <= chromEnd ) %>%
    select(-chromStart,-chromEnd) 

#-------------------------------------------------------------------------------
# 2.1 Fetch tumor pairs
#-------------------------------------------------------------------------------
Tumor_pairs <- segments %>% select(patient) %>% tidyr::separate(patient,sep = '-', into=c('Tumor1','Tumor2')) %>%
    mutate(Tumor1 = gsub('SU_T1.|SU_T2.|SU_','',Tumor1),Tumor2 = gsub('SU_T1.|SU_T2.|SU_','',Tumor2)) %>% 
    # Retrieve True Clonality
    mutate(patient1 = purrr::map_chr(Tumor1,~strsplit(.x,'_')[[1]][1])
          ,patient2 = purrr::map_chr(Tumor2,~strsplit(.x,'_')[[1]][1]),
           True_Clonality = ifelse(patient1 == patient2,'Clonal','Non-Clonal'))

# Fetch sampleIDs
Samples <- unique(c(Tumor_pairs$Tumor1,Tumor_pairs$Tumor2))



#-------------------------------------------------------------------------------
# 3.1 Reformat mutations
#-------------------------------------------------------------------------------
Filtered_mutations <- Filtered_mutations %>% 
    mutate(MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
           SampleID = paste0(SampleID,'_',Region),
           Mutation_present = as.integer(Mutation_present),
           MutationID = paste0(Hugo_Symbol,'_',MutationID)) %>% 
    select(-c(AAChange,NucleotideChange))



mutations_broad <-
    # Create fields
    Filtered_mutations %>% 
    select(SampleID,MutationID,Mutation_present) %>%
    unique() %>% 
    # Get broad format and fill missing values with 0 (not present)
    tidyr::pivot_wider(id_cols =c(MutationID),
                       values_from = Mutation_present,
                       names_from = SampleID,
                       values_fill = 0)


# fill in missing samples with 0's
Missing_Samples <- Samples[!Samples %in% colnames(mutations_broad)]
for(sample in Missing_Samples){
    mutations_broad[sample] = 0
}
# fetch numeric matrix
mutation_matrix <-
    mutations_broad %>%
    tibble::column_to_rownames(var= 'MutationID') %>%
    as.matrix()


# Fetch oncoKB annotations
oncogenic_driver_mutations <- read.delim('https://raw.githubusercontent.com/tgac-vumc/Mutation_Clonality/refs/heads/main/manifest/OncoKB_drivers.txt', stringsAsFactors = F) %>% 
    # Exclude EGFR resistance mutations for clonality testing
    filter(!(grepl('EGFR',Hugo_Symbol) & grepl('T790M',HGVSp))) %>% 
    mutate(MutationID =  paste0(Hugo_Symbol,'_',paste(gsub('chr', '', Chromosome),Start_Position, Reference_Allele,Tumor_Seq_Allele1, sep = ' '))) %>%
    pull(MutationID) %>% unique()

# manually add mutation with missing annotation
oncogenic_driver_mutations <- c(oncogenic_driver_mutations, 'EGFR_7 55242464 A NA')

#-------------------------------------------------------------------------------
# 3.1 Define MCalgorithm
#-------------------------------------------------------------------------------
MCalgorithm <- function(comparison) {
    # Fetch samplenames
    sample1 <- comparison[,'Tumor1']
    sample2 <- comparison[,'Tumor2']
    # Fetch mutations
   
    data1 <- mutation_matrix[,sample1]
    data2 <- mutation_matrix[,sample2]

    #Fetch oncogenic driver mutations
    #oncogenic1 <- data1[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) %in% oncogenic_driver_mutations]
    #oncogenic2 <- data2[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) %in% oncogenic_driver_mutations]
    
    oncogenic1 <- data1[names(data1) %in% oncogenic_driver_mutations]
    oncogenic2 <- data2[names(data2) %in% oncogenic_driver_mutations]
    
    # Fetch other mutations
    data_no_drivers1 <- data1[!names(data1) %in% names(oncogenic1)]
    data_no_drivers2 <- data2[!names(data2) %in% names(oncogenic2)]

    shared_mutations <- paste(names(data1)[which(data1 == 1 & data2 == 1)],collapse=' --- ')
    
    #-------------------------------------------------------------------------------
                                        # Call clonality
    #-------------------------------------------------------------------------------
    # If any oncogenic mutation is found in any sample
    if(any(oncogenic1 == 1) | any(oncogenic2 == 1)){
        # if samples do not have identical driver mutations it is non-clonal
        if(!identical(oncogenic1,oncogenic2)){
            if(any(grepl('KRAS',names(oncogenic1[oncogenic1 == 1]))) | any(grepl('KRAS',names(oncogenic2[oncogenic2 == 1])))){
                Clonality <- 'Probably Non-Clonal'
                reason <- 'Different KRAS mutation status'
            }else{
                Clonality <- 'Non-Clonal'
                reason = 'Different oncogenic driver mutations found'
            }
            return(
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = reason,Shared_mutations = shared_mutations)
            )
            
        }

    }
    
    # If same oncogenic mutation, or both wild type 
    if(identical(oncogenic1,oncogenic2)){
        # next check if there is more or equal than 1 shared mutation (excluding driver mutations)
        match_bool <- rowSums(matrix(c(data_no_drivers1,data_no_drivers2),ncol=2)) == 2
        n_match <- sum(match_bool)
        # If there is any match mutation it is clonal
        if(any(match_bool)){
            Clonality <- 'Clonal'
            if(all(oncogenic1 == 0) & all(oncogenic2 == 0)){
                reason <- paste0('Wildtype oncogenic but shared mutations' )
            }else{
                reason <- paste0('Same oncogenic driver and shared mutations')
            }
            return(
                data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = reason,Shared_mutations = shared_mutations)
            )
        }else{
            # In the case of no matching mutations:
            TP53_1 <- data1[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            TP53_2 <- data2[purrr::map_chr(names(data1),~strsplit(.x,'_')[[1]][1]) == 'TP53']
            if(all(TP53_1 == 0) & all(TP53_2 == 0)){
                Clonality <- 'Inconclusive'
                return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = Clonality, reason = 'No shared non-oncogenic mutations found and TP53 wildtype',Shared_mutations = shared_mutations))
            }else{
                 Clonality <- 'Probably Non-Clonal'
                 return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC =Clonality, reason = 'No shared non-oncogenic mutations found and different TP53 mutations',Shared_mutations = shared_mutations))
            }
        }
    }
    return(data.frame(Sample1=sample1,Sample2=sample2,True_clonality = comparison$True_Clonality, Clonality_MC = '????', reason = 'Unknown',Shared_mutations = shared_mutations))
}

#-------------------------------------------------------------------------------
# 3.2 Perform Clonality testing 
#-------------------------------------------------------------------------------
# Initialize dataframes
MC_Clonalities <- data.frame()

# Iterate over comparisons
for(i in 1:nrow(Tumor_pairs)){
    # Perform molecular classification algorithm
    MC_Clonalities <-rbind(MC_Clonalities, MCalgorithm(Tumor_pairs[i,]))
}

#----------------------------------------------------------------------
# 3.1 Write to file
#-------------------------------------------------------------------------------
write.table( MC_Clonalities, file = output, sep = '\t', row.names = F,quote = F )
