#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Mutation_clonality_AUMC.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Perform clonality classification using mutations in AUMC
# (based on https://github.com/tgac-vumc/Mutation_Clonality/)
#
# a) Determine True Clonality based on WES using germline filter setp
# b) Filter by panel regions
# c) Manually curate and apply molecular classification (MC) algorithm (https://github.com/tgac-vumc/Mutation_Clonality/scripts/Call_Clonality_MCalgorithm.R)
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
    input_segments <- 'output/CNA_Clonality_statistics_AUMC.txt'
    input_mutations <- 'data/Mutations/Mutations_AUMC.txt'
    input_panel <- 'reference/InhouseLungPanel.bed'
    output <- 'output/Mutations_Clonality_statistics_AUMC.txt'
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
mutations <- read.delim(input_mutations,stringsAsFactors = F)
segments <- read.delim(input_segments, stringsAsFactors = F)
panel <- read.delim(input_panel, header = F, stringsAsFactors = F)

#-------------------------------------------------------------------------------
# 1.2 Deterimine true clonality using number of high confidence somatic variants
#-------------------------------------------------------------------------------
True_clonality <-
    mutations %>%
    group_by(patient) %>%
    # select variants on diploid regions, except Patient7 with tumor ploidy of 4
    filter((ploidy == 2) | (patient == 'Patient7'& ploidy == 4),!is.na(SomaticProbability)) %>%
    # Determine top 50% of high confidence somatic variants
    mutate(IsSomatic = SomaticProbability >= quantile(SomaticProbability,0.5)) %>%
    # Select variants in intersection
    filter(Present == 'Intersection')%>%
    # Count number of somatic variants
    group_by(patient,IsSomatic) %>%
    summarise(N= dplyr::n()) %>%
    tidyr::pivot_wider(names_from = 'IsSomatic',values_from = 'N',names_prefix='Somatic:') %>%
    replace(is.na(.),0)  %>%
    ungroup() %>% 
    select(patient, colnames(.)[3]) %>%
    # Determine true clonality using Liu et al cutoff
    mutate(True_clonality = ifelse(`Somatic:TRUE` <= 2,'Non-Clonal','Clonal')) %>%
    select(patient,True_clonality)

#-------------------------------------------------------------------------------
# 1.3 Filter mutations by panel
#-------------------------------------------------------------------------------
Samples <- unique(mutations$sample)
colnames(panel)[1:3] <- c('chr','chromStart','chromEnd')

# Join panel regions and filter by panel
Filtered_mutations <-
    mutations %>%
    mutate(
        chr = purrr::map_chr(mut,~strsplit(.x,'_')[[1]][1]),
        start = purrr::map_chr(mut,~strsplit(.x,'_')[[1]][2]),
        ref = purrr::map_chr(mut,~strsplit(.x,'_')[[1]][3]),
        var = purrr::map_chr(mut,~strsplit(.x,'_')[[1]][4])) %>%
    inner_join(panel) %>%
    filter(start >= chromStart & start <= chromEnd ) %>%
    select(-chromStart,-chromEnd) 

#-------------------------------------------------------------------------------
# 3.1 Reformat mutations
#-------------------------------------------------------------------------------
Filtered_mutations <-
    Filtered_mutations %>%
    mutate(MutationID = paste(gsub('chr', '', chr),start, ref,var, sep = ' '),
           SampleID = sample,
           Mutation_present = as.integer(1),
           MutationID = paste0(Gene,'_',MutationID)) %>%
    arrange(sample) %>%
    unique() %>%
    # Manual curation,
    filter(
        # Remove suspected germline variant
        !(patient == 'Patient4' & mut == "chr17_7578406_C_T" )) # AF 0.647 and 0.882


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

# Read OncoKB annotations
oncogenic_driver_mutations <- read.delim('https://raw.githubusercontent.com/tgac-vumc/Mutation_Clonality/refs/heads/main/manifest/OncoKB_drivers.txt', stringsAsFactors = F) %>% 
    # Exclude EGFR resistance mutations for clonality testing
    filter(!(grepl('EGFR',Hugo_Symbol) & grepl('T790M',HGVSp))) %>% 
    mutate(MutationID =  paste0(Hugo_Symbol,'_',paste(gsub('chr', '', Chromosome),Start_Position, Reference_Allele,Tumor_Seq_Allele1, sep = ' '))) %>%
    pull(MutationID) %>% unique()

#-------------------------------------------------------------------------------
# 3.1 Define Tumor pairs
#-------------------------------------------------------------------------------
Tumor_pairs <- data.frame(Tumor1 = Samples[grepl('_1$',Samples)],
                          Tumor2 = Samples[grepl('_2$',Samples)],
                          True_Clonality = True_clonality$True_clonality)


Tumor1 <- 'Patient1_1'
Tumor2 <- 'Patient1_2'


#-------------------------------------------------------------------------------
# 3.1 Define MCalgorithm
#-------------------------------------------------------------------------------
MCalgorithm <- function(comparison) {
    # Fetch samplenames
    sample1 <- as.character(comparison[,'Tumor1'])
    sample2 <- as.character(comparison[,'Tumor2'])
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
