configfile: "config.yaml"
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
datasets = config["all"]["datasets"]
#--------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand('plots/Clonality_classification_{dataset}.pdf',dataset=datasets)
        
#+++++++++++++++++++++++++++++++++++++ 1 CLONALITY STATISTICS  ++++++++++++++++++++++++++++++++++++
# 1.1 Calculate LogLikeLihood and 
"""
rule Calculate_Clonality_statistics:
    input:
        Segments = data_dir + 'SupplementaryData_{dataset}.Rds',
        SampleTable = 'reference/SampleTables/SampleTable_{dataset}.csv'
    output:
        Clonality_statistics = output_dir + 'Clonality_statistics_{dataset}.txt'
    conda:
       "envs/Clonality.yaml"
    script:
        "scripts/Calculate_Clonality_statistics.R"
"""
#+++++++++++++++++++++++++++++++++++ 2 CLONALITY CLASSIFICATION  ++++++++++++++++++++++++++++++++++
# 2.1 Train GMM using TRACERx clonality statistics
rule Train_GMM:
    input:
        Clonality_statistics = output_dir + 'Clonality_statistics_TRACERx.txt',
        SampleTable = 'reference/SampleTables/SampleTable_TRACERx.csv'
    output:
        GMM_model = output_dir + 'GMM_model.Rds',
        GMM_summary = 'plots/GMM_summary.pdf'
    conda:
        'envs/GMM.yaml'
    script:
        'scripts/Train_GMM.R'
# ------------------------------------------------------------------------------------------------        
# 2.2 Classify clonality
rule Clonality_Classification:
    input:
        Clonality_statistics = output_dir + 'Clonality_statistics_{dataset}.txt',
        GMM_model = output_dir + 'GMM_model.Rds',
        SampleTable = 'reference/SampleTables/SampleTable_{dataset}.csv'
    output:
        Clonality_classification = 'plots/Clonality_classification_{dataset}.pdf',
    conda:
        'envs/GMM.yaml'
    script:
        'scripts/Clonality_classification.R'
        
