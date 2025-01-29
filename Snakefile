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
        
#+++++++++++++++++++++++++++++++++++ 1 CNA CLONALITY STATISTICS  ++++++++++++++++++++++++++++++++++
# 1.1 Calculate LogLikeLihood and 
rule Calculate_Clonality_statistics:
    input:
        Segments = data_dir + 'CNA/SupplementaryData_{dataset}.Rds',
        SampleTable = 'reference/SampleTables/SampleTable_{dataset}.csv'
    output:
        Clonality_statistics = output_dir + 'CNA_Clonality_statistics_{dataset}.txt'
    conda:
       "envs/Clonality.yaml"
    script:
        "scripts/Calculate_Clonality_statistics.R"


rule Determine_NGSpanel_clonality:
    input:
        Segments = output_dir + 'Clonality_statistics_{dataset}.txt',
        Mutations = data_dir + 'Mutations/Mutations_{dataset}.txt',
        Panel = 'reference/InhouseLungPanel.bed'
    output:
        Clonality_MC = output_dir + 'Mutations_Clonality_statistics_{dataset}.txt'
    conda:
        'envs/R.yaml'
    script:
        "scripts/Mutation_clonality_{dataset}.R"
        

#+++++++++++++++++++++++++++++++++++ 2 CLONALITY CLASSIFICATION  +++++++++++++++++++++++++++++++++
# 2.1 Classify clonality and plot
rule Clonality_Classification:
    input:
        Clonality_statistics = output_dir + 'CNA_Clonality_statistics_{dataset}.txt',
        Clonality_MC = output_dir + 'Mutations_Clonality_statistics_{dataset}.txt',
        SampleTable = 'reference/SampleTables/SampleTable_{dataset}.csv'
    output:
        Clonality_classification = 'plots/Clonality_classification_{dataset}.pdf',
        SankeyPlot = 'plots/SankeyPlot_{dataset}.pdf',
    conda:
        'envs/R.yaml'
    script:
        'scripts/Clonality_classification.R'
        
