**VIBE: “Visualization and Exploration of Bulk mRNA Expression data**


VIBE revolutionizes the process of selecting tumor types for monoclonal and bispecific antibodies in the drug discovery journey. VIBE's intuitive user interface and interactive visualization empower researchers and clinicians to make data-driven decisions by exploring gene expression patterns in bulk RNA sequencing data.

With VIBE, researchers can delve into the expression of specific genes or gene pairs within large transcriptomic datasets such as the publicly available TCGA and GTEx databases as well as in house or commercially acquired datasets. Moreover, VIBE can compute composite scores for gene sets, such as pathway-associated genes, thereby offering a comprehensive overview of gene interaction and function within the pathway of interest. The resultant data can be visualized through intuitive boxplots and heatmaps, providing an accessible and informative platform for the exploration and interpretation of complex genomic data. 

***Install VIBE package***


install.packages("devtools")  
devtools::install_github("genmab/VIBE")

***Implementation***

![VIBE_Implementation_Overview](https://github.com/genmab/VIBE/assets/139466232/33319a98-3493-44d7-845e-dd89993aebf5)



***Implementation and overview of the statistics and visualization capabilities of the VIBE package.***


(A) The df_harmonize() function harmonizes dummy data, creating a structured data frame with appropriate column names and format for VIBE functionalities. The harmonize_df_pathway() function combines genes into pathways for VIBE visuals. 
(B) The get_thresholds() function defines gene-specific thresholds (mean, median or 75% quantile (Q3)) for expression analysis and indication selection. 
(C) The calculate_percentage_per_quadrant() function extends get_thresholds() for dual-gene analysis, displaying the percentage of samples in quadrants, correlation, and statistical significance (Spearman correlation). 
(D) plot_expression_box() presents a boxplot with sample distribution, thresholds, percentage satisfying thresholds, and user-defined grouping variables and groups (i.e. “intr”) generated when harmonizing the dataframe. 
(E) plot_expression_scatter() plots a scatterplot with user-selected grouping variables and groups and statistics calculated by the calculate_percentage_per_quadrant() function. 
(F) heatmap_sample_expression() visualizes gene expression via heatmaps, allowing grouping variable selection, group plots, gene splitting, and ordering of genes. 
(G) heatmap_samples_above_median() represents the percentage of samples satisfying the threshold for all the genes with similar functionalities as mentioned in heatmap_sample_expression(). 
(H) plot_expression_box_split() function visualizes the expression of Gene 1 grouped by additional grouping columns with statistics (% above median for each group, the fold change (FC) and Kruskal-Wallis p-value (KW p)) above. 

***Tutorial***
Follow our vignettes to visualize your favourite datasets.


***Contributors***
Indu Khatri (inkh@genmab.com)  
Saskia van Asten (saas@genmab.com) 
Leandro Moreno (lemo@genmab.com)  
Francis Blokzijl (frbl@genmab.com)  
Iris Kolder (iko@genmab.com)  


***Contact***
Iris Kolder (iko@genmab.com) 
