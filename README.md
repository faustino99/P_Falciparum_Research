# P.-falciparum-research

This project uses Gene Set Enrichment Analysis (GSEA) to conduct a meta-analysis on published gene expression data for the P. falciparum parasite, one of the most prevalent malaria-transmitting organisms in the world. The following describes what the files in this repository do:  


This folder contains two main scripts: **Gene_Set_Script.R** and **GSEA_heatmap_script.R**  

**Gene_Set_Script.R** creates genesets and gene set collections (.gmt files) for the P.falciparum parasite based on localization and peak gene expression during the asexual blood cycle and gametocyte stages. The localization genesets and the .gmt file of all these genesets are located in the localization_gene_data folder. the gene expression gene sets and corresponding .gmt files are located in the **gene_expression_data** folder.  

  The following files contain the data needed to make the genesets:  
    – **Plasmodium_falciparum.csv** - Data table of P.falciparum genes and their organelle locations within the parasite, created using data from the APILOC database.  
    – **Pf3D7_gene_aliases.txt** - Data table of most recent gene symbols for each gene used to update outdated gene names from Plasmodium_falciparum.csv. Created using data from PlasmoDB database.  
    – **3D7_Winzeler_Gametocyte.txt** - Data table of gene expression during the gametocyte cycle. This data is referenced as Winzeler's genesets for all the scripts.  
    – **profiles.diff3.txt** - Data table of gene expression during the gametocyte cycle and the blood stages. This data is references as Su's genesets.  
    – **profiles.diff2.txt** - Data table of gene expression during the blood stages. Referred to as Sunnenberg's genesets.  
    – **profiles.diff.txt** - Data of gene expression during the blood stages. Referred to as Bartfai's genesets.  
    transmission expression.txt - Data of gene expression during the blood stages referred to as Derisi's genesets.  
    
    Note: all of the gene expression datasets were taken from PlamoDB.  
    
**Gene_Set_Script.R** also creates an **Images** folder that contains graphs of the sizes of each gene expression geneset over time. These graphs are good indicators of cell activity and show which peek expression cutoff (mean+1.25 or mean+1.5 SD) has the most genesets within the 50-500 gene size range that is most suitable for GSEA to work (500 and 50 are marked by dotted lines in each graph).  
  
  
This folder also contains files that allow GSEA to be run with two data sets that are comparing the genomic profile of P.falciparum when treated with a specific drug:  
  
  The first experiment (Baum J, Maier AG, Good RT, Simpson KM et al. PLoS Pathog 2005 Dec;1(4):e37.) analyzes the genomic profiles of 3 different variants of the parasite when treated and not treated with chloroquine (CQ).  
    – **pfgsm2532RNA_3.cls** is the cls file to be used for GSEA  
    – **PfGSM2532RNA.gct** is the gct file to be used for GSEA  
  
  The second experiment (Tarr SJ, Nisbet RE, Howe CJ. Mol Biochem Parasitol 2011 Sep;179(1):37-41. ) analyzes the genomic profile of the parasite when treated and not treated with the Thiostrepton antibiotic.  
    – **GSM71079x.cls** & **GSM71079x.gct** are the cls and gct files need for this experiment for GSEA  
    
  **Gene_Set_Script.R** creates a chip file for GSEA called **pfchip.chip** that is needed for any experiment using the Affymetrix Plasmodium/Anopheles Genome Array. The two experiments above use this array.  
    – The **Pf_microarray_genesymbols.csv** file (taken from Affymetrix website) is used to create the chip file.  
      
  **Gene_Set_Script.R** also creates .gmt files for each gene set collection made, the file format required for GSEA.  
    
    
**GSEA_heatmap script.R** has a function that allows for a heatmap to be made from the GSEA results of each geneset collection that reflects the Normalized Enrichment Score (NES) value for each geneset-phenotype comparison made.  

  The **GSEA_Results** folder contains heatmap images of the GSEA results from the first experiment when run with each geneset collection.  
  
  The **GSEA_Result_2** folder contains heatmap images of the GSEA results from the second experiment when run with each geneset collection.  
    
  An analysis of the results from this research is discussed in depth in **Research_Paper.pdf**
