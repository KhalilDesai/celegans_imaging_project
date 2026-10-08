## Zenodo EMBD public data repository ##

### ABOUT ### 

This folder consists of data generated for the C. elegans embryonic development project that is described here:

2023- Automated Profiling of Gene Function During Embryonic Development Rebecca A. Green1,2†, Renat N. Khaliullin1†, Zhiling Zhao1, Stacy. D. Ochoa1, Jeffrey M. Hendel1, Tiffany-Lynn Chow1, HongKee Moon1, Ronald J. Biggs1, Arshad Desai1 and Karen Oegema1. (submitted to Cell in March 2023)

### DATA  ###
This folder contains:
A. Index file: Use this to look up the EMBD#### unique identifier for a gene of interest.

B. All_genes_mySQL_table:
Contains all fields from the mySQL genes table, which includes descriptive information (gene descriptions, orthologs, etc), primary data gathered through automated analysis, manual scoring data, etc.

C. For each RNAi condition there are 7 files available:
1.EMBD####_Comp: composite movies of all analyzed data per condition
2.EMBD####_tophits: movies of top hits per condition
3.EMBD####_reciprocal_hits: movies of reciprocal hits per condition
4.EMBD####_UMAP_PAD_opt...: UMAP generated with optimal parameters
5.EMBD####_UMAP_PAD_GLS...:UMAP generated with germ-layer parameters
6.EMBD####_UMAP_PAD_MS...:MAP generated with morphogenesis parameters
7. EMBD####_hits_matrix: csv file with top hits listed for optimal, GLS, MS parameter sets


### CONTACT ###
Rebecca Green regreen@ucsd.edu