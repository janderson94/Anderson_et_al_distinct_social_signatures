![](./misc/baboons.png)

# Github repo for Anderson et al. 2021 (_in review_)

## Github repository for results and figures from Anderson et al. "Distinct gene regulatory signatures of dominance rank and social bond strength in wild baboons".

#### This github repo contains 5 major chunks:
##### 1. __preprocessing.R__ -- This script contains pre-processing of raw counts. This code can be skipped, as residual gene expression is also provided within this repo.
##### 2. __linear_models.R__ -- This script contains the linear mixed models run to obtain gene-by-gene estimates of rank and social bond strength effects on gene expression and gene expression responses. 
##### 3. __predictions.R__ -- This script contains elastic net regularization code for predicting dominance rank from gene expression values. 
##### 4. __gene_set_enrichment_analyses.R__ -- This script contains code for Gene Set Enrichment Analyses. The original GSEA scripts which I wrote a wrapper on top of are available at https://github.com/GSEA-MSigDB/GSEA_R. 
##### 5. __figures.R__ -- uploaded soon.
