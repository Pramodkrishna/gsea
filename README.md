# gsea
Gene Set Enrichment Analysis (GSEA) is a powerful method for interpreting gene expression data. Instead of looking at individual genes it focuses on groups of genes that share a common biological function, chromosomal location, or regulation. This implementation of GSEA follows the description provided in Subramanian et al., Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, PNAS 102 (43): 15545-15550, 2005.

The source code is written in Python and resides in `gsea.py`. The three most important functions are:
* `gsea` which computes normalized enrichment scores and the corresponding _p_-values based on the data in gene expression profiles,
* `plot` which displays Fig. 1 from the original paper: 
  * a heat map of gene expression profiles grouped by class,
  * gene expression correlation profile,
  * gene set member position in the profiles, 
  * the running sum (see the original paper),
* `my_gsea` which computes estimates of normalized enrichment scores and the corresponding _p_-values without gene expression profiles. A user-provided set of interesting genes is used instead.

Besides the source code three example input files are provided:
* `leukemia.txt` contains gene expression profiles for two types of leukemia: acute lymphocytic leukemia (ALL) and acute myelogenous leukemia (AML),
* `pythways.txt` lists _a priori_ defined sets of genes enoding products in metabolic pathways,
* `my_genes.txt` contains a set of interesting genes from the literature.
