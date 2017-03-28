# gsea

## Introduction
Gene Set Enrichment Analysis (GSEA) is a powerful method for interpreting gene expression data. Instead of looking at individual genes it focuses on groups of genes that share a common biological function, chromosomal location, or regulation. This implementation of GSEA follows the description provided in Subramanian et al., Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, PNAS 102 (43): 15545-15550, 2005.

The source code is written in Python and resides in `gsea.py`. The three most important __functions__ are:
* `gsea` which computes normalized enrichment scores and the corresponding _p_-values based on the data in gene expression profiles,
* `plot` which displays Fig. 1 from the original paper: 
  * a heat map of gene expression profiles grouped by class,
  * gene expression correlation profile,
  * gene set member position in the profiles, 
  * the running sum (see the original paper).
* `my_gsea` which computes estimates of normalized enrichment scores and the corresponding _p_-values without gene expression profiles. A user-provided set of interesting genes is used instead.

Besides the source code three example __input files__ are provided:
* `leukemia.txt` contains gene expression profiles for two types of leukemia: acute lymphocytic leukemia (ALL) and acute myelogenous leukemia (AML),
* `pythways.txt` lists _a priori_ defined sets of genes enoding products in metabolic pathways,
* `my_genes.txt` contains a set of interesting genes from the literature.

`figure_1.pdf` shows an __example of output__ of the `plot` function for MAP00480_Glutathione_metabolism, one of the high scoring gene sets.

## Analytical approach to computation of enrichment score
Let's compute the enrichment score of a gene set containing a single gene. If this gene is important, i.e. at the top or at the bottom of the ranked gene list, the corresponding enrichment score ±1. If, on the other hand, this gene is not important and appears in the middle of the list, the corresponding enrichment score is going to be near 1/2. After permuting gene expression profile classes, this single-gene set can appear anywhere in the ranked list, yielding an average enrichment score of 3/4.

## Computation of enrichment score from gene expression profiles
After computing Pearson correlation coefficients from gene expression profiles and the corresponding classes, the running sum (see the original paper) is calculated, from which the enrichment score is obtained. This implementation scales linearly with the number of genes _N_ in gene expression profiles. It turns out that this approach provides the most efficient solution for small datasets, like the ones given in the example input files.

For large datasets, another approach is more efficient. The explicit computation of the running sum can be avoided by realizing that the extreme value of the running sum along the genes in gene expression profiles can occur only at positions of genes in the gene set for which the enrichment score is being calculated. Therefore, the running sum needs to be computed only at those positions, thereby reducing the overall computational effort. In this case gene expression profiles need to be structured as dicts instead of lists, resulting in computational efficiency of O(1) instead of O(_N_) in terms of the number of genes in gene expression profiles.

## Computation of enrichment score without gene expression profiles
