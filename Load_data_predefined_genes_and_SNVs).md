---
editor_options: 
  markdown: 
    wrap: 72
---

# Prepare data

The following is a Rmarkdown tutorial for running similar models found in the article:

Engvall, K., Uvdal, H., Björn, N., Åvall-Lundqvist, E., & Gréen, H. (2024). Prediction models of persistent taxane-induced peripheral neuropathy among breast cancer survivors using whole-exome sequencing. NPJ Precision Oncology, 8(1). https://doi.org/10.1038/s41698-024-00594-x


In this script we prepare the literature data and our own data for
further analysis.

# Predefined SNVs and genes

Here we load the predefined literature data found in KEGG correlated to
the Taxane pathway, Axon Guidance, Axon Regeneration and Pathways of
Neurodegeneration. Based on a literature search we also load previous
variants and genes correlated to Taxane-Induced Peripheral Neuropathy
(TIPN).

```{r}
taxane_pathway <- read.table("./337/predefined_snvs_and_genes/taxane_pathway.txt", sep = "\t", stringsAsFactors = F, header = T)

prev_variants <- read.table("./337/predefined_snvs_and_genes/prev_variants.txt", sep = "\t", stringsAsFactors = F, header = T)
prev_genes <- read.table("./337/predefined_snvs_and_genes/prev_genes.txt", sep = "\t", stringsAsFactors = F, header = T)

axon_guidance <- read.table("./337/predefined_snvs_and_genes/axon_guidance_clean.txt", sep = "\t", stringsAsFactors = F, header = T)
axon_regeneration <- read.table("./337/predefined_snvs_and_genes/axon_regeneration_clean.txt", sep = "\t", stringsAsFactors = F, header = T)
neurodegeneration <- read.table("./337/predefined_snvs_and_genes/pathways_of_neurodegeneration_clean.txt", sep = "\t", stringsAsFactors = F, header = T)

```

## Genes and variants in our data

Here we upload the common and rare genes found in the gene/region-based
association assay in our own cohort data.

```{r}
common_and_rare_variants_in_genes <- read.table("./337/gene_tests/common_and_rare_variants_in_genes_20210810.txt", 
                                                sep = "\t", header = T, stringsAsFactors = F) #gjord med gene_regions.R
common_variants_in_genes <- read.table("./337/gene_tests/common_variants_in_genes.txt", 
                                                sep = "\t", header = T, stringsAsFactors = F) #gjord med gene_regions.R

```

Here we look through the data. Which genes and variants are overlapping
in our common and rare variant set?

```{r}
common_and_rare_variants_in_genes
common_and_rare_variants_in_genes[, c(1,5)]

dim(taxane_pathway)
sum(taxane_pathway$gene %in% common_and_rare_variants_in_genes$gene_name)

dim(prev_variants)
sum(prev_variants$variant %in% common_and_rare_variants_in_genes$variant)

dim(prev_genes)
sum(prev_genes$gene %in% common_and_rare_variants_in_genes$gene_name)

dim(axon_guidance)
sum(axon_guidance$gene %in% common_and_rare_variants_in_genes$gene_name)

dim(axon_regeneration)
sum(axon_regeneration$gene %in% common_and_rare_variants_in_genes$gene_name)

dim(neurodegeneration)
sum(neurodegeneration$gene %in% common_and_rare_variants_in_genes$gene_name)
```

## Extract variants from predefined genes

Here we extract the rare and common variants found in each predefined
literature gene and variant set.

```{r}
taxane_pathway_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% taxane_pathway$gene), c(1,5)]

taxane_pathway_genes_and_common_variants <- common_variants_in_genes[which(common_variants_in_genes$gene_name %in% taxane_pathway$gene), c(1,5)]

prev_genes_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% prev_genes$gene), c(1,5)]

axon_guidance_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% axon_guidance$gene), c(1,5)]
axon_regeneration_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% axon_regeneration$gene), c(1,5)]
neurodegeneration_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% neurodegeneration$gene), c(1,5)]
```

In case we want to save this results for later:

```{r}
write.csv(neurodegeneration_genes_and_variants, "./337/predefined_snvs_and_genes/cpdb predef variable.csv")
```

## Prepare and load our own cohort data

Read in SNVs results from association analyses:

```{r}

numb_snvs <- read.table("./HP/raw_snv_result_tables_20211031/numb_snps_table3.txt", header = T, stringsAsFactors = F, sep = "\t") ting_snvs <- read.table("./HP/raw_snv_result_tables_20211031/ting_snps_table3.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_snvs <- read.table("./HP/raw_snv_result_tables_20211031/cramps_snps_table3.txt", header = T, stringsAsFactors = F, sep = "\t") jar_snvs <- read.table("./HP/raw_snv_result_tables_20211031/jar_snps_table3.txt", header = T, stringsAsFactors = F, sep = "\t") weak_snvs <- read.table("./HP/raw_snv_result_tables_20211031/weak_snps_table3.txt", header = T, stringsAsFactors = F, sep = "\t") small_snvs <- read.table("./HP/raw_snv_result_tables_20211031/small_snps_suptable1.txt", header = T, stringsAsFactors = F, sep = "\t")` 



```

Read in the gene results from gene tests:

```{r}
numb_genes <- read.table("./HP/raw_gene_result_tables_20211031/numb_genes_0.001_table4.txt", header = T, stringsAsFactors = F) ting_genes <- read.table("./HP/raw_gene_result_tables_20211031/ting_genes_0.001_table4.txt", header = T, stringsAsFactors = F) cramps_genes <- read.table("./HP/raw_gene_result_tables_20211031/cramps_genes_0.001_table4.txt", header = T, stringsAsFactors = F) jar_genes <- read.table("./HP/raw_gene_result_tables_20211031/jar_genes_0.001_table4.txt", header = T, stringsAsFactors = F) weak_genes <- read.table("./HP/raw_gene_result_tables_20211031/weak_genes_0.001_table4.txt", header = T, stringsAsFactors = F) small_genes <- read.table("./HP/raw_gene_result_tables_20211031/small_genes_0.001_suptable2.txt", header = T, stringsAsFactors = F)
```

Read ConsensusPathwayDataBase (CPDB) results. Here we load CPDB results
based on 50 SNVs, optimized number of SNVs based on FDR, 120 SNVs, 210
SNVs, 400 SNVs and common genes in the taxane pathway:

```{r}
#Based on 50 SNVs 
numb_50cpdb <- read.table("./HP/pathways/numb_50snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_50cpdb <- read.table("./HP/pathways/ting_50snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_50cpdb <- read.table("./HP/pathways/cramps_50snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") jar_50cpdb <- read.table("./HP/pathways/jar_50snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") weak_50cpdb <- read.table("./HP/pathways/weak_50snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

#Based on optimized SNVs from FDR 
numb_cpdb <- read.table("./HP/pathways/numb_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_cpdb <- read.table("./HP/pathways/ting_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_cpdb <- read.table("./HP/pathways/cramps_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") jar_cpdb <- read.table("./HP/pathways/jar_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") weak_cpdb <- read.table("./HP/pathways/weak_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") small_cpdb <- read.table("./HP/pathways/small_snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

#Based on 120 SNVs 
numb_120cpdb <- read.table("./HP/pathways/numb_120snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_120cpdb <- read.table("./HP/pathways/ting_120snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_120cpdb <- read.table("./HP/pathways/cramps_120snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") jar_120cpdb <- read.table("./HP/pathways/jar_120snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") weak_120cpdb <- read.table("./HP/pathways/weak_120snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

#Based on 210 SNVs 
numb_210cpdb <- read.table("./HP/pathways/numb_210snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_210cpdb <- read.table("./HP/pathways/ting_210snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_210cpdb <- read.table("./HP/pathways/cramps_210snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") jar_210cpdb <- read.table("./HP/pathways/jar_210snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") weak_210cpdb <- read.table("./HP/pathways/weak_210snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

#Based on 400 SNVs 
numb_400cpdb <- read.table("./HP/pathways/numb_400snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_400cpdb <- read.table("./HP/pathways/ting_400snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") cramps_450cpdb <- read.table("./HP/pathways/cramps_450snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") jar_450cpdb <- read.table("./HP/pathways/jar_450snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") weak_450cpdb <- read.table("./HP/pathways/weak_450snps_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

#Based on common taxane pathway 
taxane_pathway_genes_and_common_variants_cpdb <- read.table("./337/predefined_snvs_and_genes/taxane_common_pathway_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") ting_210_taxcom_cpdb <- read.table("./HP/pathways/ting_210snps_taxcommon_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t") numb_71_taxcom_cpdb <- read.table("./HP/pathways/numb_snps_taxcommon_cpdb.txt", header = T, stringsAsFactors = F, sep = "\t")

```

Read in a summarized version of the published Meta analysis/review both
significant SNVs for each symptom and all SNVs found for each symptom in
the review:

```{r}
#only siginifcant SNVs
numb_signreview <- read.table("./HP/raw_snv_result_tables_20220926/numb_snps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
ting_signreview <- read.table("./HP/raw_snv_result_tables_20220926/ting_snps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
cramps_signreview <- read.table("./HP/raw_snv_result_tables_20220926/cramps_snps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
jar_signreview <- read.table("./HP/raw_snv_result_tables_20220926/jar_snps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
weak_signreview <- read.table("./HP/raw_snv_result_tables_20220926/weak_snps_review.txt", header = T, stringsAsFactors = F, sep = "\t")

#all SNVs found
numb_allreview <- read.table("./HP/raw_snv_result_tables_20220926/numb_allsnps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
ting_allreview <- read.table("./HP/raw_snv_result_tables_20220926/ting_allsnps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
cramps_allreview <- read.table("./HP/raw_snv_result_tables_20220926/cramps_allsnps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
jar_allreview <- read.table("./HP/raw_snv_result_tables_20220926/jar_allsnps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
weak_allreview <- read.table("./HP/raw_snv_result_tables_20220926/weak_allsnps_review.txt", header = T, stringsAsFactors = F, sep = "\t")
```
