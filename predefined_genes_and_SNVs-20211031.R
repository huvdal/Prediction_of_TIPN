####################predefined_snvs_and_genes
taxane_pathway <- read.table("./337/predefined_snvs_and_genes/taxane_pathway.txt", sep = "\t", stringsAsFactors = F, header = T)

prev_variants <- read.table("./337/predefined_snvs_and_genes/prev_variants.txt", sep = "\t", stringsAsFactors = F, header = T)
prev_genes <- read.table("./337/predefined_snvs_and_genes/prev_genes.txt", sep = "\t", stringsAsFactors = F, header = T)

axon_guidance <- read.table("./337/predefined_snvs_and_genes/axon_guidance_clean.txt", sep = "\t", stringsAsFactors = F, header = T)
axon_regeneration <- read.table("./337/predefined_snvs_and_genes/axon_regeneration_clean.txt", sep = "\t", stringsAsFactors = F, header = T)
neurodegeneration <- read.table("./337/predefined_snvs_and_genes/pathways_of_neurodegeneration_clean.txt", sep = "\t", stringsAsFactors = F, header = T)

####################genes and variants in our data
common_and_rare_variants_in_genes <- read.table("./337/gene_tests/common_and_rare_variants_in_genes_20210810.txt", 
                                                sep = "\t", header = T, stringsAsFactors = F) #gjord med gene_regions.R
common_variants_in_genes <- read.table("./337/gene_tests/common_variants_in_genes.txt", 
                                                sep = "\t", header = T, stringsAsFactors = F) #gjord med gene_regions.R

###################################################
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

##################################################
#extract variants from predefined genes
taxane_pathway_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% taxane_pathway$gene), c(1,5)]

taxane_pathway_genes_and_common_variants <- common_variants_in_genes[which(common_variants_in_genes$gene_name %in% taxane_pathway$gene), c(1,5)]

prev_genes_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% prev_genes$gene), c(1,5)]

axon_guidance_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% axon_guidance$gene), c(1,5)]
axon_regeneration_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% axon_regeneration$gene), c(1,5)]
neurodegeneration_genes_and_variants <- common_and_rare_variants_in_genes[which(common_and_rare_variants_in_genes$gene_name %in% neurodegeneration$gene), c(1,5)]

##################################################

#write.csv(neurodegeneration_genes_and_variants, "./337/predefined_snvs_and_genes/cpdb predef variable.csv")

