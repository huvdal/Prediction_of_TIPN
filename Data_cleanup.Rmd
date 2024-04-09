---
editor_options: 
  markdown: 
    wrap: 72
---

# Clean up data for building and testing logistic regression models for prediction of TIPN

## Load SNVs and gene results

To be able to perform both building and testing we need a lot of data,
all loaded below.

Packages:

```{r}
library(SNPRelate)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ROCR)
library(plotly)
library(ggfortify)
library(ggforce)
library(dplyr)

library(see)
library("shiny")

library(caret)
library(boot)
library (ROCit)
library(MLmetrics)

```

## Load training data

```{r}

#read phenotype data for training data
training_phenotype <- read.table("./337/training_pheno_20211031.txt", sep = "\t", header = T, stringsAsFactors = F)

#read covariate data for training data
training_covariables <- read.table("./337/training_covar_20211031.txt", sep = "\t", header = T, stringsAsFactors = F)

#open
genofile_training <- snpgdsOpen("./337/plink/training/common_and_rare_hwe_target_337_training.gds")

#info
genofile_training

#snps
snps_training <- read.gdsn(index.gdsn(genofile_training, "snp.id"))

#samples
samples_training <- read.gdsn(index.gdsn(genofile_training, "sample.id"))

#alleles
#the reference/non-reference alleles
alleles_training <- read.gdsn(index.gdsn(genofile_training, "snp.allele")) 
```

There are possible values stored in the variable genmat: 0, 1, 2 and
other values. "0" indicates two B alleles, "1" indicates one A allele
and one B allele, "2" indicates two A alleles, and other values indicate
a missing genotype.The function returns a matrix with values 0, 1, 2
representing the number of reference alleles.

```{r}
#genotypes 
training_genotypes <- read.gdsn(index.gdsn(genofile_training, "genotype"))
dim(training_genotypes)
table(training_genotypes) 
colnames(training_genotypes) <- snps_training

#close
snpgdsClose(genofile_training)
```

## Load testing data

```{r}
#read phenotype data for testing data
testing_phenotype <- read.table("./337/testing_pheno_20211031.txt", sep = "\t", header = T, stringsAsFactors = F)

#read covariate data for testing data
testing_covariables <- read.table("./337/testing_covar_20211031.txt", sep = "\t", header = T, stringsAsFactors = F)

#open
genofile_testing <- snpgdsOpen("./337/plink/testing/common_and_rare_hwe_target_337_testing.gds")

#info
genofile_testing

#snps
snps_testing <- read.gdsn(index.gdsn(genofile_testing, "snp.id"))

#samples
samples_testing <- read.gdsn(index.gdsn(genofile_testing, "sample.id"))

#alleles
#the reference/non-reference alleles
alleles_testing <- read.gdsn(index.gdsn(genofile_testing, "snp.allele")) 
```

There are possible values stored in the variable genmat: 0, 1, 2 and
other values. "0" indicates two B alleles, "1" indicates one A allele
and one B allele, "2" indicates two A alleles, and other values indicate
a missing genotype.The function returns a matrix with values 0, 1, 2
representing the number of reference alleles

```{r}
#genotypes 
testing_genotypes <- read.gdsn(index.gdsn(genofile_testing, "genotype"))
dim(testing_genotypes)
table(testing_genotypes) 
colnames(testing_genotypes) <- snps_testing

#close
snpgdsClose(genofile_testing)
```

## Check variables and their correlations

The function below can be used to look for correlations between the
cofactors.

```{r}
plotting_variables<- function(variable,factor1,factor2) {
  variable[which(variable == -9)] <- "Missing"
  variable[which(variable == 0)] <- as.character(factor1)
  variable[which(variable  == 1)] <- as.character(factor2)

  return(variable)
}

check_variables <- function(training_genotypes,training_phenotype,
                            testing_genotypes, testing_phenotype) {

  train_data <- cbind(training_phenotype[1], training_phenotype[4:11],training_phenotype[27],training_phenotype[23],
                      training_phenotype[31],training_phenotype[39],training_phenotype[43],training_phenotype[47])

  #checking continous variables in set
  continuous <-select_if(train_data, is.numeric)
  summary(continuous) #-> factors age,BMI,paclitaxel,Docetaxel and Taxan_type with other scales

  #variable distribution
  library(ggplot2)
  ##symptoms
  p1<-ggplot(continuous, aes(x = Numbness_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p2<-ggplot(continuous, aes(x = Tingling_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p3<-ggplot(continuous, aes(x = Cramps_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p4<-ggplot(continuous, aes(x = Jar_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p5<-ggplot(continuous, aes(x = Weakness_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p6<-ggplot(continuous, aes(x = Maximal_A1_34_C0_12)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p<-ggarrange(p1,p2,p3,p4,p5,p6,labels = c("A", "B", "C", "D","E","F"))
  print(p)
  
  ##prediction factors
  p1<-ggplot(continuous, aes(x = Age)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p2<-ggplot(continuous, aes(x = Diabetes)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p3<-ggplot(continuous, aes(x = BMI_survey)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p4<-ggplot(continuous, aes(x = Paclitaxel_mgm2)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p5<-ggplot(continuous, aes(x = Docetaxel_mgm2)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p6<-ggplot(continuous, aes(x = Taxan_type)) +
    geom_density(alpha = .2, fill = "#FF6666")
  p<-ggarrange(p1,p2,p3,p4,p5,p6,labels = c("A", "B", "C", "D","E","F"))
  print(p)

  ##variable plotting
  train_data$Diabetes<-plotting_variables(train_data$Diabetes,"No Diabetes","Diabetes")
  train_data$Numbness_A1_34_C0_12<-plotting_variables(train_data$Numbness_A1_34_C0_12,"No Numbness","Numbness")
  train_data$Tingling_A1_34_C0_12<-plotting_variables(train_data$Tingling_A1_34_C0_12,"No Tingling","Tingling")
  train_data$Cramps_A1_34_C0_12<-plotting_variables(train_data$Cramps_A1_34_C0_12,"No Cramps","Cramps")
  train_data$Jar_A1_34_C0_12<-plotting_variables(train_data$Jar_A1_34_C0_12,"No Jar diff","Jar diff")
  train_data$Weakness_A1_34_C0_12<-plotting_variables(train_data$Weakness_A1_34_C0_12,"No Weakness","Weakness")
  train_data$Maximal_A1_34_C0_12<-plotting_variables(train_data$Maximal_A1_34_C0_12,"No Maximal symp","Max symp")


  ##standardize the variables
  train_data_rescale <- train_data[3:8] %>%
    mutate_if(is.numeric, funs(as.numeric(scale(.))))
  train_data_rescale<- cbind(train_data[1:2],train_data_rescale,train_data[9:15])
  print(head(train_data_rescale))

  ##summary statistics
  print(ggplot(train_data_rescale, aes(x = Tingling_A1_34_C0_12, fill = Numbness_A1_34_C0_12)) +
    geom_bar(position = "fill") +
    theme_classic())

  ### box plot
  print(ggplot(train_data_rescale, aes(x = Jar_A1_34_C0_12, y = Age)) +
    geom_boxplot() +
    stat_summary(fun.y = mean,
                 geom = "point",
                 size = 3,
                 color = "steelblue") +
    theme_classic())

  ### Plot distribution maximal toxicity level by BMI
  print(ggplot(train_data_rescale, aes(x = BMI_survey)) +
    geom_density(aes(color = Maximal_A1_34_C0_12), alpha = 0.5) +
    theme_classic())

  ##anova-test for diff between groups
  anova <- aov(BMI_survey~Diabetes, train_data_rescale)
  print(summary(anova))

  #Non-linearity
  library(ggplot2)
  print(ggplot(train_data_rescale, aes(x = Age, y = BMI_survey,color = Diabetes)) +
    geom_point(size = 0.5) +
    stat_smooth(method = 'lm',formula = y~poly(x, 2),se = TRUE) +
    theme_classic())

  #correlation
  library(GGally)
  train_data_rescale <- cbind(train_data_rescale[1:8],training_phenotype[11],training_phenotype[27],training_phenotype[23],
                      training_phenotype[31],training_phenotype[39],training_phenotype[43],training_phenotype[47])

  # Convert data to numeric
  corr <- data.frame(lapply(train_data_rescale, as.integer)) #??
  # Plot the graph
  print(ggcorr(train_data_rescale[3:15],  #corr[3:15]
              method = c("pairwise", "spearman"),
              nbreaks = 6,
              hjust = 0.8,
              label = TRUE,
              label_size = 3,
              color = "grey50"))

}


check_variables(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype)

```

Check for NAs and significant correlations needed to take under
consideration.

```{r}
 #train
# View(training_genotypes) #snvs
# View(training_phenotype) #phenotypes
 training_data<-training_phenotype[1:17]
 training_data$test<-(rep(0,237))
# #View(training_data)

  #test
# View(testing_genotypes) #snvs
# View(testing_phenotype) #phenotypes
 testing_data<-testing_phenotype[1:17]
 testing_data$test<-(rep(1,100))
 
 high_symp<-c(3,4)
 sum(testing_phenotype$Tingling_of_toes_or_feet==4)
 
 phenotype_data<-rbind(training_data,testing_data)
# #View(phenotype_data)
# 
 ###Check missing data
print(paste("Dim of data:",dim(phenotype_data)))

print(paste("NA in Diabetes:",NA %in% phenotype_data$Diabetes))
print(paste("NA in BMI:",NA %in% phenotype_data$BMI_survey))
print(paste("NA in Age:",NA %in% phenotype_data$Age))
print(paste("NA in Taxane_type:",NA %in% phenotype_data$Taxan_type))
print(paste("NA in Docetaxel_mgm2:",NA %in% phenotype_data$Docetaxel_mgm2))
print(paste("NA in Paclitaxel_mgm2:",NA %in% phenotype_data$Paclitaxel_mgm2))
print(paste("NA in Time_months:",NA %in% phenotype_data$Time_months))
print(paste("NA in Tingling_of_toes_or_feet:","-9" %in% phenotype_data$Tingling_of_toes_or_feet))
print(paste("NA in Numbness_of_toes_or_feet:","-9" %in% phenotype_data$Numbness_of_toes_or_feet))
print(paste("NA in Cramps_in_feet:","-9" %in% phenotype_data$Cramps_in_feet))
print(paste("NA in Difficulty_opening_a_jar:","-9" %in% phenotype_data$Difficulty_opening_a_jar))
print(paste("NA in Difficulty_getting_up_weakness_legs:","-9" %in% phenotype_data$Difficulty_getting_up_weakness_legs))

#dummy variables
phenotype_data$Tax_type1 <- ifelse(phenotype_data$Taxan_type == "1", 1, 0)
phenotype_data$Tax_type2 <- ifelse(phenotype_data$Taxan_type == "2", 1, 0)
phenotype_data$Tax_type3 <- ifelse(phenotype_data$Taxan_type == "3", 1, 0)

phenotype_data$high_ting <- ifelse(phenotype_data$Tingling_of_toes_or_feet == "3" | phenotype_data$Tingling_of_toes_or_feet == "4" , 1, 0)
phenotype_data$high_numb <- ifelse(phenotype_data$Numbness_of_toes_or_feet == "3" | phenotype_data$Numbness_of_toes_or_feet == "4" , 1, 0)
phenotype_data$high_cramps <- ifelse(phenotype_data$Cramps_in_feet == "3" | phenotype_data$Cramps_in_feet == "4" , 1, 0)
phenotype_data$high_jar <- ifelse(phenotype_data$Difficulty_opening_a_jar == "3" | phenotype_data$Difficulty_opening_a_jar == "4" , 1, 0)
phenotype_data$high_weak <- ifelse(phenotype_data$Difficulty_getting_up_weakness_legs == "3" | phenotype_data$Difficulty_getting_up_weakness_legs == "4" , 1, 0)

# #View(phenotype_data)
# 
library(tidyverse)
library(ggpubr)
library(rstatix)

#Weltch t-test 
stat.test <- rbind(phenotype_data %>%
    t_test(Age ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Tax_type1 ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Tax_type2 ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Tax_type3 ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Paclitaxel_mgm2 ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Docetaxel_mgm2 ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(BMI_survey ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Diabetes ~ test) %>%
    add_significance(),
  phenotype_data %>%
    t_test(Time_months ~ test) %>%
    add_significance())

stat.test


#_______________
library(foreign)
library("readxl")
library(dplyr)

#Tukey test /ANOVA
phenotype_aov <- aov(test~Age+Tax_type1+Tax_type2+Tax_type3+Paclitaxel_mgm2+Docetaxel_mgm2+BMI_survey+Diabetes+Time_months+
                       high_numb+high_ting+high_cramps+high_jar+high_weak, data = phenotype_data)
print(summary(phenotype_aov, type=3,correction = c("auto", "GG", "HF", "none"),intercept=TRUE))


```

## Similarities in the data being removed by the modeling

Function preparing the sets with the correct phenotype and SNVs etc in
training and test when running a model. Some of the literature sets will
arise similarities of SNVs, not affecting the model when being included.
They can be found in the following function:

```{r}

#############training and testing sets
make_training_and_testing_set <- function(training_genotypes,training_phenotype, 
                                          testing_genotypes, testing_phenotype,
                                          phenotype, variants_to_test) {
  variants_to_test <- unique(variants_to_test)
  #print(length(variants_to_test))
  X_training <- training_genotypes[, variants_to_test] #select snvs to use
  Y_training <- training_phenotype[, phenotype] #select phenotype
  
  X_testing <- testing_genotypes[, variants_to_test] #select snvs to use
  Y_testing <- testing_phenotype[, phenotype] #select phenotype
  
  #column (variant) without variation? In that case exclude from analysis
  exclude_non_varying_variants <- c()
  for(i in 1:dim(X_training)[2]) {
    if(length(unique(X_training[, i])) == 1) {
      exclude_non_varying_variants <- c(exclude_non_varying_variants, i)
    }
  }
  if(length(exclude_non_varying_variants) > 0) {
    print(paste(length(exclude_non_varying_variants), "Variants needed exclusion from the total set of", dim(X_training)[2], "variants."))
    X_training <- X_training[, -exclude_non_varying_variants]
    X_testing <- X_testing[, -exclude_non_varying_variants]
    print(paste(dim(X_training)[2], "variants remain."))
    
  } else {
    print("No variants needed exclusion")
  }
  
  variant_names <- colnames(X_training) #due to some problems with algorithms not being able to funtion with : in names : is changed to p for position
  
  for(i in 1:length(variant_names)) {
    variant_names[i] <- sub(':', 'p', variant_names[i])
  }
  
  colnames(X_training) <- variant_names
  colnames(X_testing) <- variant_names
  
  return(list(X_training, Y_training, X_testing, Y_testing))
}

plotting_colors<- function(phenotype) {
  phenotype[which(phenotype == -9)] <- "Missing"
  phenotype[which(phenotype == 0)] <- "Low toxicity"
  phenotype[which(phenotype  == 1)] <- "High toxicity"
  
  return(phenotype)
}


extract_SNVs_model <- function(training_genotypes,training_phenotype, 
                               testing_genotypes, testing_phenotype,
                               phenotype_text, variants_to_test, variants_text, covariate) {
  
  data <- make_training_and_testing_set(training_genotypes,training_phenotype, 
                                        testing_genotypes, testing_phenotype,
                                        phenotype_text, variants_to_test)
  #data <- make_training_and_testing_set(training_genotypes,training_phenotype, 
  #                                      testing_genotypes, testing_phenotype,
  #                                      "Numbness_A1_34_C0_12", numb_snvs$SNP)
  
  
  X_training <- data[[1]] #snv
  Y_training <- data[[2]] #phenotype


  if (covariate==TRUE){
    train_data_rescale <- training_phenotype[5:11] %>%
      mutate_if(is.numeric, funs(as.numeric(scale(.))))
    #head(train_data_rescale)
    
    #combine x and y for log reg models
    Training <- as.data.frame(cbind(Y_training,X_training,
                                    train_data_rescale$Age,train_data_rescale$Diabetes,
                                    train_data_rescale$BMI_survey,train_data_rescale$Taxan_type))
    colnames(Training) <- c("Y",colnames(X_training),"Age", "Diabetes","BMI", "Taxane_type")
  } else {
    #combine x and y for log reg models
    Training <- as.data.frame(cbind(Y_training, X_training))
    colnames(Training) <- c("Y", colnames(X_training))
  }
  
  
  NA.SNV<-c()
  #fit log reg model
  if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
    training_exclude <- which(Y_training == -9)
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    #print(summary(fit))
  } else { #om inga patienter har missing fenotyp
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training, family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    #print(summary(fit))

  }
  SNVs<-c(names(summary(fit)$coefficients[,1]))
  
  #print(summary(fit)$coefficients[,1])
  paste("Model based on: ",length(summary(fit)$coefficients[,1]), " SNVs")
  return(SNVs)
}

```

Find all sets where we find similarities in the set. Here we also remove
the similar SNVs:

### Taxane pathway

```{r}
####ReArrange Taxane pathway to numbness----
new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Numbness_A1_34_C0_12",taxane_pathway_genes_and_common_variants_cpdb$SNP, "71SNV-tests CPDB all", covariate=TRUE)
new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]  
print(paste(length(taxane_pathway_genes_and_common_variants_cpdb$SNP)))
length(new.list.variables)
new.list.variables<-gsub("p", ":", new.list.variables)
taxane_pathway_genes_and_common_variants_cpdb<-taxane_pathway_genes_and_common_variants_cpdb[c(taxane_pathway_genes_and_common_variants_cpdb$SNP %in% new.list.variables), c(1:4)]
length(taxane_pathway_genes_and_common_variants_cpdb$SNP)
print(taxane_pathway_genes_and_common_variants_cpdb)
```

### Numbness in feet

```{r}
####ReArrange NUMBNESS----

sets_1<-list(numb_cpdb,numb_120cpdb,numb_210cpdb,numb_400cpdb,numb_signreview,numb_allreview)
sets_2 <-list()
for (set in sets_1) {
  new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Numbness_A1_34_C0_12", set$SNP, "", covariate=TRUE)
  new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]
  print(length(set$SNP))
  print(length(new.list.variables))
  new.list.variables<-gsub("p", ":", new.list.variables)
  set<-set %>%
    filter(set$SNP %in% new.list.variables)
  print(length(set$SNP))
  print(set)
  sets_2<-append(sets_2,list(set))
}

numb_cpdb<-sets_2[[1]]
numb_120cpdb<-sets_2[[2]]
numb_210cpdb<-sets_2[[3]]
numb_400cpdb<-sets_2[[4]]
numb_signreview<-sets_2[[5]]
numb_allreview<-sets_2[[6]]

```

### Tingling in feet

```{r}
####ReArrange Tingling----

sets_1<-list(ting_cpdb,ting_120cpdb,ting_210cpdb,ting_400cpdb,ting_signreview,ting_allreview)
sets_2 <-list()
for (set in sets_1) {
  new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Tingling_A1_34_C0_12", set$SNP, "", covariate=TRUE)
  new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]
  print(length(set$SNP))
  print(length(new.list.variables))
  new.list.variables<-gsub("p", ":", new.list.variables)
  set<-set %>%
    filter(set$SNP %in% new.list.variables)
  print(length(set$SNP))
  print(set)
  sets_2<-append(sets_2,list(set))
}

ting_cpdb<-sets_2[[1]]
ting_120cpdb<-sets_2[[2]]
ting_210cpdb<-sets_2[[3]]
ting_400cpdb<-sets_2[[4]]
ting_signreview<-sets_2[[5]]
ting_allreview<-sets_2[[6]]

```

### Cramps in feet

```{r}
####ReArrange Cramps----

sets_1<-list(cramps_cpdb,cramps_50cpdb,cramps_120cpdb,cramps_210cpdb,cramps_450cpdb,cramps_signreview,cramps_allreview)
sets_2 <-list()
for (set in sets_1) {
  new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Cramps_A1_34_C0_12", set$SNP, "", covariate=TRUE)
  new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]
  print(length(set$SNP))
  print(length(new.list.variables))
  new.list.variables<-gsub("p", ":", new.list.variables)
  set<-set %>%
    filter(set$SNP %in% new.list.variables)
  print(length(set$SNP))
  print(set)
  sets_2<-append(sets_2,list(set))
}

cramps_cpdb<-sets_2[[1]]
cramps_50cpdb <- sets_2[[2]]
cramps_120cpdb<-sets_2[[3]]
cramps_210cpdb<-sets_2[[4]]
cramps_450cpdb<-sets_2[[5]]
cramps_signreview<-sets_2[[6]]
cramps_allreview<-sets_2[[7]]

```

### Difficulty opening a jar

```{r}
####ReArrange JAR----

sets_1<-list(jar_cpdb,jar_50cpdb,jar_120cpdb,jar_210cpdb,jar_450cpdb,jar_signreview,jar_allreview)
sets_2 <-list()
for (set in sets_1) {
  new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Jar_A1_34_C0_12", set$SNP, "", covariate=TRUE)
  new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]
  print(length(set$SNP))
  print(length(new.list.variables))
  new.list.variables<-gsub("p", ":", new.list.variables)
  set<-set %>%
    filter(set$SNP %in% new.list.variables)
  print(length(set$SNP))
  print(set)
  sets_2<-append(sets_2,list(set))
}

jar_cpdb<-sets_2[[1]]
jar_50cpdb <- sets_2[[2]]
jar_120cpdb<-sets_2[[3]]
jar_210cpdb<-sets_2[[4]]
jar_450cpdb<-sets_2[[5]]
jar_signreview<-sets_2[[6]]
jar_allreview<-sets_2[[7]]

```

### Difficulty/weakness climbing stairs

```{r}
####ReArrange WEAKNESS----

sets_1<-list(weak_cpdb,weak_120cpdb,weak_210cpdb,weak_450cpdb,weak_signreview,weak_allreview)
sets_2 <-list()
for (set in sets_1) {
  new.list.variables <- extract_SNVs_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                         "Weakness_A1_34_C0_12", set$SNP, "", covariate=TRUE)
  new.list.variables<-new.list.variables[! new.list.variables %in% c("(Intercept)","Age", "Diabetes", "BMI", "Taxane_type")]
  print(length(set$SNP))
  print(length(new.list.variables))
  new.list.variables<-gsub("p", ":", new.list.variables)
  set<-set %>%
    filter(set$SNP %in% new.list.variables)
  print(length(set$SNP))
  print(set)
  sets_2<-append(sets_2,list(set))
}

weak_cpdb<-sets_2[[1]]
weak_120cpdb<-sets_2[[2]]
weak_210cpdb<-sets_2[[3]]
weak_400cpdb<-sets_2[[4]]
weak_signreview<-sets_2[[5]]
weak_allreview<-sets_2[[6]]

```

## Outliers?

A function for visual settings and other settings for colors of the
toxicity groups. PCA plot function:

```{r}
#Color setting----
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#colors (high,low,missing)
my_groups <- c("Low toxicity","High toxicity","Missing")
my_colors <- c("cyan","darkcyan","lightgrey")
my_dark_colors<-sapply(my_colors, darken)
my_cols<-c("Low toxicity"=as.character(my_colors[[1]]),"High toxicity"=as.character(my_colors[[2]]),"Missing"=as.character(my_colors[[3]]))
my_dark_cols<-c("Low toxicity"=as.character(my_dark_colors[[1]]),"High toxicity"=as.character(my_dark_colors[[2]]),"Missing"=as.character(my_dark_colors[[3]]))

tox_groups <- c("Low toxicity","High toxicity") #,"Missing")
tox_colors <- c("cyan","darkcyan") #,"lightgrey")
tox_dark_colors<-c("darkslategrey","black")#sapply(tox_colors, darken)
tox_cols<-c("Low toxicity"=as.character(tox_colors[[1]]),"High toxicity"=as.character(tox_colors[[2]])) #,"Missing"=as.character(tox_colors[[3]]))
tox_dark_cols<-c("Low toxicity"=as.character(tox_dark_colors[[1]]),"High toxicity"=as.character(tox_dark_colors[[2]])) #,"Missing"=as.character(tox_dark_colors[[3]]))

tox_test_groups <- c("Low toxicity","High toxicity") #,"Missing")
tox_test_colors <- c("mediumpurple1","mediumpurple4") #,"lightgrey")
tox_test_dark_colors<-c("midnightblue","black") #sapply(tox_test_colors, darken)
tox_test_cols<-c("Low toxicity"=as.character(tox_test_colors[[1]]),"High toxicity"=as.character(tox_test_colors[[2]])) #,"Missing"=as.character(tox_colors[[3]]))
tox_test_dark_cols<-c("Low toxicity"=as.character(tox_test_dark_colors[[1]]),"High toxicity"=as.character(tox_test_dark_colors[[2]])) #,"Missing"=as.character(tox_dark_colors[[3]]))


type_groups <- c("train","test") #,"Missing")
type_colors <- c("skyblue","cornflowerblue") #,"lightgrey")
type_dark_colors<-sapply(type_colors, darken)
type_cols<-c("train"=as.character(type_colors[[1]]),"test"=as.character(type_colors[[2]])) #,"Missing"=as.character(type_colors[[3]]))
type_dark_cols<-c("train"=as.character(type_dark_colors[[1]]),"test"=as.character(type_dark_colors[[2]])) #,"Missing"=as.character(type_dark_colors[[3]]))


PCA_and_plot <- function(training_genotypes,training_phenotype, 
                         testing_genotypes, testing_phenotype,
                         phenotype, variants_to_test, variants) {
  
  
  
  data <- make_training_and_testing_set(training_genotypes,training_phenotype, 
                                        testing_genotypes, testing_phenotype,
                                        phenotype, variants_to_test)
 
  X_training <- data[[1]]
  Y_training <- data[[2]]
  X_testing <- data[[3]]
  Y_testing <- data[[4]]
  
  #plotting colors
  training_plotting_color <- plotting_colors(Y_training)
  testing_plotting_color <- plotting_colors(Y_testing)
  
  
  #PCA
  pca_prcomp <- prcomp(X_training, scale. = T, center = T)
  summary(pca_prcomp)
  var_explained <- pca_prcomp$sdev^2/sum(pca_prcomp$sdev^2)
  
  pcatab <- data.frame(sample.id = samples_training,
                       EV1 = pca_prcomp$x[,1],    # the first eigenvector
                       EV2 = pca_prcomp$x[,2],    # the second eigenvector
                       stringsAsFactors = FALSE)
  
  pca_training_plot <- ggplot(data=pcatab, aes(x=EV1, y=EV2, fill=training_plotting_color, key = training_phenotype$Lopnummer)) + #, alpha=0.5) +
    geom_point(size = 2, colour="black", stroke = 0.2, pch=21) + 
    labs(title = paste("PCA: For", phenotype, "\nusing", variants),
         x = paste("PC 1: ", round(var_explained[1]*100,digits = 2), "%", sep=""),
         y = paste("PC 2: ", round(var_explained[2]*100,digits = 2), "%", sep=""), fill = phenotype) +
    scale_fill_manual(values=my_cols) + theme_cowplot(12) + theme(axis.text=element_text(size=10),
                                                                                          axis.title=element_text(size=10),
                                                                                          plot.title=element_text(size=10,face="bold"),
                                                                                          legend.title=element_text(size=10,face="bold"),
                                                                                          legend.text=element_text(size=10), 
                                                                                          axis.text.x = element_text(colour = "black"),
                                                                                          axis.text.y = element_text(colour = "black"))
  
  
  ggplotly(pca_training_plot)
  #if needed here add the prediction of the testing PCA points based on the pca above and plot them as well (not needed for initial analyses)
}

PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Numbness_A1_34_C0_12", numb_snvs$SNP, "SNV-tests")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Numbness_A1_34_C0_12", numb_cpdb$SNP, "SNV-tests from CPDB")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
              "Numbness_A1_34_C0_12", taxane_pathway_genes_and_variants$variant, "Taxane pathway")

PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
              "Tingling_A1_34_C0_12", ting_snvs$SNP, "SNV-tests")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Tingling_A1_34_C0_12", ting_cpdb$SNP, "SNV-tests from CPDB")

PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
              "Cramps_A1_34_C0_12", cramps_snvs$SNP, "SNV-tests")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Cramps_A1_34_C0_12", cramps_cpdb$SNP, "SNV-tests from CPDB")

PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
              "Jar_A1_34_C0_12", jar_snvs$SNP, "SNV-tests")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Jar_A1_34_C0_12", jar_cpdb$SNP, "SNV-tests from CPDB")

PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
            "Weakness_A1_34_C0_12", weak_snvs$SNP, "SNV-tests")
PCA_and_plot(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
             "Weakness_A1_34_C0_12", weak_cpdb$SNP, "SNV-tests from CPDB")



```
