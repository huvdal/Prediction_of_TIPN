---
editor_options: 
  markdown: 
    wrap: 72
---

# Build and test logistic regression models for prediction of TIPN

The following is a Rmarkdown tutorial for running similar models found
in the article:

Engvall, K., Uvdal, H., Björn, N., Åvall-Lundqvist, E., & Gréen, H.
(2024). Prediction models of persistent taxane-induced peripheral
neuropathy among breast cancer survivors using whole-exome sequencing.
NPJ Precision Oncology, 8(1). https://doi.org/10.1038/s41698-024-00594-x

## Packages:

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

library(umap)
library(ggbeeswarm)
library(tidyverse)

library(ISLR) #contains Hitters dataset
library(rpart) #for fitting decision trees
library(rpart.plot)
library(tree)
library(caret)

```

# Build logistic regression models

In previous scripts we have loaded in the cohort and literature data as
well as formatted and cleaned them up for building and testing the
logistic regression model. Below we create a few functions for
optimizing and visualizing the building and later validation of the
models.

## Functions for the prediction plot shown as a violin plot

Below as some functions how to visualize the prediction plot. Either in
one violin per train/test set or in splitted violins.

```{r}
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
```

## AUC ROC curve function

Function for the AUC plot for the model building and validation.

```{r}
auc_plotting <- function(preds, truth, phenotype_text, variants_text) {
  evaluate_preds <- prediction(preds, truth)
  # calculate probabilities for TPR/FPR for predictions
  perf <- performance(evaluate_preds,"tpr","fpr")
  auc <- performance(evaluate_preds,"auc")@y.values[[1]] # shows calculated AUC for model
  
  p <- ggplot(data = as.data.frame(cbind(perf@x.values[[1]], perf@y.values[[1]]))) +
    geom_line(aes(x=V1, y=V2)) + 
    geom_segment(aes(x=0, y=0, xend=1, yend=1), linetype = "dashed", color="black") +  
    labs(x = "False positive rate", y = "True positive rate")  +
    coord_cartesian(ylim = c(0, 1))  + 
    ggtitle("ROC-curve") +
    annotate("text", x=0.75, y=0.25, label = paste("AUC = ", round(auc, 4) * 100, "%", sep = ""), size = 8/.pt) + 
    #ggtitle(paste("ROC-curve: For", phenotype_text, "training using", variants_text, "\nAUC = ", round(auc, 3) * 100, "%", sep ="")) +
    theme_cowplot(12) + theme(axis.text=element_text(size=10),
                              axis.title=element_text(size=10),
                              plot.title=element_text(size=10,face="bold"),
                              legend.title=element_text(size=10,face="bold"),
                              legend.text=element_text(size=10), 
                              axis.text.x = element_text(colour = "black"),
                              axis.text.y = element_text(colour = "black"))
  
  return(p)
}
```

## Accuracy curve function

Function for the accuracy plot for the model building and validation.
This plot make it possible to see when we have the highest accuracy.

```{r}
accuracy_plot <- function(preds, truth, phenotype_text, variants_text, cutoff) {
  acc.perf = performance(prediction(preds, truth), measure = "acc")
  data <- as.data.frame(cbind(acc.perf@x.values[[1]], acc.perf@y.values[[1]]))
  #data <- data[-which(data$V1 > 1), ] #tabort det som ?r med cutoff st?rre ?n 1 because it does not make sense
  accuracy.pos<-which(abs(round(acc.perf@x.values[[1]],4)-cutoff)== min(abs(round(acc.perf@x.values[[1]],4)-cutoff)))
  #accuracy.pos_lo<-which(abs(round(acc.perf@x.values[[1]],4)-cutoff_low)== min(abs(round(acc.perf@x.values[[1]],4)-cutoff_low)))
  p <- ggplot(data = data) + #as.data.frame(cbind(acc.perf@x.values[[1]], acc.perf@y.values[[1]]))) +
    geom_line(aes(x=V1, y=V2), color="black") +
    labs(x = "Cutoff", y = "Accuracy")  +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) + 
    #geom_vline(xintercept=cutoff_low, linetype = "dashed", color="black") +
    geom_vline(xintercept=cutoff, linetype = "dashed", color="black") +
    ggtitle("Accuracy curve") +
    annotate("text", x=0.5, y=0.15, label = paste("Max accuracy = ", round(acc.perf@y.values[[1]][accuracy.pos], 4) * 100, "%",
                                                   "\nCutoff = ", cutoff, #round(acc.perf@x.values[[1]][which.max(acc.perf@y.values[[1]])]
                                                   "\nSensitivity = ", round(performance(prediction(preds, truth), measure = "sens")@y.values[[1]][accuracy.pos], 4) * 100, "%",
                                                   "\nSpecificity = ", round(performance(prediction(preds, truth), measure = "spec")@y.values[[1]][accuracy.pos], 4) * 100, "%",
                                                   sep = ""), size = 8/.pt) +
    theme_cowplot(12) + theme(axis.text=element_text(size=10),
                              axis.title=element_text(size=10),
                              plot.title=element_text(size=10,face="bold"),
                              legend.title=element_text(size=10,face="bold"),
                              legend.text=element_text(size=10), 
                              axis.text.x = element_text(colour = "black"),
                              axis.text.y = element_text(colour = "black"))
  accuracies<-order(acc.perf@y.values[[1]], decreasing=TRUE)
  #print(accuracies) #1 (0.99972), 18 (0.82447), 23(0.64527) or combo 30(0.5575)
  
  #Find specific accuracies?
  # for (set_accuracy in accuracies[c(1,18,23)]){
  #   cutoff<-round(acc.perf@x.values[[1]][set_accuracy], 5)
  #   print(paste("Max accuracy = ", acc.perf@y.values[[1]][set_accuracy] * 100, "%",
  #               "Cutoff = ", cutoff,
  #               "Sensitivity = ", round(performance(prediction(preds, truth), measure = "sens")@y.values[[1]][set_accuracy], 4) * 100, "%",
  #               "Specificity = ", round(performance(prediction(preds, truth), measure = "spec")@y.values[[1]][set_accuracy], 4) * 100, "%"))
  # }
  
  return(p)
}
```

## Cut-off functions for different accuracies

Functions for cut-offs for different accuracies based on the accuracy
plot above. These are also involved in the prediction plot to see where
the cut-off will appear in the predictions.

```{r}

#selects the cut-off based on the highest accuracy in the whole model prediction
best_cutoff <- function(preds, truth) {
  acc.perf = performance(prediction(preds, truth), measure = "acc")
  return(round(acc.perf@x.values[[1]][which.max(acc.perf@y.values[[1]])], 5))
}

#selects the cut-off based on an accuracy input
high_cutoff <- function(preds, truth) {
  acc.perf = performance(prediction(preds, truth), measure = "acc")
  acc.cut<-max(acc.perf@x.values[[1]][which(round(acc.perf@y.values[[1]],4)==0.8128)]) 

  return(round(acc.cut,5))
}

low_cutoff <- function(preds, truth) {
  acc.perf = performance(prediction(preds, truth), measure = "acc")
  acc.cut<-min(acc.perf@x.values[[1]][which(round(acc.perf@y.values[[1]],3)==0.753)])#A1 0.753 C2numb 0.834 C2ting 0.804 eller 0.783

  
  return(round(acc.cut,5))
  
}

```

## Prediction plot

Function for building the prediction model as a violin plot
investigating the prediction distribution in the reported patient
groups. Cut-off is included to investigate how much of the correctly
predicted patients will be included corresponding to the accuracy.

```{r}
prediction_plot <- function(preds, truth, phenotype_text, variants_text, title, cutoff) { #train_test_factor
  
  x <- as.data.frame(cbind(preds, truth))
  
  colnames(x) <- c("probablility", "tox_level")
  x$tox_level <- plotting_colors(x$tox_level)
  #x<-x[!x$tox_level == "Missing",]
  
  if(sum(x$tox_level == "Missing") > 0) { #om n?gon har missing fenotype
    
    exclude <- which(x$tox_level == "Missing")
    
    #x <- cbind(x, train_test_factor)
    p <- ggplot(data = x,aes(x=tox_level,y=probablility, group = tox_level, fill=tox_level)) +
      geom_violin(scale="width", draw_quantiles = c(0.7, 0.28),linetype = "dashed") +
      geom_violin(scale="width",fill="transparent",draw_quantiles = 0.5) +
      #geom_boxplot(outlier.shape = NA) + #shape=train_test_factor),
      geom_sina(aes(colour = tox_level, group = tox_level), alpha=1,scale=F, 
                method="density", 
                maxwidth = .6) +  
      #scale_shape_manual(values=c(4, 1), name  ="Type of data: ",
      #                   labels=c("Testing", "Training")) +
      #geom_hline(yintercept=cutoff_low, linetype = "dashed", color="black", size=1) +
      geom_hline(yintercept=cutoff, linetype = "dashed", color="black", size=1) +
      
      scale_x_discrete(limits=tox_groups) +
      scale_color_manual(values=tox_dark_cols) +
      scale_fill_manual(values=tox_cols, name = "Toxicity: ") +
      labs(x = NULL, y = "Probability of toxicity")  +
      coord_cartesian(ylim = c(0, 1)) + ggtitle(title) +
      # ggtitle(paste("Logistic regression prediction model:\nFor", phenotype_text, "training using", variants_text)) + 
      theme_cowplot(12) +
      theme(legend.position = "none") +
      theme(axis.text=element_text(size=12),                                                                                      
            axis.title=element_text(size=12),
            plot.title=element_text(size=12,face="bold"),
            legend.title=element_text(size=12,face="bold"),
            legend.text=element_text(size=12), 
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"))
    
  } else { #omingen fenotype ?r missing
    #x <- cbind(x, train_test_factor)
    p <- ggplot(data = x, aes(x=tox_level,y=probablility, group = tox_level, fill=tox_level)) +
      geom_violin(scale="width", draw_quantiles = c(0.7, 0.28), linetype = "dashed") +
      geom_violin(scale="width",fill="transparent",draw_quantiles = 0.5) +
      #geom_boxplot(outlier.shape = NA, aes(x=tox_level,y=probablility, group = tox_level, fill=tox_level)) + #shape=train_test_factor), 
      geom_sina(aes(colour = tox_level, group = tox_level), alpha=1,scale=F, 
                method="censity", 
                maxwidth = .6) + 
      #scale_shape_manual(values=c(4, 1), name  ="Type of data: ",
      #                   labels=c("Testing", "Training")) +
      #geom_hline(yintercept=cutoff_low, linetype = "dashed", color="black", size=1) +
      geom_hline(yintercept=cutoff, linetype = "dashed", color="black", size=1) +
      
      scale_x_discrete(limits=tox_groups) +
      scale_color_manual(values=tox_dark_cols) +
      scale_fill_manual(values=tox_cols, name = "Toxicity: ") +
      labs(x = NULL, y = "Probability of toxicity")  +
      coord_cartesian(ylim = c(0, 1)) + ggtitle(title) +
      # ggtitle(paste("Logistic regression prediction model:\nFor", phenotype_text, "training using", variants_text)) + 
      theme_cowplot(12) +
      theme(legend.position = "none") +
      theme(axis.text=element_text(size=12),                                                                                      
            axis.title=element_text(size=12),
            plot.title=element_text(size=12, face="bold"),
            legend.title=element_text(size=12, face="bold"),
            legend.text=element_text(size=12), 
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"))
  }
  return(p)
}
```

# Final model functions

Final function calling on all the previous functions to build and/or
test the model.

The first function below only build the model in the training data.

```{r}
make_log_reg_model <- function(training_genotypes,training_phenotype, 
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

  
  if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
    training_exclude <- which(Y_training == -9)
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    
    roc_auc_plot <- auc_plotting(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text, cutoff[[1]])
    print(table(preds[-training_exclude] >= best_cutoff(preds[-training_exclude], Y_training[-training_exclude]), Y_training[-training_exclude] == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training[-training_exclude] == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff[[1]])
  } else { #om inga patienter har missing fenotyp
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training, family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff<- best_cutoff(preds, Y_training)
    
    roc_auc_plot <- auc_plotting(preds, Y_training, phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds, Y_training, phenotype_text, variants_text, cutoff[[1]])
    print(table(preds >= best_cutoff(preds, Y_training), Y_training == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff[[1]])
  }
  if (covariate==TRUE) {
    p <- ggarrange(log_reg_pred_plot, accuracy_plot, roc_auc_plot, ncol = 3, nrow = 1,labels = c("A", "B", "C"))
    annotate_figure(p, top = text_grob(paste("Training: For", phenotype_text, "using", variants_text, " with covariates")))
  } else {
    p <- ggarrange(log_reg_pred_plot, accuracy_plot, roc_auc_plot, ncol = 3, nrow = 1,labels = c("A", "B", "C"))
    annotate_figure(p, top = text_grob(paste("Training: For", phenotype_text, "using", variants_text)))
  }
  
}
```

## Comparison of prediction models

A comparative function regarding prediction distribution, accuracy and
AUC ROC. If one want to combine the prediction plot a function below can
be found.

```{r}
prediction_joint_plot <- function(preds_train, Y_training, preds_test, Y_testing, phenotype_text, variants_text, covariate) {
  
  y_train <- as.data.frame(cbind(preds_train, Y_training))
  y_test <-as.data.frame(cbind(preds_test, Y_testing))
  
  colnames(y_train) <- c("probablility", "tox_level")
  y_train$tox_level <- plotting_colors(y_train$tox_level)
  y_train$type<- "train"
  y_train<-y_train[!y_train$tox_level=="Missing",]
  
  colnames(y_test) <- c("probablility", "tox_level")
  y_test$tox_level <- plotting_colors(y_test$tox_level)
  y_test$type<- "test"
  y_test<-y_test[!y_test$tox_level=="Missing",]
  
  
  my_data<- rbind(y_train,y_test)
  my_data <- my_data %>%
    mutate(type=factor(type,levels=c("train","test"))) %>%
    mutate(tox_level=factor(tox_level,levels=c("Low toxicity","High toxicity")))
    
  #violin split
  if (covariate==TRUE) {
    p <- ggplot(my_data, aes(x=tox_level,y=probablility, color=type,fill=type)) + #,fill=tox_level  
      geom_violin(scale="width", draw_quantiles = c(0.5)) +
      scale_color_manual(values = type_dark_cols) +
      geom_sina() + 
      scale_fill_manual(values = type_cols) + #tox_cols
      labs(x = NULL, y = "Probability of toxicity")  +
      coord_cartesian(ylim = c(0, 1)) + ggtitle(paste("Train & Test: For", phenotype_text, "using", variants_text, " with covariates")) 
  }else{
    p <- ggplot(my_data, aes(x=tox_level,y=probablility, color=type,fill=type)) + #,fill=tox_level  
      geom_violin(scale="width", draw_quantiles = c(0.5)) +
      scale_color_manual(values = type_dark_cols) +
      geom_sina() + 
      scale_fill_manual(values = type_cols) + #tox_cols
      labs(x = NULL, y = "Probability of toxicity")  +
      coord_cartesian(ylim = c(0, 1)) + ggtitle(paste("Train & Test: For", phenotype_text, "using", variants_text)) 
  }
  p<-p + theme_cowplot(12) +
    theme(axis.text=element_text(size=10),                                                                                      
          axis.title=element_text(size=10),
          plot.title=element_text(size=10,face="bold"),
          legend.title=element_text(size=10,face="bold"),
          legend.text=element_text(size=10), 
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))
  

}

make_log_reg_model_A <- function(training_genotypes,training_phenotype, 
                               testing_genotypes, testing_phenotype,
                               phenotype_text, variants_to_test, variants_text, covariate, plot) {
  
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
  
 
  if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
    training_exclude <- which(Y_training == -9)
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    
    roc_auc_plot <- auc_plotting(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text, cutoff[[1]])
    print(table(preds[-training_exclude] >= best_cutoff(preds[-training_exclude], Y_training[-training_exclude]), Y_training[-training_exclude] == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training[-training_exclude] == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff[[1]])
  } else { #om inga patienter har missing fenotyp
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training, family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff<- best_cutoff(preds, Y_training)
    
    roc_auc_plot <- auc_plotting(preds, Y_training, phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds, Y_training, phenotype_text, variants_text, cutoff[[1]])
    print(table(preds >= best_cutoff(preds, Y_training), Y_training == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff[[1]])
  }
  
  ##for joint plot
  #pred_joint_plot<-prediction_joint_plot(preds_train, Y_training, preds_test, Y_testing, phenotype_text, variants_text, covariate)
  #ggarrange(pred_joint_plot, roc_auc_plot, ncol = 2, nrow = 1)
  
  if (covariate==TRUE) {
    if (plot=="prediction") {
      p<-ggarrange(log_reg_pred_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text), bold(' with covariates')))))
    } else if (plot=="accuracy") {
      p<-ggarrange(accuracy_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text), bold(' with covariates')))))
    } else if (plot=="auc") {
      p<-ggarrange(roc_auc_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text), bold(' with covariates')))))
    }
  } else {
    if (plot=="prediction") {
      p<-ggarrange(log_reg_pred_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text)))))
    } else if (plot=="accuracy") {
      p<-ggarrange(accuracy_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text)))))
    } else if (plot=="auc") {
      p<-ggarrange(roc_auc_plot, ncol = 1, nrow = 1)
      annotate_figure(p, top = text_grob(substitute(paste(bold(variants_text)))))
    }
  }
}

```

Example of calling the function and running several types of comparision
either in prediction, accuracy or AUC ROC.

```{r}
#Prediction comparison B1 with different FDR cut-offs (p-values 0.00001-0.005)
p <- ggarrange(make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_50cpdb$SNP, "49SNV-tests", covariate=TRUE, "prediction"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_cpdb$SNP, "71SNV-tests", covariate=TRUE, "prediction"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_120cpdb$SNP, "120SNV-tests", covariate=TRUE, "prediction"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_210cpdb$SNP, "214SNV-tests", covariate=TRUE, "prediction"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_400cpdb$SNP, "398SNV-tests", covariate=TRUE, "prediction"),
              ncol = 5, nrow = 1,labels = c("A", "B", "C","D","E"))
annotate_figure(p, top = text_grob(paste("Training: Numbness using SNV-test with different FDR cut-offs (p-values 0.00001-0.005)")))
length(numb_50cpdb$SNP)
length(numb_cpdb$SNP)
length(numb_120cpdb$SNP)
length(numb_210cpdb$SNP)
length(numb_400cpdb$SNP)

#Accuracy comparison B2 with different CPDB cut-offs
p <- ggarrange(make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_cpdb$SNP, "71SNV-tests CPDB all", covariate=TRUE, "accuracy"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.5], "63SNV CPDB p<0.5", covariate=TRUE, "accuracy"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.1], "42SNV CPDB p<0.1", covariate=TRUE, "accuracy"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.05], "32SNV CPDB p<0.05", covariate=TRUE, "accuracy"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.02], "21SNV CPDB p<0.02", covariate=TRUE, "accuracy"),
              ncol = 5, nrow = 1,labels = c("A", "B", "C","D","E"))
annotate_figure(p, top = text_grob(paste("Training: Numbness using 71SNV-test with different CPDB cut-offs")))
length(numb_cpdb$SNP)
length(numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.5])
length(numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.1])
length(numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.05])
length(numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.02])


#AUC ROC comparison A1 and A2
p <- ggarrange(make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_signreview$SNP, "A1 8SNV Review p<0.05", covariate=TRUE, "auc"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", numb_allreview$SNP, "A1 26SNV Review", covariate=TRUE, "auc"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                    "Numbness_A1_34_C0_12", taxane_pathway_genes_and_common_variants_cpdb$SNP,"A2 KEGG Taxane pathway common variants", covariate=TRUE, "auc"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", neurodegeneration_genes_and_variants$variant, "A2 KEGG Neurodegeneration", covariate=TRUE, "auc"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", axon_guidance_genes_and_variants$variant, "A2 KEGG Axon guidance", covariate=TRUE, "auc"),
              make_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                   "Numbness_A1_34_C0_12", axon_regeneration_genes_and_variants$variant, "A2 KEGG Axon regeneration", covariate=TRUE, "auc"),
              ncol = 6, nrow = 1,labels = c("A", "B", "C","D","E", "F"))
annotate_figure(p, top = text_grob(paste("Training: Numbness A1 and A2 comparison")))
length(numb_signreview$SNP)
length(numb_allreview$SNP)
length(taxane_pathway_genes_and_common_variants_cpdb$SNP)
length(neurodegeneration_genes_and_variants$SNP)
length(axon_guidance_genes_and_variants$SNP)
length(axon_regeneration_genes_and_variants$SNP)
  
```

## B1

Run the training of a model based on variants and cofactors based on
SNV/INDEL association analysis using PLINK on the symptoms (B1).

```{r}
#SNV/INDEL association analysis with lowest FDR using PLINK on the symptom Numbness in feet (B1)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Numbness_A1_34_C0_12", numb_cpdb$SNP, "SNV-tests", covariate=TRUE)
length(numb_cpdb$SNP)
length(unique(numb_cpdb$GeneName_GRCh38_103_gtf))

#SNV/INDEL association analysis with lowest FDR using PLINK on the symptom Tingling in feet (B1)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Tingling_A1_34_C0_12", ting_cpdb$SNP, "79SNV-tests", covariate=TRUE)
length(ting_cpdb$SNP)
length(unique(ting_cpdb$GeneName_GRCh38_103_gtf))

#SNV/INDEL association analysis with lowest FDR using PLINK on the symptom Cramps in feet (B1)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Cramps_A1_34_C0_12", cramps_cpdb$SNP, "69SNV-tests", covariate=TRUE)
length(cramps_cpdb$SNP)
length(unique(cramps_cpdb$GeneName_GRCh38_103_gtf))

#SNV/INDEL association analysis with lowest FDR using PLINK on the symptom Difficulty opening a jar (B1)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Jar_A1_34_C0_12", jar_cpdb$SNP, "33SNV-tests", covariate=TRUE)
length(jar_cpdb$SNP)
length(unique(jar_cpdb$GeneName_GRCh38_103_gtf))

#SNV/INDEL association analysis with lowest FDR using PLINK on the symptom Difficulty climing stairs (B1)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Weakness_A1_34_C0_12", weak_cpdb$SNP, "46SNV-tests", covariate=TRUE)
length(weak_cpdb$SNP)
length(unique(weak_cpdb$GeneName_GRCh38_103_gtf))
```

## B2

Run the training of a model based on variants and cofactors based on
SNV/INDEL association analysis using PLINK with additional
gene/region-based associaltion analysis and CPDB over-representation
enrichment analysis on the symptoms (B2).

```{r}
#SNV/INDEL association analysis with lowest FDR using PLINK and CPDB over-representation enrichment analysis p<0.05 on the symptom Numbness in feet (B2)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Numbness_A1_34_C0_12", numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.005], "71SNV-tests from CPDB p<0.005", covariate=TRUE)
length(numb_cpdb$SNP[numb_cpdb$CPDB_p_value <0.005])
length(unique(numb_cpdb$GeneName_GRCh38_103_gtf[numb_cpdb$CPDB_p_value <0.005]))

#SNV/INDEL association analysis with lowest FDR using PLINK and CPDB over-representation enrichment analysis p<0.019 on the symptom Tingling in feet (B2)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Tingling_A1_34_C0_12", ting_400cpdb$SNP[ting_400cpdb$CPDB_p_value <0.019], "408SNV-tests CPDB p<0.019", covariate=TRUE)
length(ting_400cpdb$SNP[ting_400cpdb$CPDB_p_value <0.019])
length(unique(ting_400cpdb$GeneName_GRCh38_103_gtf[ting_400cpdb$CPDB_p_value <0.019]))

#SNV/INDEL association analysis with lowest FDR using PLINK and CPDB over-representation enrichment analysis p<0.2 on the symptom Cramps in feet (B2)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                 "Cramps_A1_34_C0_12", cramps_210cpdb$SNP[cramps_210cpdb$CPDB_p_value <0.2], "213SNV-tests CPDB p<0.2", covariate=TRUE)
length(cramps_210cpdb$SNP[cramps_210cpdb$CPDB_p_value <0.2])
length(unique(cramps_210cpdb$GeneName_GRCh38_103_gtf[cramps_210cpdb$CPDB_p_value <0.2]))

#SNV/INDEL association analysis with lowest FDR using PLINK and CPDB over-representation enrichment analysis p<0.16 on the symptom Difficulty opening a jar (B2)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                 "Jar_A1_34_C0_12", jar_210cpdb$SNP[jar_210cpdb$CPDB_p_value <0.16], "213SNV-tests CPDB p<0.16", covariate=TRUE)
length(jar_210cpdb$SNP[jar_210cpdb$CPDB_p_value <0.16])
length(unique(jar_210cpdb$GeneName_GRCh38_103_gtf[jar_210cpdb$CPDB_p_value <0.16]))

#SNV/INDEL association analysis with lowest FDR using PLINK and CPDB over-representation enrichment analysis p<0.075 on the symptom Difficulty climing stairs (B2)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                 "Weakness_A1_34_C0_12", weak_210cpdb$SNP[weak_210cpdb$CPDB_p_value <0.075], "213SNV-tests from CPDB <0.075", covariate=TRUE)
length(weak_210cpdb$SNP[weak_210cpdb$CPDB_p_value <0.075])
length(unique(weak_210cpdb$GeneName_GRCh38_103_gtf[weak_210cpdb$CPDB_p_value <0.075]))

  
```

## C1

In the C1 models we combine the literature data in A1 and the cohort
data in B2.

Function below to combine the sets between A1 and B2.

```{r}
#Combine A1 and B2
combine_sets <- function(review_variants, cohort_variants, cutoff) {

      combined_set<-rbind(review_variants %>%
                                          filter(!duplicated('SNP')),
                                        cohort_variants %>%
                                          filter(CPDB_p_value<cutoff & !duplicated('SNP'))) 
      length(combined_set$SNP)
      sum(combined_set$SNP %in% combined_set$SNP)
    return (combined_set)
}
```

Calling on the function and building the model investigating prediciton
distribution, accuracy and AUC ROC.

```{r}
#add nas in sign for equal df length for combination of sets
CPDB_p_value<-c(rep(NA,8))
numb_signreview_cpdb<- cbind(numb_signreview,CPDB_p_value) 
numb400_signreview_cpdb<-combine_sets(numb_signreview_cpdb,numb_400cpdb, 0.0025)

make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                   "Numbness_A1_34_C0_12", numb400_signreview_cpdb$SNP, "86SNVs review p<0.05, 400SNVtest CPDB p<0.0025", covariate=TRUE)
length(numb400_signreview_cpdb$SNP)
length(unique(numb400_signreview_cpdb$GeneName_GRCh38_103_gtf))

#One variant not significant in Tingling
ting_signreview_new <- subset(ting_signreview, !(ting_signreview$SNP %in% c("rs1045642")))
length(ting_signreview_new$SNP)
CPDB_p_value<-c(rep(NA,7))
ting_signreview_cpdb <-cbind(ting_signreview_new,CPDB_p_value)    
ting400_signreview_cpdb<-combine_sets(ting_signreview_cpdb,ting_400cpdb, 0.0019)

make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                     "Tingling_A1_34_C0_12", ting400_signreview_cpdb$SNP, "86SNVs review p<0.05, 408SNVtests CPDB p<0.019", covariate=TRUE)
  length(ting400_signreview_cpdb$SNP)
```

## C2

In the C2 models we add variable importance to the variant selection for
fine-tuning.

Below is a function to retrieve the variable importance from the
logarithmic regression building output.

```{r}
extract_LOGREG_model <- function(training_genotypes,training_phenotype, 
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
    
    
    if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
      training_exclude <- which(Y_training == -9)
      fit <- c()
      fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
      preds <- c()
      preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
      print(summary(fit))
      
      cutoff_high<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
     
      roc_auc_plot <- auc_plotting(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text)
      accuracy_plot <- accuracy_plot(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text, cutoff_high)
      print(table(preds[-training_exclude] >= best_cutoff(preds[-training_exclude], Y_training[-training_exclude]), Y_training[-training_exclude] == 1))
      #print(table(preds[-training_exclude] >= .5, Y_training[-training_exclude] == 1))
      log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff_high)
    } else { #om inga patienter har missing fenotyp
      fit <- c()
      fit <- glm(as.factor(Y)~., data = Training, family="binomial")
      preds <- c()
      preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
      print(summary(fit))
      #print(confint(fit))
      #print(coefficients(fit))
      #print(preds)
      cutoff_high<- best_cutoff(preds, Y_training)
      #cutoff_low<- low_cutoff(preds[-training_exclude], Y_training[-training_exclude])
      
      roc_auc_plot <- auc_plotting(preds, Y_training, phenotype_text, variants_text)
      accuracy_plot <- accuracy_plot(preds, Y_training, phenotype_text, variants_text, cutoff_high)
      print(table(preds >= best_cutoff(preds, Y_training), Y_training == 1))
      #print(table(preds[-training_exclude] >= .5, Y_training == 1))
      log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff_high)
    }

  variable_importance = varImp(fit, scale=FALSE)
  return(variable_importance)
}
```

Below is a function filtering the variable importance based on a cutoff
and a visualization of the variable importance in a plot.

```{r}
var_Imp <- function(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
           phenotype_text, dataframe, variants_to_test, variants_text, covariate, cutoff) {
  
  var_imp_df<- extract_LOGREG_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                       phenotype_text, variants_to_test, variants_text, covariate)
                    
  var_imp_df$labels <- factor(rownames(var_imp_df))
  var_imp_df$labels <- reorder(var_imp_df$labels, var_imp_df$Overall)
  var_imp_df$SNP <- row.names(var_imp_df)
                    
  var_imp_df<- var_imp_df %>%
    as.data.frame() %>% 
    arrange(desc(Overall)) 
                  
  var_imp_df_filt<- var_imp_df[var_imp_df$Overall>cutoff,]
  
   
    if (covariate==TRUE) {
    print(ggplot(data= var_imp_df_filt, aes(x=labels,y=Overall)) +
      geom_bar(position="dodge",stat="identity",width = 0, color = "black") + 
      coord_flip() + geom_point(color='black')+
      ggtitle(paste("Variable Importance -", phenotype_text, variants_text, " with covariates")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'))  +
      labs(x = "Variables",
           y = "Importance"))
    } else {
    print(ggplot(data= var_imp_df_filt, aes(x=labels,y=Overall)) +
      geom_bar(position="dodge",stat="identity",width = 0, color = "black") + 
      coord_flip() + geom_point(color='black')+
      ggtitle(paste("Variable Importance -", phenotype_text, variants_text)) + 
       theme(plot.title = element_text(hjust = 0.5)) +
      theme(panel.background = element_rect(fill = 'white', colour = 'black'))  +
      labs(x = "Variables",
           y = "Importance"))
    }
  
  length(var_imp_df$SNP[var_imp_df$Overall>cutoff])
  #threshold setting
  
  SNP_high_imp<-var_imp_df$SNP[var_imp_df$Overall>cutoff]
  length(SNP_high_imp)
  print(cutoff)
  variants_to_test_high_imp <- dataframe %>%
    filter(variants_to_test %in% SNP_high_imp)
  
  return(variants_to_test_high_imp)
  
}

```

Example of using the functions and retrieving the filtered model.

```{r}
numb_400cpdb_high_imp<-var_Imp(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                               "Numbness_A1_34_C0_12", numb400_signreview_cpdb, numb400_signreview_cpdb$SNP, "398SNVs, CPDB p<0.0024, review p<0.05", covariate=TRUE, 1.67)
make_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                               "Numbness_A1_34_C0_12", numb_400cpdb_high_imp$SNP, "35SNVs from FDR perm 398SNVs, CPDB p<0.0024, review p<0.05, high imp>1.67", covariate=TRUE)
            length(numb_400cpdb_high_imp$SNP)
            length(unique(numb_400cpdb_high_imp$GeneName_GRCh38_103_gtf))
            
      
```

# Validation and testing of models

Function for validating the model in the test set. Calculating rmse,
brierscore and rocit for the model as well as a comparision of the
prediction distribution, accuracy and AUC in a plot.

```{r}
test_log_reg_model <- function(training_genotypes,training_phenotype, 
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
  X_testing <- data[[3]] #snv
  Y_testing <- data[[4]] #phenotype
  
  ##checking continous variables in set
  #continuous <-select_if(training_phenotype, is.numeric)
  #summary(continuous) -> factors age,BMI,paclitaxel,Docetaxel and Taxan_type with other scales
  
  
  
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
  
  if (covariate==TRUE){
    test_data_rescale <- testing_phenotype[5:11] %>%
      mutate_if(is.numeric, funs(as.numeric(scale(.))))
    #head(train_data_rescale)
    
    #combine x and y for log reg models
    Testing <- as.data.frame(cbind(Y_testing,X_testing,
                                   test_data_rescale$Age,test_data_rescale$Diabetes,
                                   test_data_rescale$BMI_survey,test_data_rescale$Taxan_type))
    colnames(Testing) <- c("Y",colnames(X_testing),"Age", "Diabetes","BMI", "Taxane_type")
  } else {
    #combine x and y for log reg models
    Testing <- as.data.frame(cbind(Y_testing,X_testing))
    colnames(Testing) <- c("Y", colnames(X_testing))
  }
  
  #print(Training$Y)
  #Testing <- as.data.frame(cbind(Y_testing, X_testing))
  #colnames(Testing) <- c("Y", colnames(X_testing))
  
  #fit log reg model
  if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
    training_exclude <- which(Y_training == -9)
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
    pred.prob <- predict(fit,type='response')
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    #cutoff<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    cutoff_high<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    cutoff_low<- low_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    
    roc_auc_plot <- auc_plotting(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text,cutoff_low)
    print(table(preds[-training_exclude] >= cutoff_low, Y_training[-training_exclude] == 1))
    print(table(preds[-training_exclude] >= cutoff_high, Y_training[-training_exclude] == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff_low)
  } else { #om inga patienter har missing fenotyp
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training, family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    #cutoff<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    cutoff_high<- best_cutoff(preds, Y_training)
    cutoff_low<- low_cutoff(preds, Y_training)
    
    roc_auc_plot <- auc_plotting(preds, Y_training, phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds, Y_training, phenotype_text, variants_text, cutoff_low)
    print(table(preds >= cutoff_low, Y_training == 1))
    #printcutoff_low
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression",cutoff_low)
  }
  
  
  #calculate RMSE
  print(sqrt(mean((Training$Y - preds)^2)))
  
  brierScore <- mean((preds-Training$Y)^2)
  print(brierScore)
  # brierScore <- mean(fit$residuals^2)
  # print(brierScore)
  
  ## make the score and class
  class <- fit$y
  # score = log odds
  score <- qlogis(fit$fitted.values)
  
  ## rocit object
  rocit_emp <- rocit(score = score, 
                     class = class, 
                     method = "emp")
  rocit_bin <- rocit(score = score, 
                     class = class, 
                     method = "bin")
  rocit_non <- rocit(score = score, 
                     class = class, 
                     method = "non")
  print(summary(rocit_emp))
  print(summary(rocit_bin))
  print(summary(rocit_non))
  
  print(ciAUC(rocit_emp))
  print(ciAUC(rocit_bin,delong = TRUE))
  print(ciAUC(rocit_bin,delong = TRUE,logit = TRUE))
  # set.seed(200)
  # ciAUC_boot<-ciAUC(rocit_non,level = 0.9, nboot = 200)
  # print(ciAUC_boot)
  
  #testdata
  if(sum(Y_testing == -9) > 0) {
    testing_exclude <- which(Y_testing == -9)
    preds <- c()
    preds <- predict.glm(fit, newdata = Testing[, -1], type = "response")
    #print(summary(fit))
    roc_auc_plot_test <- auc_plotting(preds[-testing_exclude], Y_testing[-testing_exclude], phenotype_text, variants_text)
    accuracy_plot_test <- accuracy_plot(preds[-testing_exclude], Y_testing[-testing_exclude], phenotype_text, variants_text, cutoff_low)
    print(table(preds[-testing_exclude] >= cutoff_low, Y_testing[-testing_exclude] == 1))
    print(table(preds[-testing_exclude] >= cutoff_high, Y_testing[-testing_exclude] == 1))
    #print(table(preds[-testing_exclude] >= best_cutoff(preds[-testing_exclude], Y_testing[-testing_exclude]), Y_testing[-testing_exclude] == 1))
    #print(table(preds[-testing_exclude] >= .5, Y_testing[-testing_exclude] == 1))
    log_reg_pred_plot_test <- prediction_plot(preds, Y_training, phenotype_text, variants_text, "Logistic Regression", cutoff_low)
  } else {
    preds <- c()
    preds <- predict.glm(fit, newdata = Testing[, -1], type = "response")
    #print(summary(fit))
    roc_auc_plot_test <- auc_plotting(preds, Y_testing, phenotype_text, variants_text)
    accuracy_plot_test <- accuracy_plot(preds, Y_testing, phenotype_text, variants_text, cutoff_low)
    print(table(preds >= cutoff_low, Y_testing == 1))
    print(table(preds>= cutoff_high, Y_testing == 1))
    #print(table(preds >= best_cutoff(preds, Y_testing), Y_testing == 1))
    #print(table(preds[-training_exclude] >= .5, Y_testing == 1))
    log_reg_pred_plot_test <- prediction_plot(preds, Y_testing, phenotype_text, variants_text, "Logistic Regression",cutoff_low)
  }
  
  #calculate RMSE
  print(sqrt(mean((Y_testing - preds)^2)))
  mse <- MSE(preds,Y_testing)
  print(mse)
  
  brierScore <- mean((preds-Y_testing)^2)
  print(brierScore)
  # brierScore <- mean(fit$residuals^2)
  # print(brierScore)
  ## make the score and class
  class <- Y_testing
  # score = log odds
  score <- qlogis(preds)
  
  ## rocit object
  rocit_emp <- rocit(score = score, 
                     class = class, 
                     method = "emp")
  rocit_bin <- rocit(score = score, 
                     class = class, 
                     method = "bin")
  rocit_non <- rocit(score = score, 
                     class = class, 
                     method = "non")
  print(summary(rocit_emp))
  print(summary(rocit_bin))
  print(summary(rocit_non))
  
  print(ciAUC(rocit_emp))
  print(ciAUC(rocit_bin,delong = TRUE))
  print(ciAUC(rocit_bin,delong = TRUE,logit = TRUE))

    
  if (covariate==TRUE) {
    p <- ggarrange(log_reg_pred_plot, accuracy_plot, roc_auc_plot,log_reg_pred_plot_test, accuracy_plot_test, roc_auc_plot_test, ncol = 3, nrow = 2,labels = c("A", "B", "C", "A", "B", "C"))
    annotate_figure(p, top = text_grob(paste("Train & Test: For", phenotype_text, "using", variants_text, "with covariates")),left = text_grob("Train\n\n\nTest", face="bold"))
  } else {
    p <- ggarrange(log_reg_pred_plot, accuracy_plot, roc_auc_plot, log_reg_pred_plot_test, accuracy_plot_test, roc_auc_plot_test, ncol = 3, nrow = 2,labels = c("A", "B", "C", "A", "B", "C"))
    annotate_figure(p, top = text_grob(paste("Train & Test: For", phenotype_text, "using", variants_text)), left = text_grob("Train\n\nTest", face="bold"))
  }
  
}
```

Example with Numbness in feet for C2.

```{r}
test_log_reg_model(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                     "Numbness_A1_34_C0_12", numb_400cpdb_high_imp$SNP, "35SNVs from FDR perm 398SNVs, CPDB p<0.0025, review p<0.05, VarImp>1.67", covariate=TRUE)
                  length(numb_400cpdb_high_imp$SNP)
                  length(unique(numb_400cpdb_high_imp$GeneName_GRCh38_103_gtf))
```

## Validation visualization

Function for validating the model in the test set. Calculating rmse,
brierscore and rocit for the model as well as a comparision of the
prediction distribution in a clean plot.

```{r}
roc_auc_joint_plot <- function(roc_auc_plot_test,roc_auc_plot_train, preds_train, Y_training, preds_test, Y_testing) {
  
  train_data<-roc_auc_plot_train$data
  train_data$type<- "train"
  test_data<-roc_auc_plot_test$data
  test_data$type<- "test"
  
  evaluate_preds_train <- prediction(preds_train, Y_training)
  # calculate probabilities for TPR/FPR for predictions
  auc_train <- performance(evaluate_preds_train,"auc")@y.values[[1]] # shows calculated AUC for model
  
  evaluate_preds_test <- prediction(preds_test, Y_testing)
  # calculate probabilities for TPR/FPR for predictions
  auc_test <- performance(evaluate_preds_test,"auc")@y.values[[1]] # shows calculated AUC for model
  
  
  my_data<- rbind(train_data,test_data)
  #print(my_data)
  p<- ggplot(data=my_data, aes(x=V1, y=V2, group=type)) +
      geom_line() + 
      geom_segment(aes(x=0, y=0, xend=1, yend=1), linetype = "dashed", color="black") +  
      labs(x = "False positive rate", y = "True positive rate")  +
      coord_cartesian(ylim = c(0, 1))  + 
      ggtitle("") + #ROC-curve") +
      annotate("text", x=0.75, y=0.25, label = paste("Train AUC = ", round(auc_train, 4) * 100, "%", sep = ""), size = 8/.pt) + 
      annotate("text", x=0.75, y=0.20, label = paste("Test AUC = ", round(auc_test, 4) * 100, "%", sep = ""), size = 8/.pt) + 
      theme_cowplot(12) + theme(axis.text=element_text(size=10),
                                axis.title=element_text(size=10),
                                plot.title=element_text(size=10,face="bold"),
                                legend.title=element_text(size=10,face="bold"),
                                legend.text=element_text(size=10), 
                                axis.text.x = element_text(colour = "black"),
                                axis.text.y = element_text(colour = "black"))
    
}

test_log_reg_model_A <- function(training_genotypes,training_phenotype, 
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
  X_testing <- data[[3]] #snv
  Y_testing <- data[[4]] #phenotype
  
  ##checking continous variables in set
  #continuous <-select_if(training_phenotype, is.numeric)
  #summary(continuous) -> factors age,BMI,paclitaxel,Docetaxel and Taxan_type with other scales
  
  
  
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
  
  if (covariate==TRUE){
    test_data_rescale <- testing_phenotype[5:11] %>%
      mutate_if(is.numeric, funs(as.numeric(scale(.))))
    #head(train_data_rescale)
    
    #combine x and y for log reg models
    Testing <- as.data.frame(cbind(Y_testing,X_testing,
                                   test_data_rescale$Age,test_data_rescale$Diabetes,
                                   test_data_rescale$BMI_survey,test_data_rescale$Taxan_type))
    colnames(Testing) <- c("Y",colnames(X_testing),"Age", "Diabetes","BMI", "Taxane_type")
  } else {
    #combine x and y for log reg models
    Testing <- as.data.frame(cbind(Y_testing,X_testing))
    colnames(Testing) <- c("Y", colnames(X_testing))
  }
  
  #print(Training$Y)
  #Testing <- as.data.frame(cbind(Y_testing, X_testing))
  #colnames(Testing) <- c("Y", colnames(X_testing))
  
  #fit log reg model
  if(sum(Y_training == -9) > 0) { #om det finns n?gon patient som har missing fenotype
    training_exclude <- which(Y_training == -9)
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training[-training_exclude,], family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    pred.prob <- predict(fit,type='response')
    brierScore <- mean((pred.prob-Training$Y)^2)
    print(brierScore)
    brierScore <- mean(fit$residuals^2)
    print(brierScore)
    #print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff_high<- best_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    cutoff<- low_cutoff(preds[-training_exclude], Y_training[-training_exclude])
    
    roc_auc_plot <- auc_plotting(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds[-training_exclude], Y_training[-training_exclude], phenotype_text, variants_text, cutoff_high[[1]])
    #print(table(preds[-training_exclude] >= best_cutoff(preds[-training_exclude], Y_training[-training_exclude]), Y_training[-training_exclude] == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training[-training_exclude] == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, paste("TRAIN"), cutoff_high[[1]])
    preds_train<- preds
  } else { #om inga patienter har missing fenotyp
    fit <- c()
    fit <- glm(as.factor(Y)~., data = Training, family="binomial")
    preds <- c()
    preds <- predict.glm(fit, newdata = Training[, -1], type = "response")
    pred.prob <- predict(fit,type='response')
    brierScore <- mean((pred.prob-Training$Y)^2)
    print(brierScore)
    brierScore <- mean(fit$residuals^2)
    print(brierScore)
    #print(summary(fit))
    #print(confint(fit))
    #print(coefficients(fit))
    #print(preds)
    cutoff_high<- best_cutoff(preds, Y_training)
    cutoff<- low_cutoff(preds, Y_training)
    
    roc_auc_plot <- auc_plotting(preds, Y_training, phenotype_text, variants_text)
    accuracy_plot <- accuracy_plot(preds, Y_training, phenotype_text, variants_text, cutoff_high[[1]])
    #print(table(preds >= best_cutoff(preds, Y_training), Y_training == 1))
    #print(table(preds[-training_exclude] >= .5, Y_training == 1))
    log_reg_pred_plot <- prediction_plot(preds, Y_training, phenotype_text, variants_text, paste("TRAIN"), cutoff_high[[1]])
    preds_train<- preds
  }
  
  ## make the score and class
  class <- fit$y
  # score = log odds
  score <- qlogis(fit$fitted.values)
  
  ## rocit object
  rocit_emp <- rocit(score = score, 
                     class = class, 
                     method = "emp")
  rocit_bin <- rocit(score = score, 
                     class = class, 
                     method = "bin")
  rocit_non <- rocit(score = score, 
                     class = class, 
                     method = "non")
  print(summary(rocit_emp))
  print(summary(rocit_bin))
  print(summary(rocit_non))
  
  print(ciAUC(rocit_emp))
  print(ciAUC(rocit_bin,delong = TRUE))
  print(ciAUC(rocit_bin,delong = TRUE,logit = TRUE))
  set.seed(200)
  ciAUC_boot<-ciAUC(rocit_non,level = 0.9, nboot = 200)
  print(ciAUC_boot)
  
  #testdata
  if(sum(Y_testing == -9) > 0) {
    testing_exclude <- which(Y_testing == -9)
    preds <- c()
    preds <- predict.glm(fit, newdata = Testing[, -1], type = "response")
    #print(summary(fit))
    pred.prob <- predict(fit,type='response')
    brierScore <- mean((pred.prob-Testing$Y)^2)
    print(brierScore)
    brierScore <- mean(fit$residuals^2)
    print(brierScore)
    
    roc_auc_plot_test <- auc_plotting(preds[-testing_exclude], Y_testing[-testing_exclude], phenotype_text, variants_text)
    accuracy_plot_test <- accuracy_plot(preds[-testing_exclude], Y_testing[-testing_exclude], phenotype_text, variants_text, cutoff_high[[1]])
    #print(table(preds[-testing_exclude] >= cutoff, Y_testing[-testing_exclude] == 1))
    #print(table(preds[-testing_exclude] >= cutoff_high, Y_testing[-testing_exclude] == 1))
    #print(preds)
    log_reg_pred_plot_test <- prediction_plot(preds, Y_training, phenotype_text, variants_text,paste("TEST"), cutoff_high[[1]])
    log_reg_pred_plot_test <- log_reg_pred_plot_test + scale_x_discrete(limits=tox_test_groups) +
                                                      scale_color_manual(values=tox_test_dark_cols) +
                                                      scale_fill_manual(values=tox_test_cols, name = "Toxicity: ")
    preds_test<- preds
  } else {
    preds <- c()
    preds <- predict.glm(fit, newdata = Testing[, -1], type = "response")
    pred.prob <- predict(fit,type='response')
    brierScore <- mean((pred.prob-Testing$Y)^2)
    print(brierScore)
    brierScore <- mean(fit$residuals^2)
    print(brierScore)
    #print(summary(fit))
    roc_auc_plot_test <- auc_plotting(preds, Y_testing, phenotype_text, variants_text)
    accuracy_plot_test <- accuracy_plot(preds, Y_testing, phenotype_text, variants_text, cutoff_high[[1]])
    #print(table(preds >= cutoff, Y_testing == 1))
    #print(table(preds>= cutoff_high, Y_testing == 1))
    #print(table(preds[-training_exclude] >= .5, Y_testing == 1))
    log_reg_pred_plot_test <- prediction_plot(preds, Y_testing, phenotype_text, variants_text, paste("TEST"), cutoff_high[[1]])
    log_reg_pred_plot_test <- log_reg_pred_plot_test + scale_x_discrete(limits=tox_test_groups) +
                                                      scale_color_manual(values=tox_test_dark_cols) +
                                                      scale_fill_manual(values=tox_test_cols, name = "Toxicity: ")
    preds_test<- preds
    #print(preds)
  }
  if (sum(Y_training == -9) > 0){
    roc_auc_plot<- roc_auc_joint_plot(roc_auc_plot_test,roc_auc_plot, preds_train[-training_exclude], Y_training[-training_exclude], preds_test, Y_testing)
  }else if (sum(Y_testing == -9) > 0) {
    roc_auc_plot<- roc_auc_joint_plot(roc_auc_plot_test,roc_auc_plot, preds_train, Y_training, preds_test[-testing_exclude], Y_testing[-testing_exclude])
  }else if (sum(Y_testing == -9) > 0 | sum(Y_training == -9) > 0){
    roc_auc_plot<- roc_auc_joint_plot(roc_auc_plot_test,roc_auc_plot, preds_train[-training_exclude], Y_training[-training_exclude], preds_test[-testing_exclude], Y_testing[-testing_exclude])
  }else {
    roc_auc_plot<- roc_auc_joint_plot(roc_auc_plot_test,roc_auc_plot, preds_train, Y_training, preds_test, Y_testing)
  }
  
  ##for joint plot
  #pred_joint_plot<-prediction_joint_plot(preds_train, Y_training, preds_test, Y_testing, phenotype_text, variants_text, covariate)
  #ggarrange(pred_joint_plot, roc_auc_plot, ncol = 2, nrow = 1)
  
  if (covariate==TRUE) {
    p<-ggarrange(log_reg_pred_plot, log_reg_pred_plot_test, ncol = 2, nrow = 1)
    annotate_figure(p, top = text_grob(substitute(paste(bold('Numbness in feet using '), bold(variants_text), bold(' with covariates'))))) #paste(bold(phenotype_text), bold(' using '), bold(variants_text), bold(' with covariates')))))
  } else {
    p<-ggarrange(log_reg_pred_plot, log_reg_pred_plot_test, ncol = 2, nrow = 1)
    annotate_figure(p, top = text_grob(substitute(paste(bold('Numbness in feet using '), bold(variants_text)))))
  }
}

```

Example with Numbness in feet for C2.

```{r}
test_log_reg_model_A(training_genotypes,training_phenotype, testing_genotypes, testing_phenotype,
                                     "Numbness_A1_34_C0_12", numb_400cpdb_high_imp$SNP, "35SNVs from FDR perm 398SNVs, CPDB p<0.0025, review p<0.05, VarImp>1.67", covariate=TRUE)
                  length(numb_400cpdb_high_imp$SNP)
                  length(unique(numb_400cpdb_high_imp$GeneName_GRCh38_103_gtf))
```

Read the article to find more insight in the model building and
validation results.
