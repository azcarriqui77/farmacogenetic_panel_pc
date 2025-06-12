# Libraries
library(readxl)
library(lubridate)
library(tidyverse)
library(ggplot2)
library(caret)
library(rpart)
library(rpart.plot)
library(MLmetrics)
library(randomForest)
library(xgboost)
library(Ckmeans.1d.dp)

# Data loading
clinical_data <- read_xlsx('/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/Variables clinicas Analisis 1 BRvsMR 080425.xlsx',
                           na = "NA")
str(clinical_data)

# Columns to use
colnames(clinical_data)
selected_clinical_data <- clinical_data[-c(167), c(1, 3, 5, 6, 7, 9)]

# Varible type formatting
selected_clinical_data$`Fecha de diagnóstico` <- as_date(as.numeric(selected_clinical_data$`Fecha de diagnóstico`),
                                                         origin = "1899-12-30")
selected_clinical_data$`Fecha de diagnóstico` <- as.Date(selected_clinical_data$`Fecha de diagnóstico`)
selected_clinical_data$`Fecha de Nacimiento` <- as.Date(selected_clinical_data$`Fecha de Nacimiento`)

# Patient's age at the diagnosis moment
selected_clinical_data$`Edad` <- selected_clinical_data$`Fecha de diagnóstico` - selected_clinical_data$`Fecha de Nacimiento`
selected_clinical_data$Edad <- as.integer(selected_clinical_data$Edad / 365)

# Delete variables 'Fecha de nacimiento' and 'Fecha de diagnóstico'
selected_clinical_data <- selected_clinical_data[, -c(3, 4)]

# Correctly formatting variables
selected_clinical_data$`Tipo de respuesta` <- as.factor(selected_clinical_data$`Tipo de respuesta`)
summary(selected_clinical_data$`Tipo de respuesta`)

unique(selected_clinical_data$`ISUP diagnóstico`)
selected_clinical_data$`ISUP diagnóstico`[selected_clinical_data$`ISUP diagnóstico` == "2 / PR 5"] <- "2"
selected_clinical_data$`ISUP diagnóstico`[selected_clinical_data$`ISUP diagnóstico` == "1 / PR 3"] <- "1"
selected_clinical_data$`ISUP diagnóstico`[selected_clinical_data$`ISUP diagnóstico` == "3 / PR 5"] <- "3"
selected_clinical_data$`ISUP diagnóstico` <- as.factor(selected_clinical_data$`ISUP diagnóstico`)
summary(selected_clinical_data$`ISUP diagnóstico`)

typeof(selected_clinical_data$`PSA diagnóstico`)
selected_clinical_data$`PSA diagnóstico`[selected_clinical_data$`PSA diagnóstico` == ">20"] <- 20
selected_clinical_data$`PSA diagnóstico` <- as.numeric(selected_clinical_data$`PSA diagnóstico`)
selected_clinical_data$`PSA diagnóstico`[selected_clinical_data$`PSA diagnóstico` >= 20] <- 20

# Delete rows with NAs
selected_clinical_data <- na.omit(selected_clinical_data)
summary(selected_clinical_data$`Tipo de respuesta`)

selected_clinical_data$ID <- gsub("_", "-", selected_clinical_data$ID)

rm(clinical_data) # Already not necessary 

# Formatting clinical data colnames
colnames(selected_clinical_data) <- c("ID", "Response", "PSA", "ISUP", "Age")


# Load SNPs data
data <- read.csv("/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/data_feature_selection_GA_parellel_na_imputation_svm.tsv",
                 sep = "\t", header = TRUE)

sample_id <- data$samples_id
class <- data$class
snps <- data[, -((ncol(data)-1):ncol(data))]

# Correctly formatting variables
snps[snps == "True"] <- "1"
snps[snps == "False"] <- "0"
snps <- snps %>%
  mutate(across(everything(), as.integer))

snps$class <- as.factor(class)

# SNPs sets to be considered
set1_snps <- c("rs1076563", "rs197922", "rs7248783", "rs9915058", "rs7520", "rs1456908")
set2_snps <- c("rs4245738", "rs4252717", "rs4252718", "rs2290855", "rs2290854", "rs1870029", "rs12718244", "rs1993445",
               "rs1456908", "rs2277311", "rs555456", "rs1044569", "rs1044564", "rs9915058", "rs7520", "rs66487517", 
               "rs608690", "rs687127", "rs11170175", "rs589985")


# CLASSIFICATION ALGORITHMS. LOOKING FOR THE BEST COMBINATION
# Only clinical data
## Train/test division of dataset
set.seed(7)
trainIndex <- createDataPartition(selected_clinical_data$Response, p = 0.8, list = FALSE)
train <- selected_clinical_data[trainIndex, -c(1)]
test  <- selected_clinical_data[-trainIndex, -c(1)]

## Downsampling for balacing
summary(train$Response)
train_under <- downSample(x = train[, -which(names(train) == "Response")],  
                          y = train$Response
                          )
summary(train_under$Class)

## Random Forest
### Parameters
paramgrid <- expand.grid(
  .mtry = seq(from = 1, to = 4, length.out = 10)        # Valores para mtry
)
### Model training
rf_model <- train(Class ~ .,  
                  data = train_under,  
                  method = "rf",  
                  trControl = trainControl(method = "cv", number = 3),
                  tuneGrid = paramgrid,
                  ntree = 200)
plot(rf_model)
#### Checking with test set
predictions <- predict(rf_model, newdata = test[, -which(names(test) == "Response")])
conf_matrix <- confusionMatrix(predictions, test$Response)
conf_matrix$byClass
conf_matrix$overall['Accuracy'] 

varImpPlot(rf_model$finalModel, 
           main = "Importancia de Variables - Random Forest",
           n.var = 6)


# Set 1 of 6 SNPs + Clinical Data
## Filter SNP set
snps_colnames <- colnames(snps)
snps_colnames <- sub("_.*", "", snps_colnames)
snps_cols <- snps_colnames %in% set1_snps
sum(snps_cols)

df <- snps[, snps_cols]
df$ID <- sample_id
df$class <- class

## Delete control cases
df <- df %>%
  filter(class != "Control")

## Merge with clinical data
df <- merge(df, selected_clinical_data, by = "ID")
df$class <- as.factor(df$class)
summary(df$class)


set.seed(7)
trainIndex <- createDataPartition(df$class, p = 0.8, list = FALSE)
train <- df[trainIndex, -c(1, 12)] # Eliminamos columnas ID y Tipo de respuesta
test  <- df[-trainIndex, -c(1, 12)]

## Downsampling
summary(train$class)
train_under <- downSample(x = train[, -which(names(train) == "class")],  
                          y = train$class)
summary(train_under$Class)

## Random Forest
paramgrid <- expand.grid(
  .mtry = seq(from = 1, to = 12, length.out = 10)       
)
rf_model <- train(Class ~ .,  
                  data = train_under,  
                  method = "rf",  
                  trControl = trainControl(method = "cv", number = 3),
                  tuneGrid = paramgrid,
                  ntree = 150)
plot(rf_model)

predictions <- predict(rf_model, newdata = test[, -which(names(test) == "class")])
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$byClass
conf_matrix$overall['Accuracy'] 

varImpPlot(rf_model$finalModel, 
           main = "Importancia de Variables - Random Forest",
           n.var = 10)



# Set 2 of 20 SNPs + Clinical Data
snps_colnames <- colnames(snps)
snps_colnames <- sub("_.*", "", snps_colnames)
snps_cols <- snps_colnames %in% set2_snps
sum(snps_cols)

df <- snps[, snps_cols]
df$ID <- sample_id
df$class <- class

df <- df %>%
  filter(class != "Control")

df <- merge(df, selected_clinical_data, by = "ID")
df$class <- as.factor(df$class)
summary(df$class)

set.seed(7)
trainIndex <- createDataPartition(df$class, p = 0.8, list = FALSE)
train <- df[trainIndex, -c(1, 20)] 
test  <- df[-trainIndex, -c(1, 20)]


summary(train$class)
train_under <- downSample(x = train[, -which(names(train) == "class")],  
                          y = train$class)
summary(train_under$Class)

## Random Forest
paramgrid <- expand.grid(
  .mtry = seq(from = 1, to = 20, length.out = 20)   
)

rf_model <- train(Class ~ .,  
                  data = train_under,  
                  method = "rf",  
                  trControl = trainControl(method = "cv", number = 3),
                  tuneGrid = paramgrid,
                  ntree = 150)
plot(rf_model)

predictions <- predict(rf_model, newdata = test[, -which(names(test) == "class")])
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$byClass
conf_matrix$overall['Accuracy'] 

varImpPlot(rf_model$finalModel, 
           main = "Importancia de Variables - Random Forest",
           n.var = 10)


# LOOKING FOR THE BEST ALGORITHM FOR THE CLASSIFICATION MODEL
# CHOOSE SET1 OF 6 SNPS + CLINICAL DATA

# RF (done)

# XGBOOST
dtrain <- train_under[, -ncol(train_under)] %>%
  mutate(across(everything(), as.numeric))
dtest <- test[, -which(colnames(test) == "class")] %>%
  mutate(across(everything(), as.numeric))
dtrain <- as.matrix(dtrain)
dtest <- as.matrix(dtest)

params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  eta = 0.1, 
  max_depth = 10,  
  subsample = 1.0,  
  colsample_bytree = 1.0
)

xgb_model <- xgboost(
  params = params,
  data = dtrain,
  label = as.integer(train_under$Class) - 1,
  nrounds = 200,  
  early_stopping_rounds = 5, 
  verbose = 1  
)

predictions <- predict(xgb_model, dtest)
conf_matrix <- confusionMatrix(as.factor(round(predictions)), as.factor(as.integer(test$class)-1))
conf_matrix$byClass
## Importance of variables
importance <- xgb.importance(feature_names = colnames(train_under[, -ncol(train_under)]), model = xgb_model)
xgb.ggplot.importance(importance, top_n = 7, n_clusters = 1,
                    xlab = "Importance", ylab = "Features", rel_to_first = TRUE) +
  ggtitle("") +
  theme_minimal() +
  scale_fill_manual(values = c("aquamarine3")) +
  theme(legend.position = "none")

## SHAP values
library(shapviz)

shap_values <- shapviz(xgb_model, X_pred = dtrain)
## Feature importance
sv_importance(shap_values, 
              kind = 'bar',
              max_display = 15,
              show_numbers = TRUE)
sv_importance(shap_values, 
              kind = 'bee',
              max_display = 15,
              show_numbers = TRUE)

# ElasticNet
library(glmnet)

## Normalise train and test set using Z-score
norm_train <- scale(dtrain)
norm_test <- scale(dtest)

enet_model <- cv.glmnet(x = dtrain,
                     y = as.factor(train_under[, which(colnames(train_under) == "Class")]),
                     alpha = 0.5,
                     family = "binomial")
plot(enet_model)
enet_model$lambda.min

predictions <- predict(enet_model, 
                       newx = dtest, 
                       s = "lambda.min",
                       type = "response")
predictions <- ifelse(predictions > 0.5, 1, 0)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(as.integer(test$class)-1))
conf_matrix$byClass

coefs <- as.matrix(coef(enet_model, s = "lambda.min"))
important_features <- coefs[order(abs(coefs[,1]), decreasing = TRUE),]
as.matrix(important_features)

# Repetw the best algorithm without ISUP information
train_under <- train_under %>%
  select(-"ISUP diagnóstico")
test <- test %>%
  select(-"ISUP diagnóstico")

## Random Forest
paramgrid <- expand.grid(
  .mtry = seq(from = 1, to = 11, length.out = 10)
)

rf_model <- train(Class ~ .,  
                  data = train_under,  
                  method = "rf",  
                  trControl = trainControl(method = "cv", number = 3),
                  tuneGrid = paramgrid,
                  ntree = 150)
plot(rf_model)

predictions <- predict(rf_model, newdata = test[, -which(names(test) == "class")])
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$byClass
conf_matrix$overall['Accuracy'] 

varImpPlot(rf_model$finalModel, 
           main = "Importancia de Variables - Random Forest",
           n.var = 10)






