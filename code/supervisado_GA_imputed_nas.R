# Libraries
library(tidyverse)
library(ggplot2)
library(corrplot)
library(caret)
library(rpart)
library(rpart.plot)
library(MLmetrics)
library(randomForest)

# Load data after the feature selection from the genetic algorithm
data <- read.csv("/Users/alzorrcarri/projects_genyo/ml_farmacogen/data/data_feature_selection_GA_parellel_na_imputation_svm.tsv",
                 sep = "\t", header = TRUE)

# Separate information from the dataset
sample_id <- data$samples_id
class <- data$class
snps <- data[, -((ncol(data)-1):ncol(data))]

# Create a new dataframe with the relevant information
df <- snps

str(df)
# Correctly formatting variables
df[df == "True"] <- "1"
df[df == "False"] <- "0"
df <- df %>%
  mutate(across(everything(), as.integer))
str(df)
sum(as.numeric(lapply(df, is.integer))) == dim(df)[2] # TRUE

df$class <- as.factor(class)

# Verify everything's correct
sum(is.na(df)) # NAs = 0
summary(df$class)

#---------------------------------------------------------------------#
# Correlation study
cor_matrix <- cor(df[, -ncol(df)])
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.8, exact = TRUE)
filtered_df <- df[, -highly_correlated] # Delete highly-correlated mutations

#---------------------------------------------------------------------#
# Divide between train and test sets
set.seed(7)
trainIndex <- createDataPartition(filtered_df$class, p = 0.7, list = FALSE)
train <- filtered_df[trainIndex, ]
test  <- filtered_df[-trainIndex, ]

# Balancing train set with downsampling
summary(train$class)
train_under <- downSample(x = train[, -which(names(train) == "class")],  
                          y = train$class)
summary(train_under$Class)

# train_up <- upSample(x = train[, -which(names(train) == "class")],  
#                      y = train$class)
# summary(train_up$Class)
#---------------------------------------------------------------------#
# Decision tree
tree_control <- rpart.control(cp = 0.001,   # Reduce poda para permitir más ramas
                              minsplit = 1, # Mínimo de observaciones para dividir un nodo
                              minbucket = 1, # Mínimo de observaciones en nodos terminales
                              maxdepth = 10)  # Controla la profundidad del árbol
tree_model <- train(Class ~ .,
                    data = train_under,
                    method = "rpart",
                    trControl = trainControl(method = "cv", number = 3),
                    control = tree_control)
rpart.plot(tree_model$finalModel, box.palette = 'GnRd')

predictions <- predict(tree_model, newdata = test)
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$byClass 
conf_matrix$overall['Accuracy'] 

#-------------------------------------------------------------------#
# Random Forest
# Parameters
paramgrid <- expand.grid(
  .mtry = seq(from = 1, to = 10, length.out = 10)        
)

rf_model <- train(Class ~ .,  
                  data = train_under,  
                  method = "rf",  
                  trControl = trainControl(method = "cv", number = 3),
                  tuneGrid = paramgrid,
                  ntree = 150)
plot(rf_model)

predictions <- predict(rf_model, newdata = test)
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$byClass 
conf_matrix$overall['Accuracy'] 

varImpPlot(rf_model$finalModel, 
           main = "Importancia de Variables - Random Forest",
           n.var = 20)

#-------------------------------------------------------------------#
# SVM 
library(kernlab)

paramgrid <- expand.grid(
  sigma = seq(from = 0.005, to = 0.02, length.out = 20),  # Valores para sigma
  C = seq(from = 0.1, to = 1.0, length.out = 20)  # Valores para C
)
# Training with downsampled train set
modelo_svm <- train(
  Class ~ ., 
  data = train_under,
  method = "svmRadial",     # Kernel radial
  trControl = trainControl(method = 'cv', number = 3),
  tuneGrid = paramgrid
)
# Plotting accuracy as a function of the parameters
print(modelo_svm)
plot(modelo_svm)

predictions <- predict(modelo_svm, newdata = test)
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$overall['Accuracy'] # 0.5940594
conf_matrix$byClass # ¡Se obtienen aceptables/buenos resultados! (Comparado con lo que teníamos antes)
write.table(as.data.frame(conf_matrix$byClass),
            file = "metricas_svm_prediccion.csv", sep=",")

#--------------------------------------#
# # Train with upsampled train set 
#
# paramgrid <- expand.grid(
#   sigma = seq(from = 0.001, to = 0.02, length.out = 50),  # Valores para sigma
#   C = seq(from = 0.1, to = 1.2, length.out = 50)  # Valores para C
# )
# modelo_svm <- train(
#   Class ~ ., 
#   data = train_up,
#   method = "svmRadial",     # Kernel radial
#   trControl = trainControl(method = 'cv', number = 3),
#   tuneGrid = paramgrid
# )
#
# print(modelo_svm)
# plot(modelo_svm)
#
# predictions <- predict(modelo_svm, newdata = test)
# conf_matrix <- confusionMatrix(predictions, test$class)
# conf_matrix$overall['Accuracy']
# conf_matrix$byClass 

#-------------------------------------------------------------------#
# ElasticNet
library(glmnet)

paramgrid <- expand.grid(
  alpha = seq(0, 1, length.out = 100),
  lambda = seq(0.0001, 0.1, length.out = 10)
)

elasticnet_model <- train(Class ~ .,  
                          data = train_under,  
                          method = "glmnet",  
                          trControl = trainControl(method = "cv", number = 3),  # Validación cruzada de 10 pliegues
                          tuneGrid = paramgrid)  # Usar el paramgrid con alpha y lambda
# Plot model's results
plot(elasticnet_model) 

predictions <- predict(elasticnet_model, newdata = test)
conf_matrix <- confusionMatrix(predictions, test$class)
conf_matrix$overall['Accuracy']
conf_matrix$byClass 


#-------------------------------------------------------------------#
#-------------------------------------------------------------------#
# XGBoost
library(xgboost)
# Required format
dtrain <- train_under[, -ncol(train_under)] %>%
  mutate(across(everything(), as.numeric))
dtest <- test[, -ncol(test)] %>%
  mutate(across(everything(), as.numeric))
dtrain <- as.matrix(dtrain)
dtest <- as.matrix(dtest)

# Define model's parameters
params <- list(
  objective = "multi:softmax",  # Para clasificación multiclase
  eval_metric = "mlogloss",      # Log-loss para multiclase
  num_class = length(unique(train_under$Class)),  # Número de clases
  eta = 0.3,  # Tasa de aprendizaje
  max_depth = 10,  # Profundidad del árbol
  subsample = 1.0,  # Submuestreo para evitar overfitting
  colsample_bytree = 0.6  # Fracción de características usadas por árbol
)
# Train the model
xgb_model <- xgboost(
  params = params,
  data = dtrain,
  label = as.integer(train_under$Class) - 1,
  nrounds = 200,  
  early_stopping_rounds = 20,  
  verbose = 1  
)
# Model's prediction with test set
predictions <- predict(xgb_model, dtest)
conf_matrix <- confusionMatrix(as.factor(predictions), as.factor(as.integer(test$class)-1))
conf_matrix$byClass
write.table(as.data.frame(conf_matrix$byClass),
          file = "metricas_xgboost_prediccion.csv", sep=",")
# Variables importance
importance <- xgb.importance(feature_names = colnames(train_under[, -ncol(train_under)]), model = xgb_model)
xgb.plot.importance(importance, top_n = 50, col = "aquamarine3")
write_tsv(importance, file="xgboost_feature_importance.tsv")
save(xgb_model, file = "xgboost_model.RData")

# SHAP values
library(shapviz)
# Get SHAP values for each observation
shap_values <- shapviz(xgb_model, X_pred = dtrain)
names(shap_values) <- levels(train_under$Class)
# Feature importance
sv_importance(shap_values, 
              kind = 'bar',
              max_display = 15,
              show_numbers = TRUE)
sv_importance(shap_values, 
              kind = 'bee',
              max_display = 15,
              show_numbers = TRUE)
save(shap_values, file = "shap_xgboost.RData")

#-------------------------------------------------------------------#
# Train a XGBoost model of good vs. bad responders
# and calculate SHAP values
bg_df <- df %>%
  filter(class == "Good" | class == "Bad")
bg_df$class <- factor(bg_df$class, levels = c("Good", "Bad"))

set.seed(7)
bg_trainIndex <- createDataPartition(bg_df$class, p = 0.7, list = FALSE)
bg_train <- bg_df[bg_trainIndex, ]
bg_test  <- bg_df[-bg_trainIndex, ]
# Downsampling
summary(bg_train$class)
bg_train_under <- downSample(x = bg_train[, -which(names(bg_train) == "class")],  
                          y = bg_train$class)
summary(bg_train_under$Class)

bg_dtrain <- bg_train_under[, -ncol(bg_train_under)] %>%
  mutate(across(everything(), as.numeric))
bg_dtest <- bg_test[, -ncol(bg_test)] %>%
  mutate(across(everything(), as.numeric))
bg_dtrain <- as.matrix(bg_dtrain)
bg_dtest <- as.matrix(bg_dtest)

library(xgboost)
params <- list(
  objective = "multi:softmax",
  eval_metric = "mlogloss",      
  num_class = length(unique(bg_train_under$Class)), 
  eta = 0.3,  
  max_depth = 10, 
  subsample = 1.0, 
  colsample_bytree = 0.4  
)
bg_xgb_model <- xgboost(
  params = params,
  data = bg_dtrain,
  label = as.integer(bg_train_under$Class) - 1,
  nrounds = 200,
  early_stopping_rounds = 20,  
  verbose = 1  
)
# Predictions over test
bg_predictions <- predict(bg_xgb_model, bg_dtest)
bg_conf_matrix <- confusionMatrix(as.factor(bg_predictions), as.factor(as.integer(bg_test$class)-1))
bg_conf_matrix$byClass
write_tsv(data.frame(metric=names(bg_conf_matrix$byClass), values=bg_conf_matrix$byClass),
          file = "bg_metricas_xgboost_predicción. tsv")
# Feature importance
bg_importance <- xgb.importance(feature_names = colnames(bg_train_under[, -ncol(bg_train_under)]), 
                             model = bg_xgb_model)
xgb.plot.importance(bg_importance, top_n = 20, col = "aquamarine3")
write_tsv(bg_importance, file="bg_xgboost_feature_importance.tsv")

# SHAP
library(shapviz)
bg_shap_values <- shapviz(bg_xgb_model, X_pred = bg_dtrain)
names(bg_shap_values) <- levels(bg_train_under$Class)
sv_importance(bg_shap_values, kind = 'bee', max_display = 20, show_numbers = TRUE)
sv_importance(bg_shap_values, kind = 'bar', max_display = 20, show_numbers = TRUE)

