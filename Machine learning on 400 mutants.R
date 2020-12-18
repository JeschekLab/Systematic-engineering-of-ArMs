library(tidyverse)
library(Interpol)
library(factoextra)
library(data.table)
library(xgboost)
library(mltools)
library(rsample)
library(broom)
library(e1071)
library(keras)

#Preparing the data sets

#Load results of Sav 112X K121X screening ----
#Set working directory do location of the "Combined data 400 mutants" file
setwd("") 

combined_data <- loadWorkbook("Combined data 400 mutants.xlsx")

#>>> Metathesis
metathesis <- readWorksheet(combined_data, sheet = "Screening_metathesis") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt),
            norm_wt_sd = sd(norm_wt), .groups = "drop")  %>%
  mutate(aa112 = factor(aa112), aa121 = factor(aa121)) %>%
  select(-c("norm_wt_sd"))

#>>> Deallylation (coumarin)
coumarin <- readWorksheet(combined_data, sheet = "Screening_deallylation_coumarin") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop") %>%
  mutate(aa112 = factor(aa112), aa121 = factor(aa121)) %>%
  select(-c("norm_wt_sd"))

#>>> Deallylation (indole)
RuIndole <- readWorksheet(combined_data, sheet = "Screening_deallylation_indole") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop") %>%
  mutate(aa112 = factor(aa112), aa121 = factor(aa121)) %>%
  select(-c("norm_wt_sd"))

#>>> Hydroamination
AuIndole <- readWorksheet(combined_data, sheet = "Screening_Hydroamination") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop") %>%
  mutate(aa112 = factor(aa112), aa121 = factor(aa121)) %>%
  select(-c("norm_wt_sd"))

#>>> Hydroarylation
AuFluo <- readWorksheet(combined_data, sheet = "Screening_hydroarylation") %>%
  separate(Variant, into = c("aa112", "aa121"), sep = 1) %>%
  group_by(aa112, aa121) %>%
  summarise(norm_wt_mean = mean(norm_wt), 
            norm_wt_sd = sd(norm_wt), .groups = "drop") %>%
  mutate(aa112 = factor(aa112), aa121 = factor(aa121)) %>%
  select(-c("norm_wt_sd"))

#Split into training and test sets ----
set.seed(42)

AuIndole_split <- initial_split(AuIndole, prop = 0.799)
AuFluo_split <- initial_split(AuFluo, prop = 0.799)
RuIndole_split <- initial_split(RuIndole, prop = 0.799)
coumarin_split <- initial_split(coumarin, prop = 0.799)
metathesis_split <- initial_split(metathesis, prop = 0.799)

#Training data
AuIndole_train <- training(AuIndole_split)
AuFluo_train <- training(AuFluo_split)
RuIndole_train <- training(RuIndole_split)
coumarin_train <- training(coumarin_split)
metathesis_train <- training(metathesis_split)

#Testing data
AuIndole_test <- testing(AuIndole_split)
AuFluo_test <- testing(AuFluo_split)
RuIndole_test <- testing(RuIndole_split)
coumarin_test <- testing(coumarin_split)
metathesis_test <- testing(metathesis_split)

#Further splits for cross-validation
AuIndole_cv <- vfold_cv(AuIndole_train, v = 4)
AuFluo_cv <- vfold_cv(AuFluo_train, v = 4)
RuIndole_cv <- vfold_cv(RuIndole_train, v = 4)
coumarin_cv <- vfold_cv(coumarin_train, v = 4)
metathesis_cv <- vfold_cv(metathesis_train, v = 4)

AuIndole_folds <- AuIndole_cv %>%
  mutate(train = map(splits, ~training(.x)), validate = map(splits, ~testing(.x)))
AuFluo_folds <- AuFluo_cv %>%
  mutate(train = map(splits, ~training(.x)), validate = map(splits, ~testing(.x)))
RuIndole_folds <- RuIndole_cv %>%
  mutate(train = map(splits, ~training(.x)), validate = map(splits, ~testing(.x)))
coumarin_folds <- coumarin_cv %>%
  mutate(train = map(splits, ~training(.x)), validate = map(splits, ~testing(.x)))
metathesis_folds <- metathesis_cv %>%
  mutate(train = map(splits, ~training(.x)), validate = map(splits, ~testing(.x)))

#Training data with cross-validation folds for all reactions
training <- bind_rows("AuIndole" = AuIndole_folds,
                      "AuFluo" = AuFluo_folds,
                      "RuIndole" = RuIndole_folds,
                      "coumarin" = coumarin_folds,
                      "metathesis" = metathesis_folds,
                      .id = "reaction")


#Encoding of amino acids ----
#> Calculate PCscores ----
#(As described in Xu et al. 2020. DOI: 10.1021/acs.jcim.0c00073)

#Load AAindex data
data(AAindex)

AAindex_df <- as_tibble(AAindex) %>%
  as.data.frame(t(AAindex)) %>%
  rename("descriptor" = "name") %>%
  pivot_longer(-"descriptor", names_to = "AA", values_to = "value") %>%
  filter(descriptor != "Free_energy_in_beta-strand_region_(Munoz-Serrano,_1994)") %>% #Excluded because of multiple values per amino acid
  pivot_wider(names_from = "descriptor", values_from = "value")

#Principal component analysis
AA_PCA <- prcomp(AAindex_df[,2:532], center = TRUE, scale = TRUE)

indi <- get_pca_ind(AA_PCA)

PCA_descriptors <- bind_cols(AAindex_df$AA, as_tibble(indi$coord)) %>%
  select(1:12) %>%
  rename("aa" = "...1")

#> Load descriptors ----
#(Provide excel file "AA descriptors")

descriptor_book <- loadWorkbook("AA descriptors.xlsx")
descriptor_VHSE <- readWorksheet(descriptor_book, sheet = "VHSE")
descriptor_Zscales <- readWorksheet(descriptor_book, sheet = "Z-scales_1998")
descriptor_Barley <- readWorksheet(descriptor_book, sheet = "Barley")
descriptor_PCscores <- readWorksheet(descriptor_book, sheet = "PCscores")


#Function for applying the encoding ----
encode <- function(df, descriptor) {
  if (descriptor == 1) { 
    df_mod <- left_join(df, descriptor_VHSE, by = c("aa112" = "aminoacid")) %>%
      left_join(descriptor_VHSE, by = c("aa121" = "aminoacid"), suffix = c("_112", "_121")) %>%
      select(-c("aa112", "aa121"))
  }
  else if (descriptor == 2) { 
    df_mod <- left_join(df, descriptor_Barley, by = c("aa112" = "aminoacid")) %>%
      left_join(descriptor_Barley, by = c("aa121" = "aminoacid"), suffix = c("_112", "_121")) %>%
      select(-c("aa112", "aa121"))
  }
  else if (descriptor == 3) { 
    df_mod <- left_join(df, descriptor_PCscores, by = c("aa112" = "aminoacid")) %>%
      left_join(descriptor_PCscores, by = c("aa121" = "aminoacid"), suffix = c("_112", "_121")) %>%
      select(-c("aa112", "aa121"))
  }
  else if (descriptor == 4) { 
    df_mod <- left_join(df, descriptor_Zscales, by = c("aa112" = "aminoacid")) %>%
      left_join(descriptor_Zscales, by = c("aa121" = "aminoacid"), suffix = c("_112", "_121")) %>%
      select(-c("aa112", "aa121"))
  }
  else if (descriptor == 5) {
    df_mod <- data.table(df) %>%
      mutate(aa112 = factor(aa112, levels = c("F", "Y", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "S", "V", "A", "D", "E", "G")), 
             aa121 = factor(aa121, levels = c("F", "Y", "C", "W", "L", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "S", "V", "A", "D", "E", "G"))) %>%
      one_hot(dropUnusedLevels = FALSE, dropCols = TRUE) %>%
      as_tibble()
  }
  return(df_mod)
}

#Vectorize training data (including cross-validation folds)
train_encoded <- training %>%
  crossing(descriptor = c(1:5)) %>%
  mutate(train = map2(.x = train, .y = descriptor, ~encode(df = .x, descriptor = .y)),
         validate = map2(.x = validate, .y = descriptor, ~encode(df = .x, descriptor = .y)))


#Helper function for extracting the first element of a list
from_list <- function(list){
  x <- list[[1]]
  return(x)
}


#Gradient boosting ----

#Training data for gradient boosting
training_xg <- bind_rows("AuIndole" = AuIndole_train,
                         "AuFluo" = AuFluo_train,
                         "RuIndole" = RuIndole_train,
                         "coumarin" = coumarin_train,
                         "metathesis" = metathesis_train,
                         .id = "reaction")


#> Hyperparameter tuning via cross-validation ----
crossvalidation_xgb <- training_xg %>%
  group_by(reaction) %>%
  nest(aa = starts_with("aa"),
       data = c("norm_wt_mean")) %>%
  crossing(descriptor = c(1:5)) %>%
  mutate(aa = map2(.x = aa, .y = descriptor, ~encode(df = .x, descriptor = .y))) %>%
  crossing(eta = c(0.01, 0.05), 
           max.depth = c(10, 15),
           subsample = c(0.3, 0.5, 1),
           colsample_bytree = c(0.3, 0.667)) %>%
  rowwise() %>%
  mutate(xgbcv = list(xgb.cv(data = as.matrix(aa),
                             label = from_list(data),
                             nfold = 4,
                             nrounds = 500,
                             eta = eta,
                             max.depth = max.depth,
                             early_stopping_rounds = 10,
                             objective = "reg:squarederror",
                             subsample = subsample,
                             colsample_bytree = colsample_bytree,
                             prediction = TRUE,
                             stratified = TRUE,
                             verbose = FALSE))) %>%
  mutate(elog = list(xgbcv$evaluation_log),
         test_rmse_min = min(elog$test_rmse_mean),
         ntrees = xgbcv$best_ntreelimit,
         train_rmse = elog$train_rmse_mean[ntrees]) 

#Print best combinations of hyperparameters
crossvalidation_xgb %>% 
  group_by(reaction) %>% 
  slice_min(order_by = test_rmse_min, n = 1)


#> Training and testing of final models ----
#>>> Metathesis ----
metathesis_xgb <- metathesis_train %>%
  encode(descriptor = 3) %>%
  xgboost(data = as.matrix(select(., -c("norm_wt_mean"))),
          label = .$norm_wt_mean,
          nrounds = 125,
          eta = 0.05,
          max.depth = 15,
          objective = "reg:squarederror",
          subsample = 1,
          colsample_bytree = 0.667,
          prediction = TRUE,
          verbose = FALSE)

metathesis_xgb_results <- metathesis_test %>%
  encode(descriptor = 3) %>%
  mutate(prediction = predict(metathesis_xgb, as.matrix(select(., -c("norm_wt_mean")))),
         r2 = cor(norm_wt_mean, prediction)^2)


#>>> Deallylation (indole) ----
RuIndole_xgb <- RuIndole_train %>%
  encode(descriptor = 3) %>%
  xgboost(data = as.matrix(select(., -c("norm_wt_mean"))),
          label = .$norm_wt_mean,
          nrounds = 500,
          eta = 0.01,
          max.depth = 10,
          objective = "reg:squarederror",
          subsample = 0.5,
          colsample_bytree = 0.3,
          prediction = TRUE,
          verbose = FALSE)

RuIndole_xgb_results <- RuIndole_test %>%
  encode(descriptor = 3) %>%
  mutate(prediction = predict(RuIndole_xgb, as.matrix(select(., -c("norm_wt_mean")))),
         r2 = cor(norm_wt_mean, prediction)^2) 


#>>> Deallylation (coumarin) ----
coumarin_xgb <- coumarin_train %>%
  encode(descriptor = 3) %>%
  xgboost(data = as.matrix(select(., -c("norm_wt_mean"))),
          label = .$norm_wt_mean,
          nrounds = 136,
          eta = 0.05,
          max.depth = 15,
          objective = "reg:squarederror",
          subsample = 0.5,
          colsample_bytree = 0.667,
          prediction = TRUE,
          verbose = FALSE)

coumarin_xgb_results <- coumarin_test %>%
  encode(descriptor = 3) %>%
  mutate(prediction = predict(coumarin_xgb, as.matrix(select(., -c("norm_wt_mean")))),
         r2 = cor(norm_wt_mean, prediction)^2)


#>>> Hydroamination ----
AuIndole_xgb <- AuIndole_train %>%
  encode(descriptor = 1) %>%
  xgboost(data = as.matrix(select(., -c("norm_wt_mean"))),
          label = .$norm_wt_mean,
          nrounds = 499,
          eta = 0.01,
          max.depth = 10,
          objective = "reg:squarederror",
          subsample = 1,
          colsample_bytree = 0.667,
          prediction = TRUE,
          verbose = FALSE)

AuIndole_xgb_results <- AuIndole_test %>%
  encode(descriptor = 1) %>%
  mutate(prediction = predict(AuIndole_xgb, as.matrix(select(., -c("norm_wt_mean")))),
         r2 = cor(norm_wt_mean, prediction)^2) 

#>>> Hydroarylation ----
AuFluo_xgb <- AuFluo_train %>%
  encode(descriptor = 3) %>%
  xgboost(data = as.matrix(select(., -c("norm_wt_mean"))),
          label = .$norm_wt_mean,
          nrounds = 495,
          eta = 0.01,
          max.depth = 15,
          objective = "reg:squarederror",
          subsample = 0.5,
          colsample_bytree = 0.667,
          prediction = TRUE,
          verbose = FALSE)

AuFluo_xgb_results <- AuFluo_test %>%
  encode(descriptor = 3) %>%
  mutate(prediction = predict(AuFluo_xgb, as.matrix(select(., -c("norm_wt_mean")))),
         r2 = cor(norm_wt_mean, prediction)^2) 


#Support vector machine regression ----

#> Hyperparameter tuning by cross-validation
svm_crossvalidation <- training %>%
  crossing(descriptor = c(1:5)) %>%
  mutate(train = map2(.x = train, .y = descriptor, ~encode(df = .x, descriptor = .y)),
         validate = map2(.x = validate, .y = descriptor, ~encode(df = .x, descriptor = .y))) %>%
  crossing(cost = c(1, 2, 2^2, 2^3, 2^4, 2^6)) %>%
  rowwise() %>%
  mutate(svm = list(svm(data = train, norm_wt_mean ~ .,
                        kernel = "radial",
                        cost = cost)),
         predictions = list(predict(svm, select(validate, -c("norm_wt_mean")))),
         r2 = cor(validate$norm_wt_mean, predictions)^2,
         rmse = rmse(validate$norm_wt_mean, predictions)) %>%
  group_by(reaction, descriptor, cost) %>%
  summarise(r2_mean = mean(r2), 
            rmse_mean = mean(rmse),
            .groups = "drop")

#Print best hyperparamter sets
svm_crossvalidation %>%
  group_by(reaction) %>%
  slice_min(order_by = rmse_mean, n = 1)


#> Train and test final models ----
#>>> Metathesis ----
SVM_metathesis <- metathesis_train %>%
  encode(3) %>%
  svm(data = ., 
      norm_wt_mean ~ .,
      kernel = "radial",
      cost = 2)

metathesis_test_svm <- metathesis_test %>%
  encode(3) %>%
  mutate(predictions = predict(SVM_metathesis, select(., -c("norm_wt_mean")))) %>%
  mutate(r2 = cor(norm_wt_mean, predictions)^2)


#>>> Deallylation (indole) ----
SVM_RuIndole <- RuIndole_train %>%
  encode(3) %>%
  svm(data = ., 
      norm_wt_mean ~ .,
      kernel = "radial",
      cost = 1)

RuIndole_test_svm <- RuIndole_test %>%
  encode(3) %>%
  mutate(predictions = predict(SVM_RuIndole, select(., -c("norm_wt_mean")))) %>%
  mutate(r2 = cor(norm_wt_mean, predictions)^2)


#>>> Deallylation (coumarin) ----
SVM_coumarin <- coumarin_train %>%
  encode(1) %>%
  svm(data = ., 
      norm_wt_mean ~ .,
      kernel = "radial",
      cost = 8)

coumarin_test_svm <- coumarin_test %>%
  encode(1) %>%
  mutate(predictions = predict(SVM_coumarin, select(., -c("norm_wt_mean")))) %>%
  mutate(r2 = cor(norm_wt_mean, predictions)^2)


#>>> Hydroamination ----
SVM_AuIndole <- AuIndole_train %>%
  encode(4) %>%
  svm(data = ., 
      norm_wt_mean ~ .,
      kernel = "radial",
      cost = 4)

AuIndole_test_svm <- AuIndole_test %>%
  encode(4) %>%
  mutate(predictions = predict(SVM_AuIndole, select(., -c("norm_wt_mean")))) %>%
  mutate(r2 = cor(norm_wt_mean, predictions)^2)


#>>> Hydroarylation ----
SVM_AuFluo <- AuFluo_train %>%
  encode(3) %>%
  svm(data = ., 
      norm_wt_mean ~ .,
      kernel = "radial",
      cost = 4)

AuFluo_test_svm <- AuFluo_test %>%
  encode(3) %>%
  mutate(predictions = predict(SVM_AuFluo, select(., -c("norm_wt_mean")))) %>%
  mutate(r2 = cor(norm_wt_mean, predictions)^2)


#Artificial neural network ----

#Training data set
training_ann <- training %>%
  crossing(descriptor = c(1:5)) %>%
  mutate(train = map2(.x = train, .y = descriptor, ~encode(df = .x, descriptor = .y)),
         validate = map2(.x = validate, .y = descriptor, ~encode(df = .x, descriptor = .y))) %>%
  rowwise() %>%
  mutate(train_targets = list(train$norm_wt_mean),
         validate_targets = list(validate$norm_wt_mean),
         train = list(select(train, -c("norm_wt_mean"))),
         validate = list(select(validate, -c("norm_wt_mean"))))


#> Define model architectures ----
#(various numbers of hidden layers and nodes)

train_ANN_4x50 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 50, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_3x50 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 50, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_2x50 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 50, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_4x100 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 100, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_3x100 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 100, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_2x100 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 100, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_4x200 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 200, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_3x200 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 200, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}

train_ANN_2x200 <- function(train, train_targets, validate, validate_targets, batches = 10, epochs = 300){
  
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 200, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches,
        callbacks = c(early_stopping),
        validation_data = list(validate, validate_targets))
  
  return(model)
}



#> Hyperparameter tuning by cross-validation ----
MLP_combined <- training_ann %>%
  crossing(batches = c(10, 20)) %>%
  rowwise() %>%
  mutate(model_4X50 = list(train_ANN_4x50(train = as.matrix(from_list(.$train)),
                                          train_targets = as.matrix(from_list(.$train_targets)),
                                          validate = as.matrix(from_list(.$validate)),
                                          validate_targets = as.matrix(from_list(.$validate_targets)),
                                          batches = batches)),
         model_3x50 = list(train_ANN_3x50(train = as.matrix(from_list(.$train)),
                                          train_targets = as.matrix(from_list(.$train_targets)),
                                          validate = as.matrix(from_list(.$validate)),
                                          validate_targets = as.matrix(from_list(.$validate_targets)),
                                          batches = batches)),
         model_2x50 = list(train_ANN_2x50(train = as.matrix(from_list(.$train)),
                                          train_targets = as.matrix(from_list(.$train_targets)),
                                          validate = as.matrix(from_list(.$validate)),
                                          validate_targets = as.matrix(from_list(.$validate_targets)),
                                          batches = batches)),
         model_4x100 = list(train_ANN_4x100(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches)),
         model_3x100 = list(train_ANN_3x100(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches)),
         model_2x100 = list(train_ANN_2x100(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches)),
         model_4x200 = list(train_ANN_4x200(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches)),
         model_3x200 = list(train_ANN_3x200(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches)),
         model_2x200 = list(train_ANN_2x200(train = as.matrix(from_list(.$train)),
                                            train_targets = as.matrix(from_list(.$train_targets)),
                                            validate = as.matrix(from_list(.$validate)),
                                            validate_targets = as.matrix(from_list(.$validate_targets)),
                                            batches = batches))) %>%
  pivot_longer(cols = starts_with("model"), names_to = "model_architecture", values_to = "model") %>%
  rowwise() %>%
  mutate(val_loss_last = last(model$history$history$val_loss),
         val_loss_min = min(model$history$history$val_loss),
         epochs = which.min(model$history$history$val_loss)) %>%
  select(reaction, id, model_architecture, batches, descriptor, val_loss_last, val_loss_min, epochs)


#Best hyperparameters
MLP_combined %>%
  group_by(reaction, model_architecture, batches, descriptor) %>%
  summarise(mean_mse = mean(val_loss_min), mean_epochs = mean(epochs), .groups = "drop") %>%
  group_by(reaction) %>%
  slice_min(n = 1, order_by = mean_mse)


#> Train final models ----

#>>> Metathesis ----

#Training data set
Metathesis_train_MLP <- metathesis_train %>%
  encode(4)

#Function for model training
train_ANN_2x100_final <- function(train, train_targets, batches = 10, epochs = 300){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 100, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 100, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches)
  
  return(model)
}

#Train final model
Metathesis_MLP_final <- train_ANN_2x100_final(train = as.matrix(select(Metathesis_train_MLP, -c("norm_wt_mean"))),
                                              train_targets = as.matrix(Metathesis_train_MLP$norm_wt_mean),
                                              batches = 10,
                                              epochs = 21)

#Test model
Metathesis_test_MLP <- metathesis_test %>%
  encode(4) %>%
  mutate(predictions = predict(Metathesis_MLP_final, as.matrix(select(., -c("norm_wt_mean")))))

cor(Metathesis_test_MLP$predictions, Metathesis_test_MLP$norm_wt_mean)^2


#>>> Deallylation (indole) ----

#Training data set
RuIndole_train_MLP <- RuIndole_train %>%
  encode(1)

#Function for training final model
train_ANN_3x50_final <- function(train, train_targets, batches = 10, epochs = 300){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 50, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 50, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches)
  
  return(model)
}

#Train final model
RuIndole_MLP_final <- train_ANN_3x50_final(train = as.matrix(select(RuIndole_train_MLP, -c("norm_wt_mean"))),
                                                 train_targets = as.matrix(RuIndole_train_MLP$norm_wt_mean),
                                                 batches = 20,
                                                 epochs = 42)

#Test model
RuIndole_test_MLP <- RuIndole_test %>%
  encode(1) %>%
  mutate(predictions = predict(RuIndole_MLP_final, as.matrix(select(., -c("norm_wt_mean")))))

cor(RuIndole_test_MLP$predictions, RuIndole_test_MLP$norm_wt_mean)^2


#>>> Deallylation (coumarin) ----

#Training data set
coumarin_train_MLP <- coumarin_train %>%
  encode(4)

#Function for training final model
train_ANN_3x200_final <- function(train, train_targets, batches = 10, epochs = 300){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 200, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches)
  
  return(model)
}

#Train final model
coumarin_MLP_final <- train_ANN_3x200_final(train = as.matrix(select(coumarin_train_MLP, -c("norm_wt_mean"))),
                                            train_targets = as.matrix(coumarin_train_MLP$norm_wt_mean),
                                            batches = 20,
                                            epochs = 37)

#Test model
coumarin_test_MLP <- coumarin_test %>%
  encode(4) %>%
  mutate(predictions = predict(coumarin_MLP_final, as.matrix(select(., -c("norm_wt_mean")))))

cor(coumarin_test_MLP$predictions, coumarin_test_MLP$norm_wt_mean)^2


#>>> Hydroamination ----

#Training data set
AuIndole_train_MLP <- AuIndole_train %>%
  encode(1)

#Train final model
AuIndole_MLP_final <- train_ANN_3x200_final(train = as.matrix(select(AuIndole_train_MLP, -c("norm_wt_mean"))),
                                            train_targets = as.matrix(AuIndole_train_MLP$norm_wt_mean),
                                            batches = 10,
                                            epochs = 40)

#Test model
AuIndole_test_MLP <- AuIndole_test %>%
  encode(1) %>%
  mutate(predictions = predict(AuIndole_MLP_final, as.matrix(select(., -c("norm_wt_mean")))))

cor(AuIndole_test_MLP$predictions, AuIndole_test_MLP$norm_wt_mean)^2



#>>> Hydroarylation ----

#Training data set
AuFluo_train_MLP <- AuFluo_train %>%
  encode(2)

#Function for training the final model
train_ANN_4x200_final <- function(train, train_targets, batches = 10, epochs = 300){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 200, activation = "relu",
                input_shape = dim(from_list(train))[[2]]) %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1) %>% 
    compile(
      optimizer = "adam",
      loss = "mse")
  
  
  hist <- model %>%
    fit(x = train, 
        y = train_targets,
        epochs = epochs, 
        batch_size = batches)
  
  return(model)
}

#Train final model
AuFluo_MLP_final <- train_ANN_4x200_final(train = as.matrix(select(AuFluo_train_MLP, -c("norm_wt_mean"))),
                                          train_targets = as.matrix(AuFluo_train_MLP$norm_wt_mean),
                                          batches = 20,
                                          epochs = 22)

#Test model
AuFluo_test_MLP <- AuFluo_test %>%
  encode(2) %>%
  mutate(predictions = predict(AuFluo_MLP_final, as.matrix(select(., -c("norm_wt_mean")))))

cor(AuFluo_test_MLP$predictions, AuFluo_test_MLP$norm_wt_mean)^2


#Plot results ----

#Summary of model performances
ML_sd <- data.frame(method = c("Gradient boosting", "Neural network", "SVM"),
                    Metathesis = c(0.02, 0.02, 0),
                    Deallylation_coumarin = c(0.02, 0.03, 0),
                    Deallylation_indole = c(0.01, 0.05, 0),
                    Hydroamination = c(0.01, 0.03, 0),
                    Hydroarylation = c(0.02, 0.03, 0)) %>%
  pivot_longer(cols = -c("method"), names_to = "reaction", values_to = "sd")

ML_r2 <- data.frame(method = c("Gradient boosting", "Neural network", "SVM"),
                    Metathesis = c(0.68, 0.56, 0.70),
                    Deallylation_coumarin = c(0.68, 0.66, 0.67),
                    Deallylation_indole = c(0.71, 0.53, 0.63),
                    Hydroamination = c(0.46, 0.46, 0.45),
                    Hydroarylation = c(0.55, 0.34, 0.49)) %>%
  pivot_longer(cols = -c("method"), names_to = "reaction", values_to = "R2") %>%
  left_join(ML_sd, by = c("method", "reaction")) %>%
  mutate(lower = R2 - sd,
         upper = R2 + sd,
         method = factor(method, levels = c("SVM", "Gradient boosting", "Neural network")),
         reaction = factor(reaction),
         reaction = fct_relevel(reaction, "Metathesis"),
         reaction = fct_recode(reaction, "Deallylation\n(coumarin)" = "Deallylation_coumarin",
                               "Deallylation\n(indole)" = "Deallylation_indole"))

#Plot performance on test set

ML_r2 %>%
  ggplot(aes(reaction, R2, fill = method)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.8), width = 0.4) +
  scale_fill_manual(values = colorspace::sequential_hcl(palette = "GnBu", n = 5)) +
  theme_classic() +
  labs(y = expression(paste("R"^"2", " on test set")), x = "Reaction", fill = NULL) +
  theme(text = element_text(size = 7, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 6, r = 0, b = 0, l = 0, unit = "pt")),
        legend.text = element_text(size = 7, face = "plain"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black"), 
        plot.margin = margin(5,0,0,0,"pt"),
        legend.position = "bottom",
        legend.margin = margin(-5,0,0,0,"pt")) +
  scale_y_continuous(limits = c(0,0.8), expand = c(0,0)) 


#Plot predictions versus measured values - Gradient boosting, deallylation (indole)
#Note: Algorithm is stochastic and results will vary upon re-run
RuIndole_xgb_results %>%
  arrange(desc(prediction)) %>%
  mutate(top10 = c(rep(TRUE, 8), rep(FALSE, 72))) %>%
  ggplot(aes(norm_wt_mean, prediction, col = top10)) +
  geom_point(alpha = 0.6, shape = 16, size = 2) +
  labs(x = "Measured relative cell-specific activity", y = "Predicted relative cell-specific activity") +
  scale_x_continuous(limits = c(0, 7), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 6), expand = expansion(mult = c(0, 0))) +
  scale_colour_discrete(type = c("black", "green4")) +
  theme_bw() +
  theme(text = element_text(size = 7, colour = "black", family = "sans"),
        axis.text = element_text(size = 7, colour = "black", family = "sans"),
        strip.text = element_text(size = 7, colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(margin = margin(0,5,0,0, unit = "pt")),
        axis.title.x = element_text(margin = margin(3,0,0,0, unit = "pt")),
        legend.position = "none",
        plot.margin = margin(5,5,0,0,"pt"))
