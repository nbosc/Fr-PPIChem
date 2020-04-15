### 1
###################################################
### goal: optimization of QSAR model using SVM method
### ippi: iPPIDBv2 + Timbal
### nonIppi: bindingDB + BDM 
### desriptors: selection of 167 significants descriptors made on MOE/RDKit/Dragon
### variables are already normalized
###################################################
library(class)
library(e1071)
library(caret)
library(doMC)
source("/home/nicolasb/Scripts/roc.R")

# working directory
setwd("/home/nicolasb/Documents/PPI-Chem/Modeling/MOE+RDKIT+DRAGON")



# training set - 70% of the data - normalized in Knime
training <- read.csv(file="/home/nicolasb/data/model_input/IPPI_nonIPPI_variable_selection_from_moe_trainingSet.csv", row.names = 1, check.names = FALSE)
# test set - 30% of the data - normalized in Knime
test <- read.csv(file="/home/nicolasb/data/model_input/IPPI_nonIPPI_variable_selection_from_moe_testSet.csv", row.names = 1, check.names = FALSE)

start.time <- Sys.time() # starting time

#required for parallelization  
registerDoMC(cores = 5) # use 5 cpu cores

# settings for the training
settings <- trainControl(method="cv", number=5, repeats=1, returnData=FALSE, returnResamp="none", allowParallel=TRUE, index=createFolds(training$provider, 5), savePredictions=TRUE, summaryFunction=roc)

method <- "svmRadial"
  
# Grid generation
gammas <- c(0.0001,0.001,0.01,0.1,1)
costs <- c(0.1,1,10,100,1000,10000)
tune_grid <- expand.grid(.sigma=gammas, .C=costs)
  
# Modeling
tuning <- train(provider ~ ., data = training, method=method, trControl=settings, tuneGrid=tune_grid, metric="F1")
tuning # check if the method selects the best parameters

end.time <- Sys.time() # ending time
time.taken <- end.time - start.time
time.taken ### 
  
write.table(tuning$results, file="SVM_parameter_optimization_radial.csv", row.names=F, sep=",", quote=T)

# if the method has selected the best parameters
saveRDS(tuning$finalModel, "/home/nicolas/data/model_input/SVM_model_radial.rds")

# otherwise, build a new model with your own gamma and cost
#modelSVM <- svm(provider ~ ., data=training, type="C-classification", kernel="radial", gamma=0.01, cost=100)
#saveRDS(modelSVM, "/home/nicolas/data/model_input/SVM_model_radial.rds")

### check the inner model predictivity

start.time <- Sys.time() # starting time

train.predictions <- as.vector(predict(modelSVM, newdata = training))
internal.confusion.Matrix <- confusionMatrix(train.predictions, training$provider, positive = "IPPI")
internal.confusion.Matrix


test.predictions <- as.vector(predict(modelSVM, newdata = test))
int.confusion.Matrix <- confusionMatrix(test.predictions, test$provider, positive = "IPPI")
int.confusion.Matrix

end.time <- Sys.time() # ending time
time.taken <- end.time - start.time




# 
#
### 2
###################################################
### goal: performance estimation of SVM model build on diversity data
### ippi: iPPIDBv2 + Timbal
### nonIppi: bindingDB + BDM
### /!\ only the diverse molecule from ippi and nonIppi were kept
### descriptors: selection of 167 significants descriptors made on MOE/RDKit/Dragon
### variables are already normalized
###################################################
library(class)
library(e1071)
library(caret)
library(doMC)
source("/home/nicolasb/Scripts/roc.R")


setwd("/home/nicolasb/Documents/PPI-Chem/Modeling/MOE+RDKIT+DRAGON/diverse")


### training data
training <- read.csv(file="/home/nicolasb/data/model_input/IPPI_nonIPPI_diverse_variable_selection_from_moe_rdkit_dragon_trainingSet.csv", row.names = 1)

#test <- read.csv(file="/home/nicolasb/data/model_input/IPPI_nonIPPI_diverse_variable_selection_from_moe_rdkit_dragon_testSet.csv", row.names = 1)

#val <- read.csv(file="/home/nicolasb/data/model_input/ChemBridge_variable_selection_from_moe_rdkit_dragon_validationSet.csv", row.names = 1)


# find the best parameters

start.time <- Sys.time() # starting time

registerDoMC(cores = 5) # use 5 cpu cores

settings <- trainControl(method="cv", number=5, repeats=1, returnData=FALSE, returnResamp="none", allowParallel=TRUE, index=createFolds(training$provider, 5), savePredictions=TRUE, summaryFunction=roc)

method <- "svmRadial"

# Grid1
gammas <- c(0.0002, 0.0004, 0.0006, 0.0008, 0.001)
costs <- c(0.1, 1, 10, 100, 1000, 10000)
tune_grid <- expand.grid(.sigma=gammas, .C=costs)

# Model
tuning <- train(provider ~ ., data = training, method=method, trControl=settings, tuneGrid=tune_grid, metric="F1")


end.time <- Sys.time() # ending time
time.taken <- end.time - start.time
time.taken ### min

write.table(tuning$results, file="SVM_parameters_optimization_diverse2.csv", row.names=F, sep=",", quote=T)

### check the inner model predictivity

# 
# start.time <- Sys.time() # starting time
# 
# modelSVM <- svm(provider ~ ., data=data, type="C-classification", kernel="radial", gamma=0.01, cost=100)
# 
# train.predictions <- as.vector(predict(modelSVM, newdata = data))
# internal.confusion.Matrix <- confusionMatrix(train.predictions, data$provider, positive = "IPPI")
# internal.confusion.Matrix
# 
# 
# test.predictions <- as.vector(predict(modelSVM, newdata = test))
# int.confusion.Matrix <- confusionMatrix(test.predictions, test$provider, positive = "IPPI")
# int.confusion.Matrix
# 
# val.predictions <- as.vector(predict(modelSVM, newdata = val))
# sensitivityValidation <- length(grep("\\<IPPI\\>", val.predictions))/length(val.predictions)
# 
# resPred <- data.frame(val,val.predictions) 
# write.csv(resPred, file="/home/nicolasb/data/model_output/ChemBridge_variable_selection_from_moe_rdkit_dragon_validationSetPredictions.csv")

#	Reference
#Prediction active inactive
#  active   
#  inactive

#Sensitivity :  
#Specificity : 

#end.time <- Sys.time() # ending time
#time.taken <- end.time - start.time

