#' Remove zeros
#'
#' This function removes events with 0 expression of mandatory makers (e.g. DNA 
#' intercalator) as well as cells with 0 expression in all optional markers 
#' (e.g. CD4, CD8, EpCAM, CD31, etc.). Mandatory markers are those which should 
#' be expressed for an event to be considered valid. Optional markers are those 
#' that allow cell type identification, therefore their absolute abscence make 
#' the event invalid because it cannot be identified as a specific cell type.
#'
#' @param df An object of class \code{data.frame}.
#' @param mand_inx A vector with the column indexes for the mandadory markers.
#' @param opt_idx A vector with the column indexes for the optional markers.
#'
#' @return Returns a \code{data.frame} with the invalid events already removed.
#' @export rm_zeros()
rm_zeros <- function (df, mand_idx, opt_idx, ...){
        mand <- c()
        for (i in 1:length(mand_idx)){
            mand <- c(mand,(which(df[,mand_idx[i]] == 0)))
        }
        op <- c()
        for (i in 1:length(opt_idx)){
            op <- c(op,(which(df[,opt_idx[i]] == 0)))
        }
        # Select only the events == 0 for all 'optional' markers
        optional_combined <- as.data.frame(table(op))
        op <- as.numeric(subset(optional_combined, Freq == length(opt_idx), 
            select = op)$op)
        k <- unique(c(mand, op))
        if (length(k) == 0) {
            df <- df
        } else {
            df <- df[-k, ]
        }
        #cat('# of rows dropped: ', length(k), '\n')

    return(df)
}



#' Obtain a class-balanced sample
#'
#' This function samples the rows of a \code{data.frame} with balanced classes.
#' Useful when original training set has severely unbalanced classes.
#'
#' @param df An object of class \code{data.frame}.
#' @param sample_size Numeric. Size of the sample from each class.
#' @param class_col A character vector with the name of the column that identify 
#'        the classes.
#'
#' @return Returns a \code{data.frame} with n=\code{sample_size} rows per class
#'         randomly sampled.
#' @export BalancedSample()
BalancedSample <- function(df, sample_size = 5000, class_col = class_col){
    idx <- c(sample(which(df[class_col] == 1), sample_size, replace = T), #noise/beads
             sample(which(df[class_col] == 0), sample_size, replace = T)) #cells
    df <- df[idx,]
    return(df)
}




#' Get a list of dataframes for training
#'
#' This function takes a list of \code{data.frame} and returns a random sample 
#' of the original list to be used for for training a model. The remaining 
#' \code{data.frame} should be used to test the model.
#' 
#' @param dt_ls A \code{list} of \code{data.frame}s
#' @param s_train Numeric. The size of the training set. Any value between (0:1)
#' 
#' @return A random subset of \code{data.frame}s as a \code{list}.
#' @export TrainSetList()
TrainSetList <- function(dt_ls, s_train = 0.75) {
    train_size <- round(length(dt_ls)*s_train)
    dt_train <- sample(dt_ls, train_size)
    return(dt_train)
}



#' Generates datasets for model training
#'
#' This function generates training and test sets for model training with 
#' balanced classes using functions \code{BalancedSample} and \code{TrainSetList}.
#' The training and test datasets as well as vectors that identify the file of 
#' origin of each row are saved as an RData object.
#'
#' @param dt_ls A \code{list} of \code{data.frame}s with balanced classes.
#' @param output_path A character vector with full path names.
#' @param label A character vector with the prefix of the RData file to be outputted.
#' @param class_col A character vector with the name of the column that identify 
#'        the classes.
#' @export TrainTest()
TrainTest <- function(dt_ls, output_path = output_path, label = label, 
    class_col = class_col){
    dt_train <- TrainSetList(dt_ls)
    k <- which(names(dt_ls)%in% names(dt_train))
    dt_test <- dt_ls[-k]

    class_size <- table(dt_ls[[1]][class_col])[1]

    # Save names of files used for training and test for reference
    train_nms <- rep(names(dt_train), each=class_size*2)
    test_nms <- rep(names(dt_test), each=class_size*2)

    train_set <- do.call("rbind", c(dt_train, make.row.names = FALSE))
    test_set <- do.call("rbind", c(dt_test, make.row.names = FALSE))
    test_train_l <- list('train_set' = train_set, 'test_set' = test_set)
    
    save(train_set, test_set, train_nms, test_nms, file = file.path(output_path, 
        paste0(label, '_TrainTest.RData')))
}



#' Converts the class labels from numeric to factor.
#'
#' A function to make the class column compatible with \code{caret} training function
#' for classification. Should be used when classes are labeled as numeric and
#' need to be transform to a factor.
#' 
#' @param df An object of class \code{data.frame}.
#' @param class_col A character vector with the name of the column that identify 
#'        the classes.
#' @param name_0 A character vector with the name for the class==0
#' @param name_1 A character vector with the name for the class==1
#'
#' @return Returns a \code{data.frame} with the original class labels converted
#'         into factors.
#' @export ClassVarToFactor()
ClassVarToFactor <- function(df, class_col = class_col, name_0 = name_0, 
  name_1 = name_1) {
  idx_0 <- which(df[class_col] == 0)
  idx_1 <- which(df[class_col] == 1)

  df[idx_0, class_col] <- name_0
  df[idx_1, class_col] <- name_1

  return(df)
}





#' Tunes and trains a classification model
#'
#' This function tunes and train a model (or models) for classification. Current
#' version supports Random Forest and XGBoost algorithms. Allows parallelization.
#' The final model, feature importance, prediction results of the test dataset 
#' and confusion matrix are outputted as an RData object.
#'
#' @param TrainSet A \code{data.frame} with balanced classes.
#' @param TestSet A \code{data.frame}.
#' @param alg A character vector with the name of the classification algorithm 
#'        to be used. Options are Random Forest 'RF', XGBoost 'XGB' or both 'all'.
#' @param class_col A character vector with the name of the column that identify 
#'        the classes.
#' @param seed A numeric vector to set a seed for reproducible results.
#' @param name_0 A character vector with the name for the class==0
#' @param name_1 A character vector with the name for the class==1
#' @param label A character vector with the prefix of the RData file to be outputted.
#' @param allowParallel Logical. If \code{TRUE} allows parallel computation.
#' @param free_cores A numeric vector with the number of cores to be left free 
#'        if \code{allowParallel=TRUE}'.
#' @export TrainModel()
TrainModel <- function(TrainSet, TestSet, alg = c('all', 'RF', 'XGB'), 
    class_col = class_col, seed = 40, name_0 = name_0, name_1 = name_1,
    label = label, allowParallel = TRUE, free_cores = 4){
    # name_0 == 'cells' |  name_1 == 'debris'/'beads'/'dead'
    # Set seed
    set.seed(seed)    
    
    alg <- match.arg(alg, c('all', 'RF', 'XGB'))

    # Ask for channels for training features   
    print(as.matrix(colnames(TrainSet)))
    prompt <- "Enter the column INDICES of the training features (separated by 
              single space only, no comas allowed) \n"
    features <- as.numeric(strsplit(readline(prompt), " ")[[1]])
    features <- colnames(TrainSet)[features]
    features <- features[order(features)]

    # Classes to names
    TrainSet <- ClassVarToFactor(TrainSet, class_col = class_col, 
        name_0 = name_0, name_1 = name_1)
    TestSet <- ClassVarToFactor(TestSet, class_col = class_col, 
        name_0 = name_0, name_1 = name_1)

    TrainSet <- t_asinh(TrainSet[c(class_col, features)])
    TestSet <- t_asinh(TestSet[c(class_col, features)])

    colnames(TrainSet) <- c(class_col, features)

    colnames(TestSet) <- colnames(TrainSet)



    if(alg == 'RF' || alg == 'all'){
        print('Running Random Forest')
        # Initiating parallelization
        if(allowParallel){
        cluster <- parallel::makeCluster(parallel::detectCores() - free_cores) 
        doParallel::registerDoParallel(cluster)
        }
        # Computing training control parameters
        fitControl <- caret::trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 3,
            search = 'grid',
            savePredictions = 'final',
            classProbs = T,
            summaryFunction = twoClassSummary,
            allowParallel = allowParallel,
            returnData = FALSE)

        # Setting tune grid
        tunegrid <- expand.grid(.mtry=c(1:length(features)))

        # Fitting the model
        model_rf <- caret::train(TrainSet[,-1], TrainSet[,1],
            method ='rf',
            trControl = fitControl,
            metric = 'ROC',
            tuneGrid=tunegrid)

        # Assesing feature importance
        ftimp_rf <- caret::varImp(model_rf)

        # Test
        # Prediction on test dataset
        pred_rf <- predict(model_rf, TestSet)
        
        # Compute confusion matrix
        conf_rf <- caret::confusionMatrix(
            reference = as.factor(TestSet[,1]),
            data = pred_rf,
            mode = 'everything',
            positive = name_1)

        #save RData
        save(model_rf, ftimp_rf, pred_rf, conf_rf, file = paste0(label, '_RFmodel.RData'))
        print('Random Forest completed')
        
        # Finalizing parallelization
        if(allowParallel){
        parallel::stopCluster(cluster)
        foreach::registerDoSEQ()
        }
    }
    

    if(alg == 'XGB' || alg == 'all'){
        print('Running XGBoost')
        # Rearrangement
        X_train = xgboost::xgb.DMatrix(as.matrix(TrainSet[,-1]))
        y_train = TrainSet[,1]

        # Computing training control parameters
        xgb_trcontrol = caret::trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 3,
            search = 'grid',
            savePredictions = 'final',
            summaryFunction = twoClassSummary,
            allowParallel = allowParallel,
            classProbs = TRUE,
            verboseIter = FALSE,
            returnData = FALSE)

        # Specifying grid space
        xgbGrid <- expand.grid(nrounds = c(100,200),  # this is n_estimators in the python code above
                               max_depth = c(10, 15, 20, 25),
                               colsample_bytree = seq(0.5, 0.9, length.out = 5),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.1,
                               gamma=0,
                               min_child_weight = 1,
                               subsample = 1)

        # Fitting the model
        model_xgb = caret::train(
            X_train, y_train,
            trControl = xgb_trcontrol,
            tuneGrid = xgbGrid,
            method = "xgbTree",
            metric = 'ROC',
            nthread = detectCores() - free_cores)

        # Assesing feature importance
        ftimp_xgb <- caret::varImp(model_xgb)
        
        # Test
        X_test = xgboost::xgb.DMatrix(as.matrix(TestSet[,-1]))
        y_test = TestSet[,1]

        # Prediction on test dataset
        pred_xgb <- predict(model_xgb, X_test)
        # Compute confusion matrix
        conf_xgb <- caret::confusionMatrix(
            reference = as.factor(y_test),
            data = pred_xgb,
            mode = 'everything',
            positive = name_1)

        # Save RData
        save(model_xgb, ftimp_xgb, pred_xgb, conf_xgb, file = paste0(label, '_XGBmodel.RData'))
        print('XGBoost completed')
    }

}




#' Predict class
#'
#' This function predicts the class of the events (rows) based on a classification 
#' model previously trained and returns the row indexes of the positive class.
#'
#' @param df An object of class \code{data.frame}.
#' @param model Object of class \code{train} that contains the model to be used 
#'        for classification of the new data.
#' @param alg A character vector with the name of the classification algorithm 
#'        used to train \code{model}.
#' @param features A character vector with the features used to train \code{model}
#'        (also known as predictors).
#' @param label A character vector with the name of the positive class. (e.g. 'beads')
#'
#' @return Returns a numeric vector with the row indexes of the positive class.
#' @export predict_cl()
predict_cl <- function(df, model = model, alg = c('RF', 'XGB'), 
  features = features, label = label){
    
    alg <- match.arg(alg, c('RF', 'XGB'))

    features <- features[order(features)]

    df <- t_asinh(df[features])
    
    if(alg == 'XGB'){
        df <- xgboost::xgb.DMatrix(as.matrix(df))

    }

    pred <- predict(model, df)

    k <- which(pred == label)

    return(k)
}