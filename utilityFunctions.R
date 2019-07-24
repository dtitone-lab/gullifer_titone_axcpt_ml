# Utility funcitons used by other R scripts.

# Preprocessing -----------------------------------------------------------
findNAs <- function(data){
  nas= as.data.frame(data) %>% summarise_all(funs(sum(is.na(.)))) %>%
    gather(column,countNA) %>% arrange(desc(countNA))
  return((nas))
}

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

runPCA <- function(data, print=F){
  entropy_data <- data%>%select(participant,contains("entropy"), contains("Entropy"))
  
  #fa.parallel(entropy_data %>% select(-participant))
  
  entropy.pca<-principal(entropy_data %>%select(-participant), 2, rotate = "oblimin")
  
  if(print){
    print(loadings(entropy.pca), cutoff=.3)
  }
  scores <- data.frame(entropy=entropy.pca$scores)
  
  return(entropy.pca)
}

predictPCA <- function(testdata, traindata, pca){
  entropy_test <- testdata %>% select(participant, contains("entropy"), contains("Entropy")) 
  entropy_train <- traindata %>% select(participant,contains("entropy"), contains("Entropy"))
  
  entropy_test_pca <- predict.psych(pca, data=entropy_test %>% select(-participant),
                                    old.data = entropy_train %>%select(-participant))
  
  scores <- data.frame(entropy=entropy_test_pca)
  return(scores)
}

preprocess_steps <- function(training_set, testing_set=NULL){
  
  
  training_sdata <- training_set %>% select(participant, contains("entropy")) %>% 
    group_by(participant) %>% summarise_all(first)
  
  training_pca <- runPCA(training_sdata)
  
  # Create test set for this iteration
  # Subset all the datapoints where .folds matches the current fold
  if(!is.null(testing_set)){
    test_pca_scores <- predictPCA(testing_set, training_sdata, training_pca)
  }
  
  
  training_sdata <- data.frame(training_sdata, entropy=training_pca$scores)
  training_sdata <- training_sdata %>% select(participant, PCA_General = entropy.TC1, PCA_Work = entropy.TC2)
  
  training_set <- left_join(training_set, training_sdata, by="participant")
  
  training_sdata <- training_set %>% select(participant, aoa, L2_exposure) %>% 
    group_by(participant) %>% summarise_all(first)
  
  scales <- build_scales(dataSet = training_sdata, 
                         cols = c("aoa", "L2_exposure"), verbose = F)
  
  training_set <- fastScale(dataSet = training_set, scales = scales, verbose = F)
  
  
  if(!is.null(testing_set)){
    
    
    
    testing_set <- data.frame(testing_set, test_pca_scores)
    
    testing_set  <- testing_set %>% rename(PCA_General = entropy.TC1, PCA_Work = entropy.TC2)
    
    testing_set  <- fastScale(dataSet = testing_set, scales = scales, verbose = F)
    
    
    return(list(training_set=training_set, testing_set=testing_set))
  } else{
    return(training_set)
  }
}

axcpt_preprocess <- function(axcpt){
  axcpt <- axcpt[axcpt$Condition!="AX",]
  axcpt$RT_correct    <- ifelse(axcpt$accuracy_target == 1 & axcpt$Target.RT > 40,
                                axcpt$Target.RT, NA)
  
  axcpt$acc <- factor(axcpt$accuracy_target)
  
  axcpt$Condition <- factor(axcpt$Condition)
  axcpt$Condition <- relevel(axcpt$Condition, ref="BX")
  
  
  return(axcpt)
}

# Cross-validation --------------------------------------------------------
crossvalidate <- function(data, k, mod, dependent, dv_continuous, random = FALSE, returnRuns=FALSE){
  # data is the training set with the ".folds" column
  # k is the number of folds we have
  # model is a string describing a linear regression model formula
  # dependent is a string with the name of the score column we want to predict
  # random is a logical; do we have random effects in the model?
  
  # Initialize empty list for recording performances
  if(isTRUE(dv_continuous)){
    ncol=7
  }else{
    ncol=5
  }
  
  data <- data %>% 
    select(-contains("PCA"))
  
  performances <- matrix(ncol = ncol, nrow=k)
  
  #Get rid of PCA columns from the entire dataset, they will be recalculated for each fold
  #data <- data %>% select(-PCA_General, -PCA_Work)
  
  # One iteration per fold
  for (fold in 1:k){
    
    if(fold %% 10 == 0){
      print(paste("current iteration:",fold,"/",k))
    }
    # Create training set for this iteration
    # Subset all the datapoints where .folds does not match the current fold
    training_set <- data[data$.folds != fold,]
    testing_set <- data[data$.folds == fold,]
    
    # preprocess the data for training set, and testing set based on training
    # set
    preproc_data <- preprocess_steps(training_set, testing_set)
    
    training_set <- preproc_data$training_set
    testing_set <- preproc_data$testing_set
    
    
    ## Train model
    
    # If there is a random effect,
    # use lmer() to train model
    # else use lm()
    
    if (isTRUE(random)){
      if(isTRUE(dv_continuous)){
        # Train linear mixed effects model on training set
        model <-  lmer(mod, training_set, REML=FALSE)
      } else{
        # Train gen. linear mixed effects model on training set
        model <-  glmer(mod, training_set,  family = "binomial",
                        glmerControl(optimizer ="bobyqa", optCtrl = list(maxfun=1000000)))
      }
    } else{
      if(isTRUE(dv_continuous)){
        # Train linear model on training set
        model <-  lm(mod, training_set)
      } else{
        # Train gen. linear model on training set
        model <-  glm(mod, training_set, family = "binomial")
      }
    }
    
    if(isTRUE(dv_continuous)){
      # Predict the dependent variable in the testing_set with the trained model
      predicted_train <- predict(model, training_set)
      predicted_test  <- predict(model, testing_set, allow.new.levels=TRUE)
      
      #rmse for train too
      
      # Get the Root Mean Square Error between the predicted and the observed
      RMSE_train <- rmse(predicted_train, training_set[[dependent]])
      RMSE_test  <- rmse(predicted_test, testing_set[[dependent]])
      
      MAE_train  <- mae(predicted_train, training_set[[dependent]])
      MAE_test   <- mae(predicted_test, testing_set[[dependent]])
      
      r2_train <- cor(predicted_train, training_set[[dependent]], use="pairwise.complete.obs")^2
      r2_test  <- cor(predicted_test, testing_set[[dependent]],   use="pairwise.complete.obs")^2
      
      performances_cur <- cbind(RMSE_train, MAE_train, RMSE_test, MAE_test, r2_train, r2_test, fold)
      
    }else{
      predicted_train <- predict(model, training_set, type="response")
      predicted_test  <- predict(model, testing_set, allow.new.levels=TRUE, type="response")
      
      roc_train <- roc(response=training_set[[dependent]],
                       predictor=predicted_train)
      
      roc_test <- roc(response=testing_set[[dependent]],
                      predictor=predicted_test)
      
      # e <- cbind(roc_train$thresholds,roc_train$sensitivities+roc_train$specificities)
      # opt_t <- subset(e,e[,2]==max(e[,2]))[,1]
      # 
      # print(paste("thresholding predicted resp at", opt_t))
      
      auc_train <- as.numeric(auc(roc_train))
      auc_test  <- as.numeric(auc(roc_test))
      
      predicted_train <- ifelse(predicted_train >= 0.5,
                                levels(training_set[[dependent]])[2],
                                levels(training_set[[dependent]])[1])
      
      predicted_test  <- ifelse(predicted_test  >= 0.5,
                                levels(testing_set[[dependent]])[2],
                                levels(testing_set[[dependent]])[1])
      
      error_train <- mean(predicted_train != training_set[[dependent]])
      error_test  <- mean(predicted_test != testing_set[[dependent]])
      
      performances_cur <- cbind(error_train, error_test, auc_train, auc_test, fold)
      
      
    }
    # Add the RMSE to the performance list
    performances[fold,] <- performances_cur
    
  }
  
  performances <- data.frame(performances)
  colnames(performances) <- colnames(performances_cur)
  
  performances_mean <- matrix(colMeans(performances), ncol=ncol)
  colnames(performances_mean) <- colnames(performances_cur)
  performances_mean <- as.data.frame(performances_mean)
  
  if(isTRUE(dv_continuous)){
    performances_mean$RMSE_train_sem <- sd(performances$RMSE_train) / sqrt(k)
    performances_mean$RMSE_test_sem  <- sd(performances$RMSE_test)  / sqrt(k)
    performances_mean$MAE_train_sem  <- sd(performances$MAE_train)  / sqrt(k)
    performances_mean$MAE_test_sem   <- sd(performances$MAE_test)   / sqrt(k)
  } else{
    performances_mean$error_train_sem <- sd(performances$error_train) / sqrt(k)
    performances_mean$error_test_sem  <- sd(performances$error_test)  / sqrt(k)
    performances_mean$auc_train_sem <- sd(performances$auc_train) / sqrt(k)
    performances_mean$auc_test_sem  <- sd(performances$auc_test)  / sqrt(k)
    
  }
  
  df = attr(logLik(model), "df")
  performances      <- cbind.data.frame(performances, mod, df=df)
  performances_mean <- cbind.data.frame(performances_mean, mod, df=df) 
  
  performances_mean <- performances_mean %>% select(-fold)
  # Return the mean of the recorded RMSEs
  
  if(isTRUE(returnRuns)){
    return(list(performances=performances,performances_mean=performances_mean))
  }
  else{
    return(performances_mean)
  }
}

print_cval <- function(df){
  
  df.min <- df %>% 
    filter(lambda==lambda.min) %>% 
    mutate(loglambda = log(lambda)) %>% 
    mutate_if(is.numeric, round, 2)
    
  df.1se <- df %>% 
    filter(lambda==lambda.1se) %>% 
    mutate(loglambda = log(lambda)) %>% 
    mutate_if(is.numeric, round, 2)
  
  with(df.min,
       cat(paste0("min log λ: ", loglambda, " error metric: ", cvm, " [", cvlo, ", ", cvup, "]\n")))

                  
                  
  with(df.1se,
       cat(paste0("1se log λ: ", loglambda, " error metric: ", cvm, " [", cvlo, ", ", cvup, "]\n")))
}

# Plotting ----------------------------------------------------------------


grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

