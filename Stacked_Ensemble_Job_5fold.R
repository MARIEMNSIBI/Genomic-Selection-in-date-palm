library(readxl)
library(data.table)
library(BGLR)
library(rrBLUP)
library(randomForest)
library(xgboost)

Pheno <- read_excel("Pheno_BLUEs.xlsx")
Pheno$Variety<-as.factor(Pheno$Variety)
Pheno$FL<-as.numeric(as.character(Pheno$FL))
Pheno$R<-as.numeric(as.character(Pheno$R))
Pheno$G<-as.numeric(as.character(Pheno$G))
Pheno$B<-as.numeric(as.character(Pheno$B))
Pheno$L<-as.numeric(as.character(Pheno$L))
Pheno$a<-as.numeric(as.character(Pheno$a))
Pheno$b<-as.numeric(as.character(Pheno$b))
Pheno$RAB<-as.numeric(as.character(Pheno$RAB))
Pheno$Hue<-as.numeric(as.character(Pheno$Hue))
Pheno$ANT<-as.numeric(as.character(Pheno$ANT))
Pheno$MO<-as.numeric(as.character(Pheno$MO))
Pheno$SUC<-as.numeric(as.character(Pheno$SUC))
Pheno$GLU<-as.numeric(as.character(Pheno$GLU))
Pheno$FRU<-as.numeric(as.character(Pheno$FRU))
Pheno$RSGFS<-as.numeric(as.character(Pheno$RSGFS))
Pheno$OXA<-as.numeric(as.character(Pheno$OXA))
Pheno$CIT<-as.numeric(as.character(Pheno$CIT))
Pheno$TAR<-as.numeric(as.character(Pheno$TAR))
Pheno$MAL<-as.numeric(as.character(Pheno$MAL))
Pheno$SUI<-as.numeric(as.character(Pheno$SUI))
Pheno$FUM<-as.numeric(as.character(Pheno$FUM))


data<-fread("SNPs_sample_filtered_INPUT.csv")
data<-as.matrix(data)
Markers<-colnames(data)[2:1812002]
Individuals<- data[,1]
data <- data[,-1]
data <-apply(data, 2,as.numeric) #Numeric-matrix
rownames(data) <- Individuals

data[data==0] <- -1
data[data==1] <- 0
data[data==2] <- 1
geno <- data
meanimputed <- function(X){
  colimputed <- function(x){
    x[is.na(x)] <- mean(!is.na(x))
    return(x)
  }
  Ximp <- apply(X,2,colimputed)
  return(as.matrix(Ximp))
}
genoimputed <- meanimputed(geno)
data <- genoimputed

Pheno<-as.matrix(Pheno)
row.names(Pheno)<-Pheno[,1]

Pheno<-Pheno[,-1]
Individuals<-rownames(Pheno)
Pheno <-apply(Pheno, 2,as.numeric)
row.names(Pheno)<-Individuals

############################################################################
#                   Genomic prediction model: Stacked Ensemble             #
############################################################################

data <- as.matrix(data)
mode(data) <- "numeric"

traits = 1
cycles = 50

accuracy = matrix(nrow = cycles, ncol = traits)
mse = matrix(nrow = cycles, ncol = traits)
bias = matrix(nrow = cycles, ncol = traits)
rel_mse = matrix(nrow = cycles, ncol = traits)
bias_std = matrix(nrow = cycles, ncol = traits)
mse_std = matrix(nrow = cycles, ncol = traits)

# PCA for RF and XGBoost
data_pca <- prcomp(data)$x[, 1:50]

# Kernel matrix for RKHS
data_scaled <- scale(data, center=TRUE,scale=FALSE)
kin <- tcrossprod(data_scaled) / ncol(data)
# XGBoost grid
xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = 3,
  eta = c(0.01, 0.1, 0.3),
  gamma = 0,
  colsample_bytree = 0.6,
  min_child_weight = 1,
  subsample = c(0.6, 0.8)
)

system.time(
for (r in 1:cycles) {
  cat("Cycle", r, "\n")
   num_folds <- 5
   folds <- sample(rep(1:num_folds, length.out = nrow(Pheno)))

  fold_accuracy <- numeric(num_folds)
  fold_mse <- numeric(num_folds)
  fold_bias <- numeric(num_folds)
  fold_rel_mse <- numeric(num_folds)
  fold_bias_std <- numeric(num_folds)
  fold_mse_std <- numeric(num_folds)
  
  for (fold in 1:num_folds) {
    train <- which(folds != fold)
    test  <- which(folds == fold)
    
    Pheno_train <- Pheno[train, ]
    Pheno_valid <- Pheno[test, ]
    y_valid <- Pheno_valid[, 1]
    
    m_train <- as.matrix(data[train, ])
    m_valid <- as.matrix(data[test, ])
    mode(m_train) <- "numeric"
    mode(m_valid) <- "numeric"
    
    m_train_pca <- data_pca[train, ]
    m_valid_pca <- data_pca[test, ]
    
    # BayesB
    yNA = rep(NA, nrow(Pheno))
    yNA[train] = Pheno_train[, 1]
    fm_bayesB = BGLR(y = yNA, ETA = list(list(X = data, model = "BayesB")),
                     verbose = FALSE, nIter = 1500, burnIn = 500)
    pred_bayesB = fm_bayesB$yHat[test]
    
    # BayesC
    fm_bayesC = BGLR(y = yNA, ETA = list(list(X = data, model = "BayesC")),
                     verbose = FALSE, nIter = 1500, burnIn = 500)
    pred_bayesC = fm_bayesC$yHat[test]
    
    # BL
    fm_bl = BGLR(y = yNA, ETA = list(list(X = data, model = "BL")),
                 verbose = FALSE, nIter = 1500, burnIn = 500)
    pred_bl = fm_bl$yHat[test]
    
    # RF
    rf_model = randomForest(Pheno_train[, 1] ~ ., data = as.data.frame(m_train_pca), ntree = 100)
    pred_rf = predict(rf_model, newdata = as.data.frame(m_valid_pca))
    
    # rrBLUP
    rr_model = mixed.solve(y = Pheno_train[, 1], Z = m_train)
    u_vec = rr_model$u
    
    if (ncol(m_valid) != length(u_vec)) {
      min_dim <- min(ncol(m_valid), length(u_vec))
      m_valid_adj <- m_valid[, 1:min_dim, drop = FALSE]
      u_vec_adj <- u_vec[1:min_dim]
    } else {
      m_valid_adj <- m_valid
      u_vec_adj <- u_vec
    }
    
    pred_rr = as.vector(m_valid_adj %*% u_vec_adj) + rr_model$beta
    
    # XGBoost
    best_rmse = Inf
    best_pred_xgb = rep(NA, length(test))
    for (i in 1:nrow(xgb_grid)) {
      params <- as.list(xgb_grid[i, ])
      params$objective <- "reg:squarederror"
      params$eval_metric <- "rmse"
      
      dtrain <- xgb.DMatrix(data = as.matrix(m_train_pca), label = Pheno_train[, 1])
      dtest <- xgb.DMatrix(data = as.matrix(m_valid_pca), label = y_valid)
      
      xgb_model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = params$nrounds,
        watchlist = list(eval = dtest),
        early_stopping_rounds = 10,
        verbose = 0
      )
      
      pred_xgb = predict(xgb_model, newdata = dtest)
      rmse_score = sqrt(mean((pred_xgb - y_valid)^2))
      
      if (rmse_score < best_rmse) {
        best_rmse = rmse_score
        best_pred_xgb = pred_xgb
      }
    }
    
    # RKHS
    yNA_rkhs = rep(NA, nrow(Pheno))
    yNA_rkhs[train] = Pheno_train[, 1]
    fm_rkhs = BGLR(y = yNA_rkhs, ETA = list(list(K = kin, model = "RKHS")),
                   verbose = FALSE, nIter = 1500, burnIn = 500)
    pred_rkhs = fm_rkhs$yHat[test]
    
    # Meta-model
    pred_df = data.frame(
      BayesB = pred_bayesB,
      BayesC = pred_bayesC,
      BL = pred_bl,
      RF = pred_rf,
      rrBLUP = pred_rr,
      XGB = best_pred_xgb,
      RKHS = pred_rkhs
    )
    
    meta_model = lm(y_valid ~ ., data = pred_df)
    final_pred = predict(meta_model, newdata = pred_df)
    
    # Evaluation
    var_y = var(y_valid)
    sd_y = sd(y_valid)
    
    fold_accuracy[fold] = cor(final_pred, y_valid, use = "complete.obs")
    fold_mse[fold] = mean((y_valid - final_pred)^2)
    fold_bias[fold] = mean(y_valid - final_pred)
    fold_rel_mse[fold] = fold_mse[fold] / var_y
    fold_bias_std[fold] = fold_bias[fold] / sd_y
    fold_mse_std[fold] = fold_mse[fold] / sd_y
  }
  
  # Average across folds per cycle
  accuracy[r, 1] = mean(fold_accuracy)
  mse[r, 1]      = mean(fold_mse)
  bias[r, 1]     = mean(fold_bias)
  rel_mse[r, 1]  = mean(fold_rel_mse)
  bias_std[r, 1] = mean(fold_bias_std)
  mse_std[r, 1]  = mean(fold_mse_std)
})

write.table(accuracy, "Ensemble_5fold_cv_Ajwa_sampled_filtered_accuracy_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
write.table(mse, "Ensemble_5fold_cv_Ajwa_sampled_filtered_mse_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
write.table(bias, "Ensemble_5fold_cv_Ajwa_sampled_filtered_bias_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
write.table(rel_mse, "Ensemble_5fold_cv_Ajwa_sampled_filtered_rel_mse_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
write.table(bias_std, "Ensemble_5fold_cv_Ajwa_sampled_filtered_bias_std_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
write.table(mse_std, "Ensemble_5fold_cv_Ajwa_sampled_filtered_mse_std_FL_50iterations.txt", sep=" ", quote=F, col.names=F, row.names=T)
