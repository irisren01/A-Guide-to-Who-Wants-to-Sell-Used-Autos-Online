
generate_indices = function(n,per) {
  return(sample.int(n,round(per*n)))
}

#############################################################################
linear_model = function(newdata,test_ind) {
  
  
  model = lm(newdata[-test_ind,])
  mse = mean((predict(model,newdata[test_ind,2:ncol(newdata)]) - newdata$biddy1[test_ind])^2)
  
  return(mse)
} 
#############################################################################
covariate_selection = function(newdata2,test_ind) {
  y = newdata2$biddy1
  covariates = newdata2[,2:ncol(newdata2)]
  covariates_save = covariates
  omit = NULL
  M = matrix(nrow=ncol(covariates)-1,ncol=3)
  outer_iters = ncol(covariates) - 1
  for (outer in 1:outer_iters) {
    criterion = numeric(ncol(covariates))
    for (k in ncol(covariates):1) {
      fit = lm(data.frame(y[-test_ind],covariates[-test_ind,-k]))
      criterion[k] = summary(fit)$adj.r.squared #rsquared
      #criterion[k] = BIC(fit) #BIC
    }
    remove = which(criterion==max(criterion))[1]
    omit = c(omit,colnames(covariates)[remove])
    covariates = covariates[,-remove]
    
    
    M[outer,] = c(outer,which(criterion==max(criterion))[1],
                  criterion[which(criterion==max(criterion))[1]])
  }
  
  minimizer = which(M[,3]==max(M[,3]))[length(which(M[,3]==max(M[,3])))]
  
  removed_terms = omit[1:(minimizer)]
  
  model = lm(data.frame(y[-test_ind],
                        covariates_save[-test_ind,][,!colnames(covariates_save[-test_ind,])%in%removed_terms]))
  
  mse = mean((predict(model,covariates_save[test_ind,]) - y[test_ind])^2)
  
  return(mse)
}
#############################################################################
cv = function(lambda,newdata2,test_indices,ridge) {
  y = newdata2$biddy1
  covariates = newdata2[,2:ncol(newdata2)]
  train = data.frame(y[-test_ind],covariates[-test_ind,])
  
  sum = 0
  for (k in 1:10) {
    
    test_indices_cv = seq(floor(seq(1,nrow(train),length.out=11))[k],
                          floor(seq(1,nrow(train),length.out=11))[k+1],1)
    
    test_x = as.matrix(train[test_indices_cv,2:ncol(train)])
    train_x = as.matrix(train[-test_indices_cv,2:ncol(train)])
    
    test_y = train[,1][test_indices_cv]
    train_y = train[,1][-test_indices_cv]
    
    if (ridge) {
      fit = glmnet(train_x,train_y,alpha=0,lambda=lambda)
      sum = sum + sum((predict(fit,test_x,s=lambda) - test_y)^2)
    } else {
      fit = glmnet(train_x,train_y,alpha=1,lambda=lambda)
      sum = sum + sum((predict(fit,test_x,s=lambda) - test_y)^2)
    }
  }
  return(sum)
}
#############################################################################
ridge_lasso = function(newdata2,test_ind,lambda,ridge) {
  y = newdata2$biddy1
  covariates = newdata2[,2:ncol(newdata2)]
  if (ridge) {
    alpha = 0
  } else {
    alpha = 1
  }
  
  fit = glmnet(as.matrix(covariates[-test_ind,]),y[-test_ind],alpha=alpha,lambda=lambda)
  mse = mean((predict(fit,as.matrix(covariates[test_ind,]),s=lambda) - y[test_ind])^2)
  
  return(mse)
}
#############################################################################