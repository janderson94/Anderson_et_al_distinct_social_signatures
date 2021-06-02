###########################
#Predicting dominance rank 
###########################

#######################
#Within sex predictions
#######################

#######################
#Females
#######################

###########
#NULL
###########
library(glmnet)
library(parallel)

#Custom code to run parallel elastic net predictions acros a range of alphas
parallel_predict<-function(age_vector, counts_file, sample_type, no_cores = NULL, alphas){
  if (!is.null(no_cores)) {
    no_cores <- no_cores
  }
  else {
    no_cores <- detectCores() - 1
  }
  N <- length(age_vector)
  cl <- makeCluster(no_cores)
  clusterExport(cl = cl, varlist = ls(), envir = environment())
  clusterEvalQ(cl, library(glmnet))
  predicted <- matrix(NA, nrow = N, ncol = length(alphas))
  colnames(predicted) <- alphas
  par_out <- parLapply(cl, X = alphas, function(z) {
    for (i in 1:N) {
      if (sample_type == "GE") {
        norm_counts <- as.matrix(apply(counts_file, 2, 
                                       function(y) {
                                         return(qqnorm(log2(y/sum(y) + 1), plot = F)$x)
                                       }))
      }
      else if (sample_type == "epi") {
        norm_counts <- as.matrix(apply(counts_file, 2, 
                                       function(y) {
                                         return(qqnorm(y, plot = F)$x)
                                       }))
      }
      else {
        return("invalid sample_type")
      }
      norm_train <- norm_counts[, -i]
      norm_test <- norm_counts[, i]
      if (sample_type == "epi") {
        trainreads_norm <- as.matrix(apply(norm_train, 
                                           1, function(y) {
                                             return(qqnorm(y, plot = F)$x)
                                           }))
      }
      else if (sample_type == "GE") {
        trainreads_norm <- as.matrix(apply(norm_train, 
                                           1, function(x) {
                                             return(get_z(x))
                                           }))
      }
      else {
        return("invalid sample_type")
      }
      trainage <- age_vector[-i]
      testage <- age_vector[i]
      if (sample_type == "epi") {
        testreads_normalized <- norm_test
        for (d in 1:dim(norm_counts)[1]) {
          a <- ecdf(norm_train[d, ])
          probs <- a(norm_test[d])
          probs[probs == 1] <- 0.99
          probs[probs == 0] <- 0.01
          testreads_normalized[d] <- qnorm(probs)
        }
      }
      else if (sample_type == "GE") {
        testreads_normalized <- norm_test
        for (d in 1:dim(norm_counts)[1]) {
          mean_d <- mean(norm_train[d, ])
          sd_d <- sd(norm_train[d, ])
          testreads_normalized[d] <- (norm_test[d] - 
                                        mean_d)/sd_d
        }
      }
      else {
        return("invalid sample_type")
      }
      model <- cv.glmnet(trainreads_norm, trainage, nfolds = (N - 1), alpha = z, standardize = F)
      predicted[i, colnames(predicted) == z] <- predict(model, 
                                                        newx = t(testreads_normalized), s = "lambda.min")
    }
    return(predicted)
  })
  stopCluster(cl)
  return(par_out)
}


predicted_female_SI<-matrix(NA,nrow=length(which(meta_info_female$treatment=="NULL")),ncol=10)

for(f in 1:10){
predicted_female_SI[,f]<-parallel_predict(age_vector = meta_info_female$rank[meta_info_female$treatment=="NULL"],
                                       counts_file = resids_female[,meta_info_female$treatment=="NULL" ],sample_type = "epi",alphas = seq(.1,1,by = .1)[f]
                                        
                                        
                                        )[[1]]
print(f)
}
head(predicted_female_SI)

rsq<-1:10
for(f in 1:10){
  rsq[f]<-summary(lm(predicted_female_SI[which(meta_info_female$sname[meta_info_female$treatment=="NULL"]!="AMB_2"),f]~meta_info_female$rank[meta_info_female$treatment=="NULL" & meta_info_female$sname!="AMB_2"]))$r.squared
}

rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_female_SI[,9]~meta_info_female$rank[meta_info_female$treatment=="NULL"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_female_SI[,9]~meta_info_female$rank[meta_info_female$treatment=="NULL"]),lty=2)
summary(lm(predicted_female_SI[,9]~meta_info_female$rank[meta_info_female$treatment=="NULL"]))

#Without AMB_2
which(meta_info_female$sname[meta_info_female$treatment=="NULL"]=="AMB_2")
summary(lm(predicted_female_SI[-22,9]~meta_info_female$rank[meta_info_female$treatment=="NULL"][-22]))

###########
#LPS
###########
predicted_female<-matrix(NA,nrow=length(which(meta_info_female$treatment=="LPS")),ncol=10)
for(f in 1:10){
  predicted_female[,f]<-parallel_predict(age_vector = meta_info_female$rank[meta_info_female$treatment=="LPS"],
                                         counts_file = resids_female[,meta_info_female$treatment=="LPS" ],sample_type = "epi",alphas = seq(.1,1,by = .1)[f]
                                         
                                         
  )[[1]]
  print(f)
}
head(predicted_female)
rsq<-1:10
for(f in 1:10){
  rsq[f]<-summary(lm(predicted_female[which(meta_info_female$sname[meta_info_female$treatment=="LPS"]!="AMB_2"),f]~meta_info_female$rank[meta_info_female$treatment=="LPS" & meta_info_female$sname!="AMB_2"]))$r.squared
}

rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_female[,10]~meta_info_female$rank[meta_info_female$treatment=="LPS"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_female[,10]~meta_info_female$rank[meta_info_female$treatment=="LPS"]),lty=2)
summary(lm(predicted_female[,10]~meta_info_female$rank[meta_info_female$treatment=="LPS"]))





#######################
#Males
#######################

###########
#NULL
###########
predicted_male<-matrix(NA,nrow=length(which(meta_info_male$treatment=="NULL")),ncol=10)

for(f in 1:10){
  predicted_male[,f]<-parallel_predict(age_vector = meta_info_male$rank[meta_info_male$treatment=="NULL"],
                                         counts_file = resids_male[,meta_info_male$treatment=="NULL" ],sample_type = "epi",alphas = seq(.1,1,by = .1)[f]
                                         
                                         
  )[[1]]
  print(f)
}
head(predicted_male)
rsq<-1:10
for(f in 1:10){
  rsq[f]<-summary(lm(predicted_male[,f]~meta_info_male$rank[meta_info_male$treatment=="NULL"]))$r.squared
}
rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_male[,3]~meta_info_male$rank[meta_info_male$treatment=="NULL"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_male[,3]~meta_info_male$rank[meta_info_male$treatment=="NULL"]),lty=2)
summary(lm(predicted_male[,3]~meta_info_male$rank[meta_info_male$treatment=="NULL"]))

###########
#LPS
###########
predicted_male_SI<-matrix(NA,nrow=length(which(meta_info_male$treatment=="LPS")),ncol=10)

for(f in 1:10){
  predicted_male_SI[,f]<-parallel_predict(age_vector = meta_info_male$rank[meta_info_male$treatment=="LPS"],
                                       counts_file = resids_male[,meta_info_male$treatment=="LPS" ],sample_type = "epi",alphas = seq(.1,1,by = .1)[f]
                                       
                                       
  )[[1]]
  print(f)
}

head(predicted_male_SI)
rsq<-1:10
for(f in 1:10){
  rsq[f]<-summary(lm(predicted_male_SI[,f]~meta_info_male$rank[meta_info_male$treatment=="LPS"]))$r.squared
}
rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_male_SI[,1]~meta_info_male$rank[meta_info_male$treatment=="LPS"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_male_SI[,1]~meta_info_male$rank[meta_info_male$treatment=="LPS"]),lty=2)
summary(lm(predicted_male_SI[,1]~meta_info_male$rank[meta_info_male$treatment=="LPS"]))


#################
#Across sex predictions
#################
predicted_female_across<-matrix(NA,nrow=length(which(meta_info_female$treatment=="NULL")),ncol=10)

for(f in seq(10)){
  model<-cv.glmnet(x = t(resids_male[,meta_info_male$treatment=="NULL"]),
                   y=meta_info_male$rank[meta_info_male$treatment=="NULL"],alpha=seq(.1,1,by=.1)[f],nfolds = length(meta_info_male$rank[meta_info_male$treatment=="NULL"]))
  
  predicted_female_across[,f]<-predict(object = model,newx = t(resids_female[,meta_info_female$treatment=="NULL"]),s = "lambda.min")
  print(f)
}
rsq<-1:10

for(f in 1:10){
  rsq[f]<-summary(lm(predicted_female_across[,f]~meta_info_female$rank[meta_info_female$treatment=="NULL"]))$r.squared
}

rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_female_across[,1]~meta_info_female$rank[meta_info_female$treatment=="NULL"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_female_across[,1]~meta_info_female$rank[meta_info_female$treatment=="NULL"]),lty=2)
summary(lm(predicted_female_across[,1]~meta_info_female$rank[meta_info_female$treatment=="NULL"]))

#exclude AMB_2
which(meta_info_female$sname[meta_info_female$treatment=="NULL"]=="AMB_2")
plot(predicted_female_across[-22,1]~meta_info_female$rank[meta_info_female$treatment=="NULL" ][-22],xlab="Rank",ylab="Predicted rank",pch=20,main="Train on males, predict females")
summary(lm(predicted_female_across[-22,1]~meta_info_female$rank[meta_info_female$treatment=="NULL" ][-22]))
abline(lm(predicted_female_across[-22,1]~meta_info_female$rank[meta_info_female$treatment=="NULL"][-22]),lty=2)

#Train on females, predict males
#Train on everyone but AMB_2
predicted_male_across<-matrix(NA,nrow=length(which(meta_info_male$treatment=="NULL")),ncol=10)

for(f in seq(10)){
  model<-cv.glmnet(x = t(resids_female[,meta_info_female$treatment=="NULL" & meta_info_female$sname!="AMB_2"]),
                   y=meta_info_female$rank[meta_info_female$treatment=="NULL"& meta_info_female$sname!="AMB_2"],alpha=seq(.1,1,by=.1)[f],
                   nfolds=length(meta_info_female$rank[meta_info_female$treatment=="NULL"& meta_info_female$sname!="AMB_2"]))
  
  predicted_male_across[,f]<-predict(object = model,newx = t(resids_male[,meta_info_male$treatment=="NULL"]),s = "lambda.min")
  print(f)
}
rsq<-1:10

for(f in 1:10){
  rsq[f]<-summary(lm(predicted_male_across[,f]~meta_info_male$rank[meta_info_male$treatment=="NULL"]))$r.squared
}

rsq

par(mfrow=c(1,1),pty="s")
plot(predicted_male_across[,5]~meta_info_male$rank[meta_info_male$treatment=="NULL"],xlab="Rank",ylab="Predicted rank",pch=20)
abline(lm(predicted_male_across[,5]~meta_info_male$rank[meta_info_male$treatment=="NULL"]),lty=2)
summary(lm(predicted_male_across[,5]~meta_info_male$rank[meta_info_male$treatment=="NULL"]))



