###################################################
#Linear mixed effects models for males and females
###################################################
meta_info<-read.table("./Data/meta_info.txt",header=T)
meta_info$treatment<-relevel(meta_info$treatment,ref="NULL")
meta_info$sname<-gsub(gsub(meta_info$sra_sname,pattern="_NULL",replacement = ""),pattern="_LPS",replacement = "")

#Reading in our K matrix of genetic correlations between samples
#Repeat samples from the same individuals will have genetic correlations of 1
rel<-read.table("./Data/K_matrix.txt",header=T)
length(which(colnames(rel)==meta_info$sra_sname))

############################################
############  Male model ##################
############################################
meta_info_male<-meta_info[which(meta_info$sex=="M"),]

#Rescale treatment values for balanced nesting
meta_info_male$treatment_resc<-as.numeric(meta_info_male$treatment)-1.5

#Male model: Gene expression ~ Intercept + LPS effect + Age (in Null samples) + Age (in LPS samples) + Dominance rank (in Null samples) + Dominance rank (in LPS samples)
design<-cbind.data.frame(rep(1,length(meta_info_male[,1])),meta_info_male$treatment_resc,scale(meta_info_male$age),scale(meta_info_male$age),scale(meta_info_male$rank),scale(meta_info_male$rank))
colnames(design)<-c("intercept","LPS","age_null","age_lps","rank_null","rank_lps")

#Nesting age and rank within treatment
design[meta_info_male$treatment=="LPS",3]<-0
design[meta_info_male$treatment=="NULL",4]<-0
design[meta_info_male$treatment=="LPS",5]<-0
design[meta_info_male$treatment=="NULL",6]<-0

z<-which(complete.cases(design))

rel<-as.matrix(rel)
rel_male<-rel[which(meta_info$sex=="M"),which(meta_info$sex=="M")]
resids_male<-resids[,which(meta_info$sex=="M")]

library(doParallel)
#Can change depending on the number of cores available
clus <- makeCluster(7)
registerDoParallel(cores=7) 
design<-as.matrix(design)

clusterExport(clus,varlist=c("resids_male","design",'z',"rel_male"),envir=environment())
EMMA_RNA_nested=t(parApply(clus,resids_male,1,function(y){
  library(EMMREML)
  emma=emmreml(y=y,X=design,Z=diag(length(z)),K=rel_male,varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
stopCluster(clus)

dim(EMMA_RNA_nested)
par(mfrow=c(2,3),cex=1.25)
for(f in 2:6){
  hist(EMMA_RNA_nested[,f+12],breaks=100,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
}

par(mfrow=c(1,3),cex=1.25)
for(f in c(2,5,6)){
  hist(EMMA_RNA_nested[,f+12],breaks=100,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
}

EMMA_RNA_nested_males<-as.data.frame(EMMA_RNA_nested)

library(qvalue)
length(which(qvalue(EMMA_RNA_nested_males[,17])$qvalues<.1))
length(which(qvalue(EMMA_RNA_nested_males[,18])$qvalues<.1))


############################################
############ Female model ##################
############################################

#Female model: Gene expression ~ Intercept + LPS effect + Age (in Null samples) + Age (in LPS samples) + Dominance rank (in Null samples) + 
#Dominance rank (in LPS samples) + Social Bond Strength (in NULL samples) + Social Bond Strength (in LPS samples)
meta_info_female<-meta_info[meta_info$sex=="F",]
meta_info_female$treatment_resc<-as.numeric(meta_info_female$treatment)-1.5

design<-cbind.data.frame(rep(1,length(meta_info_female[,1])),meta_info_female$treatment_resc,scale(meta_info_female$age),scale(meta_info_female$age),scale(meta_info_female$rank),scale(meta_info_female$rank),scale(meta_info_female$DSI_F),scale(meta_info_female$DSI_F))

#Again nesting each of our effects within treatment
design[meta_info_female$treatment=="LPS",3]<-0
design[meta_info_female$treatment=="NULL",4]<-0
design[meta_info_female$treatment=="LPS",5]<-0
design[meta_info_female$treatment=="NULL",6]<-0
design[meta_info_female$treatment=="LPS",7]<-0
design[meta_info_female$treatment=="NULL",8]<-0

colnames(design)<-c("intercept","LPS","age_null","age_lps","rank_null","rank_lps","DSI_null","DSI_LPS")

clus <- makeCluster(7)
registerDoParallel(cores=7) 
design<-as.matrix(design)
resids_female<-resids[,which(meta_info$sex=="F")]
rel_female<-rel[which(meta_info$sex=="F"),which(meta_info$sex=="F")]

z<-which(complete.cases(design))
clusterExport(clus,varlist=c("resids_female","design",'z',"rel_female"),envir=environment())

EMMA_RNA_nested=t(parApply(clus,resids_female,1,function(y){
  library(EMMREML)
  emma=emmreml(y=y,X=design[z,],Z=diag(length(z)),K=rel_female[z,z],varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))

dim(EMMA_RNA_nested)
par(mfrow=c(2,4),cex=1.25,pty="s")
for(f in 1:8){
  hist(EMMA_RNA_nested[,f+16],breaks=250,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
  
}

par(mfrow=c(1,4),cex=1.25,pty="s")
for(f in 5:8){
  hist(EMMA_RNA_nested[,f+16],breaks=250,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
  
}

EMMA_RNA_nested_females<-as.data.frame(EMMA_RNA_nested)

#Excluding AMB_2
#Excluding missing data and SADIE only first
z<-which(complete.cases(design) & meta_info_female$sname!="AMB_2")
clusterExport(clus,varlist=c("resids_female","design",'z',"rel_female"),envir=environment())
EMMA_RNA_nested=t(parApply(clus,resids_female[,z],1,function(y){
  library(EMMREML)
  emma=emmreml(y=y,X=design[z,],Z=diag(length(z)),K=rel_female[z,z],varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))

dim(EMMA_RNA_nested)
par(mfrow=c(2,4),cex=1.25,pty="s")
for(f in 1:8){
  hist(EMMA_RNA_nested[,f+16],breaks=250,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
  
}

par(mfrow=c(1,4),cex=1.25,pty="s")
for(f in 5:8){
  hist(EMMA_RNA_nested[,f+16],breaks=250,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
  
}

EMMA_RNA_nested_females_no_AMB_2<-as.data.frame(EMMA_RNA_nested)

length(which(qvalue(EMMA_RNA_nested_females_no_AMB_2[,21])$qvalues<.1))
length(which(qvalue(EMMA_RNA_nested_females_no_AMB_2[,22])$qvalues<.1))

length(which(qvalue(EMMA_RNA_nested_females_no_AMB_2[,23])$qvalues<.1))
length(which(qvalue(EMMA_RNA_nested_females_no_AMB_2[,24])$qvalues<.1))

#
keep<-which(apply(EMMA_RNA_nested_females_no_AMB_2[,9:16],1,function(x){length(which(x>0 & x<1))})==8)

cor.test(EMMA_RNA_nested_females_no_AMB_2[keep,7]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,15]),
         EMMA_RNA_nested_females_no_AMB_2[keep,8]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,16]))

cor.test(EMMA_RNA_nested_females_no_AMB_2[keep,7]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,15]),
         EMMA_RNA_nested_females_no_AMB_2[keep,5]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,13]))

cor.test(EMMA_RNA_nested_females_no_AMB_2[keep,8]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,16]),
         EMMA_RNA_nested_females_no_AMB_2[keep,6]/sqrt(EMMA_RNA_nested_females_no_AMB_2[keep,14]))


############################################
######## Female fold-change model ##########
############################################
tmp<-cbind.data.frame(meta_info_female$sname)
colnames(tmp)<-c("sname")

tmp$count<-NA

for(f in 1:90){
  tmp$count[f]<-length(which(tmp$sname[f] == meta_info_female$sname))
}

tmp<-unique(tmp[tmp$count==2,])

#44 females with paired measures
fc<-matrix(NA,nrow=dim(resids_female)[1],ncol=44)
colnames(fc)<-tmp$sname
fc<-as.data.frame(fc)
rownames(fc)<-rownames(resids)

for(f in 1:44){
  fc[,f]<-
    resids_female[rownames(resids_female)%in% rownames(fc),meta_info_female$sname==tmp$sname[f] &
                    meta_info_female$treatment=="LPS"]- resids_female[rownames(resids_female)%in% rownames(fc),meta_info_female$sname==tmp$sname[f] &
                                                                        meta_info_female$treatment=="NULL"]
  
  
}

t1<-unique(cbind.data.frame(meta_info_female$sname,meta_info_female$rank,meta_info_female$DSI_F,meta_info_female$age))
colnames(t1)<-c("sname","rank","dsi","age")
tmp$order<-1:44

tmp2<-merge(tmp,t1,by=c("sname"))
tmp2<-tmp2[order(tmp2$order),]

rm(t1,tmp)

design<-cbind.data.frame(rep(1,length(tmp2$age)),tmp2$age,tmp2$rank,tmp2$dsi)
colnames(design)<-c("intercept","age","rank","DSI")

library(doParallel)
clus <- makeCluster(7)
registerDoParallel(cores=7) 
design<-as.matrix(design)

z<-which(meta_info_female$sname %in% tmp2$sname & meta_info_female$treatment=="NULL" & meta_info_female$sname!="AMB_2")
z2<-which(tmp2$sname!="AMB_2")

clusterExport(clus,varlist=c("resids_female","design","z2",'z',"rel_female"),envir=environment())
EMMA_RNA_nested=t(parApply(clus,fc[,z2],1,function(y){
  library(EMMREML)
  emma=emmreml(y=y,X=design[z2,],Z=diag(length(z2)),K=rel_female[z,z],varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))


EMMA_RNA_nested_females_fc_model<-EMMA_RNA_nested

dim(EMMA_RNA_nested)
par(mfrow=c(2,2),cex=1.25,pty="s")
for(f in 1:4){
  hist(EMMA_RNA_nested[,f+8],breaks=250,col="Steel Blue",xlab="Pvalue",main=colnames(design)[f])
  
}

rm(tmp2,f,z,z2,keep,fc,clus,design)

