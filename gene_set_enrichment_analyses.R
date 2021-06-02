################################################################################################
############################# Gene Set Enrichment Analyses #####################################
################################################################################################
library(doParallel)
library(parallel)
library(qusage)

########################################################################
#Source GSEA code from the BROAD institute. 
#https://github.com/GSEA-MSigDB/GSEA_R
#I do not claim to have made or maintain the "GSEA.EnrichmentScore code. 
########################################################################

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
    
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    norm.tag <- 1/sum.correl.tag
    norm.no.tag <- 1/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
      # ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
    } else {
      # ES <- min.ES
      ES <- signif(min.ES, digits = 5)
      arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
  }

############################################################
#Read in the MSIG database Hallmark pathways
#List of 50 different hallmark pathways with associated genes
############################################################
hallmark<-read.gmt("Data/h.all.v7.2.symbols.gmt")

########################################################
#Run males rank, female rank, and female DSI all at once
########################################################

###############################
#Effects within the Null samples
###############################
hall_enrich<-matrix(NA,nrow=50,ncol=6)
hall_enrich<-as.data.frame(hall_enrich)
colnames(hall_enrich)<-c("male_rank_ES","male_rank_p",
                         "female_rank_ES","female_rank_p",
                         "female_DSI_ES","female_DSI_p")

for(f in 1:50){
  #Male rank
  #Retain genes that converged in our model (indicated by genes with reasonable SE(betas))
  keep<-which(apply(EMMA_RNA_nested_males[,7:12],1,function(x){length(which(x>0))})==6)
  o<-EMMA_RNA_nested_males[keep,5]/sqrt(EMMA_RNA_nested_males[keep,11])
  names(o)<-rownames(EMMA_RNA_nested_males)[keep]
  o<-o[order(o)]
  hall_enrich[f,1]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  temp<-1:10000
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
    })
  stopCluster(clus)
  
  #The pvalue is the % of abs(observed ES)< abs(permuted ES)  
  hall_enrich[f,2]<-length(which(abs(hall_enrich[f,1]) < abs(temp2)))/10000
  
  #Female rank
  o<-EMMA_RNA_nested_females_no_acid[keep,5]/sqrt(EMMA_RNA_nested_females_no_acid[keep,13])
  names(o)<-rownames(EMMA_RNA_nested_females_no_acid)[keep]
  o<-o[order(o)]
  hall_enrich[f,3]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(observed ES)< abs(permuted ES)  
  hall_enrich[f,4]<-length(which(abs(hall_enrich[f,3]) < abs(temp2)))/10000
  
  #Female DSI
  o<-EMMA_RNA_nested_females_no_acid[keep,7]/sqrt(EMMA_RNA_nested_females_no_acid[keep,15])
  names(o)<-rownames(EMMA_RNA_nested_females_no_acid)[keep]
  o<-o[order(o)]
  hall_enrich[f,5]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich[f,6]<-length(which(abs(hall_enrich[f,5]) < abs(temp2)))/10000
  
  
  print(f)
}
rownames(hall_enrich)<-names(hallmark)

head(hall_enrich)

#write.table(hall_enrich,"~/Desktop/hallmark_enrichment_effects_null.txt",quote=F,row.names=T,col.names=T,sep="\t")

###############################
#Effects within the LPS samples
###############################
hall_enrich_lps<-matrix(NA,nrow=50,ncol=6)
hall_enrich_lps<-as.data.frame(hall_enrich_lps)
colnames(hall_enrich_lps)<-c("male_rank_ES","male_rank_p",
                         "female_rank_ES","female_rank_p",
                         "female_DSI_ES","female_DSI_p")

for(f in 1:50){
  #Male rank
  #Retain genes that converged in our model (indicated by genes with reasonable SE(betas))
  keep<-which(apply(EMMA_RNA_nested_males[,7:12],1,function(x){length(which(x>0))})==6)
  o<-EMMA_RNA_nested_males[keep,6]/sqrt(EMMA_RNA_nested_males[keep,12])
  names(o)<-rownames(EMMA_RNA_nested_males)[keep]
  o<-o[order(o)]
  hall_enrich_lps[f,1]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  temp<-1:10000
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich_lps[f,2]<-length(which(abs(hall_enrich_lps[f,1]) < abs(temp2)))/10000
  
  #Female rank
  #Retain genes that converged in our model (indicated by genes with reasonable SE(betas))
  keep<-which(apply(EMMA_RNA_nested_females_no_acid[,9:16],1,function(x){length(which(x>0))})==8)
  o<-EMMA_RNA_nested_females_no_acid[keep,6]/sqrt(EMMA_RNA_nested_females_no_acid[keep,14])
  names(o)<-rownames(EMMA_RNA_nested_females_no_acid)[keep]
  o<-o[order(o)]
  hall_enrich_lps[f,3]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich_lps[f,4]<-length(which(abs(hall_enrich_lps[f,3]) < abs(temp2)))/10000
  
  #Female DSI
  o<-EMMA_RNA_nested_females_no_acid[keep,8]/sqrt(EMMA_RNA_nested_females_no_acid[keep,16])
  names(o)<-rownames(EMMA_RNA_nested_females_no_acid)[keep]
  o<-o[order(o)]
  hall_enrich_lps[f,5]<-GSEA.EnrichmentScore(gene.list = names(o),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = o,weighted.score.type = 1)$ES
  
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich_lps[f,6]<-length(which(abs(hall_enrich_lps[f,5]) < abs(temp2)))/10000
  
  
  print(f)
}
rownames(hall_enrich_lps)<-names(hallmark)

head(hall_enrich_lps)

#write.table(hall_enrich_lps,"~/Desktop/hallmark_enrichment_effects_lps.txt",quote=F,row.names=T,col.names=T,sep="\t")


###################################
#Effects from the interaction model
###################################
hall_enrich_interaction<-matrix(NA,nrow=50,ncol=4)
hall_enrich_interaction<-as.data.frame(hall_enrich_interaction)
colnames(hall_enrich_interaction)<-c("female_rank_ES","female_rank_p",
                             "female_DSI_ES","female_DSI_p")


temp<-1:10000
for(f in 1:50){
  #Female rank
  #Retain genes that converged in our model (indicated by genes with reasonable SE(betas))
  keep<-which(apply(EMMA_RNA_nested_females_fc_model[,5:8],1,function(x){length(which(x>0))})==4)
  o<-EMMA_RNA_nested_females_fc_model[keep,3]/sqrt(EMMA_RNA_nested_females_fc_model[keep,7])
  names(o)<-rownames(EMMA_RNA_nested_females_fc_model)[keep]
  o<-o[order(o)]
  hall_enrich_interaction[f,1]<-GSEA.EnrichmentScore(gene.list = names(o),
                                             gene.set = unlist(hallmark[[f]]),
                                             correl.vector = o,weighted.score.type = 1)$ES
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich_interaction[f,2]<-length(which(abs(hall_enrich_interaction[f,1]) < abs(temp2)))/10000
  
  #Female DSI
  o<-EMMA_RNA_nested_females_fc_model[keep,4]/sqrt(EMMA_RNA_nested_females_fc_model[keep,8])
  names(o)<-rownames(EMMA_RNA_nested_females_fc_model)[keep]
  o<-o[order(o)]
  hall_enrich_interaction[f,3]<-GSEA.EnrichmentScore(gene.list = names(o),
                                             gene.set = unlist(hallmark[[f]]),
                                             correl.vector = o,weighted.score.type = 1)$ES
  
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark","keep",'o',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:10000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(keep),size= length(keep),replace=FALSE)
    o_perm<-o[rand]
    names(o_perm)<-names(o)
    
    o_perm<-o_perm[order(o_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(o_perm), gene.set = unlist(hallmark[[f]]),correl.vector = o_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(permuted ES) > abs(observed ES)
  hall_enrich_interaction[f,4]<-length(which(abs(hall_enrich_interaction[f,3]) < abs(temp2)))/10000
  
  
  print(f)
}
rownames(hall_enrich_interaction)<-names(hallmark)

head(hall_enrich_interaction)

#write.table(hall_enrich_interaction,"~/Desktop/SI_tables/hallmark_enrichment_effects_interaction.txt",quote=F,row.names=T,col.names=T,sep="\t")

