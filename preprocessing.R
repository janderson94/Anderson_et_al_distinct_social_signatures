#################################################################################################################################################
#Code for pre-processing gene expression data for downstream models. Can skip, as filtered, residual gene expression is also provided.
#################################################################################################################################################

#Read in files
meta_info<-read.table("./Data/meta_info.txt",header=T)

#Get sum of non-overlapping exon length from GTF for RPKM calculations
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("./Data/panubis1.gtf",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))

counts<-read.table("./Data/counts.txt",header=T)

genes<-as.data.frame(exonic.gene.sizes[rownames(exonic.gene.sizes)%in% rownames(counts),])
rownames(genes)<-rownames(exonic.gene.sizes)[rownames(exonic.gene.sizes)%in% rownames(counts)]
genes<-as.data.frame(genes[order(rownames(genes)),1])

tmp<-rownames(exonic.gene.sizes)[rownames(exonic.gene.sizes)%in% rownames(counts)]
rownames(genes)<-tmp[order(tmp)]

length(which(rownames(genes)==rownames(counts)))
colnames(genes)<-"length"
rm(tmp,txdb,exonic.gene.sizes,exons.list.per.gene)

########################################################################
#RPKM normalization and filtering
########################################################################
total<-colSums(counts)

#Median rpkm>1 within condition
rpkm<-10^9*counts/genes$length
rpkm2<-rpkm

for (i in 1:length(colnames(rpkm))){
  rpkm2[,i]<-rpkm[,i]/total[i]
}

mean_null<-apply(rpkm2[,meta_info$treatment=="NULL"],1,median)
mean_lps<-apply(rpkm2[,meta_info$treatment=="LPS"],1,median)

# IL6 demonstrates succesfull LPS stimulation and differences in estimated expression in baseline and LPS-stimulated samples
mean_null[names(mean_lps)=="IL6"]
mean_lps[names(mean_lps)=="IL6"]

rpkm3<-subset(rpkm2,mean_null>2 | mean_lps>2)

#This leaves us with 10,281 genes

########################################################################
#Regress out known batch effects and first 3 PCS from flow cytometry
########################################################################
design <- as.matrix(model.matrix(~as.factor(meta_info$year)+meta_info$Flow_PC1+meta_info$Flow_PC2+meta_info$Flow_PC3))
counts2<-counts[rownames(counts)%in% rownames(rpkm3),]

library(edgeR)
dge <- DGEList(counts=counts2)
dge <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge,design,plot=FALSE)
fit <-lmFit(v,design)
resids<-residuals.MArrayLM(fit,v)

#write.table(resids,"Data/residual_gene_expression.txt",quote=F,row.names=T,col.names=T,sep="\t")

rm(rpkm,rpkm2,rpkm3,counts,counts2,genes,fit,dge,design,v,i,f,mean_lps,mean_null,total)

#Residual gene expression can now be modeled in the section "linear_models"


