#set wd 
setwd("/data/gent/vo/000/gvo00056/vsc42347/InternshipLeuven/H1_fastq/9_RnBeads")

#create output directory where you want to save plots
dir.create("/data/gent/vo/000/gvo00056/vsc42347/InternshipLeuven/H1_fastq/9_RnBeads/heatmap")

#load packages
if(!require(RnBeads)){
    install.packages("RnBeads")
    library(RnBeads)
}

#BiocManager::install("RnBeads")
#library(RnBeads)

BiocManager::install("RnBeads.hg19") 
library(RnBeads.hg19) 

rnb.get.assemblies()

#specify directory with files
data.dir <- "/data/gent/vo/000/gvo00056/vsc42347/InternshipLeuven1/Methylation_data_Huiwen/File2JM/bismarkcall/" 

#specify output driectory
report.dir = "/kyukon/data/gent/vo/000/gvo00056/vsc42347/InternshipLeuven/H1_fastq/9_RnBeads/heatmap" 

#load data and annotation file located in wd
data.source = c(data.dir, "mt5.csv") 

#create rnb object
rnb.options(import.bed.style = "bismarkCov")
rnb <- rnb.run.import(data.source = data.source, data.type = "bed.dir", dir.reports = report.dir) 

#check rnb object
rnb

#filter out CpGs on sex chromosomes
rnb.set.filtered <- rnb.execute.sex.removal(rnb$rnb.set) 

#check filtered dataset
rnb.set.filtered

#impute missing values
rnb.filtered.impute <- rnb.execute.imputation(rnb.set.filtered$dataset, method = "mean.cpgs", update.ff = TRUE) 

#extract promoter methylation info
X.impute <- meth(rnb.filtered.impute, type="promoters") 

#create heatmap
clusterings.promoters <- rnb.execute.clustering(rnb.filtered.impute, region.type="promoters") 

cresult <- clusterings.promoters[[7]]@result 

attr(cresult, "class") <- "hclust" 

cresult <- as.dendrogram(cresult) 

pdf(file="promoter_heatmap.pdf") 

heatmap.2(X.impute, Rowv=TRUE, Colv=cresult, dendrogram="both", scale="none", trace="none") 

dev.off() 
