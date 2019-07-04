rm(list=ls())

library(minfi)
library(RColorBrewer)
library(sva)
library(plyr)

### 450K ####

#methylation data
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/450K_M.matrix_pOOBAHfiltered.Rdata")

#phenotype data
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/pheno450.long.Rdata")

#make sure the order matches up
M.want = M.qc2[,which(colnames(M.qc2) %in% pheno.forAnalysis$Array)]
table(colnames(M.want) == pheno.forAnalysis$Array)
order = c()
for(i in 1:ncol(M.want)){
	order = c(order, which(pheno.forAnalysis$Array==colnames(M.want)[i]))
}
pheno.forAnalysis = pheno.forAnalysis[order,]
table(colnames(M.want)==pheno.forAnalysis$Array)
pheno.forAnalysis$row = sapply(strsplit(pheno.forAnalysis$Sentrix_Position, split="C", fixed=TRUE), "[[", 1)

batch = paste(pheno.forAnalysis$Sample_Plate, pheno.forAnalysis$row, sep="_")
mod = model.matrix(~as.factor(group4), data=pheno.forAnalysis)

M.sesame.batch = ComBat(M.want, batch, mod)

save(M.sesame.batch, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesame450K.batchAdj.Mmatrix.Rdata")

### EPIC ####

rm(list=ls())

library(minfi)
library(RColorBrewer)
library(sva)
library(plyr)

#methylation data
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/EPIC_M.matrix_pOOBAHfiltered.Rdata")

#phenotype data
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/EPICpheno.long.Rdata")


#make sure the order matches up
M.want = M.qc2[,which(colnames(M.qc2) %in% pheno.forAnalysis$rgName)]
table(colnames(M.want) == pheno.forAnalysis$rgName)
order = c()
for(i in 1:ncol(M.want)){
	order = c(order, which(pheno.forAnalysis$rgName==colnames(M.want)[i]))
}
pheno.forAnalysis = pheno.forAnalysis[order,]
table(colnames(M.want)==pheno.forAnalysis$rgName)

pheno.forAnalysis$row = sapply(strsplit(as.character(pheno.forAnalysis$Array), split="C", fixed=TRUE), "[[", 1)

batch = paste(pheno.forAnalysis$Sample_Plate, pheno.forAnalysis$row, sep="_")
mod = model.matrix(~as.factor(group4), data=pheno.forAnalysis)

M.sesame.batch = ComBat(M.want, batch, mod)

save(M.sesame.batch, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesameEPIC.batchAdj.Mmatrix.Rdata")
