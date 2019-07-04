############
# CB: 450K #
############

rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_crossSectionalCellTypeAdjModels.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/pheno450.long.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesame450K.batchAdj.Mmatrix.Rdata")

#subset for t1d CWB
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$T1Dgroup==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$Array==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$Array)

#select only CWB times 
pheno.t1d.CWB = pheno.t1d[which(pheno.t1d$casegroup2=="CWB"),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.CWB$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d.CWB$Array)

#remove technical replicates;
pheno.t1d.CWB2 = pheno.t1d.CWB[!duplicated(pheno.t1d.CWB$key),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.CWB2$Array)]
table(colnames(M.norm.batch.2)==pheno.t1d.CWB2$Array)

#check subjects run in longitudinal 
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")

table(pheno.t1d.CWB2$IVYOmicsMethylationMaster_ID %in% long450.subjects)
table(pheno.t1d.CWB2$IVYOmicsMethylationMaster_ID %in% longEPIC.subjects)

#look into the 1 subject not in either longitudinal analysis
pheno.t1d.CWB2[-which(pheno.t1d.CWB2$IVYOmicsMethylationMaster_ID %in% long450.subjects),]

CB450.arrays = data.frame(pheno.t1d.CWB2[,c(1,2)], "yes")
colnames(CB450.arrays) = c("key", "CB.450.array", "CB.450.analysis")
save(CB450.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/CB.450K.arrays.Rdata")

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.CWB2$T1Dgroup)
t = as.factor(pheno.t1d.CWB2$t.y)
ID = as.factor(pheno.t1d.CWB2$ID)
gender=as.factor(pheno.t1d.CWB2$SEX.x)
clinage=pheno.t1d.CWB2$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender)
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
#test  = M.norm.batch.2[1:10,]
#CB.450.results  = data.frame(t(apply(test, 1, getCrossSection)))
#summary(MMtest)
####YES!!!!! It works!!!!! Do on full dataset; 

CB.450.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

CB.450.results$FDR_group = p.adjust(CB.450.results$groupT1D.control_pvalue, method="BH")
CB.450.results$FDR_gender = p.adjust(CB.450.results$genderMale_pvalue, method="BH")

CB.450.results$Mcase = CB.450.results$X.Intercept.
CB.450.results$Mcontrol = CB.450.results$X.Intercept. + CB.450.results$groupT1D.control
 
CB.450.results$Bcontrol = getBeta(CB.450.results$Mcontrol)
CB.450.results$Bcase = getBeta(CB.450.results$Mcase)

CB.450.results$dfBeta = CB.450.results$Bcase - CB.450.results$Bcontrol

save(CB.450.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/CB.450K.Rdata")

###############
# PreSV: 450K #
###############
rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#Did this for the preSV cross-section
#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_crossSectionalCellTypeAdjModels.Rdata")


load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/pheno450.long.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesame450K.batchAdj.Mmatrix.Rdata")

#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$T1Dgroup==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$Array==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$Array)

#remove technical replicates;
pheno.t1d.2 = pheno.t1d[!duplicated(pheno.t1d$key),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.2$Array)]
table(colnames(M.norm.batch.2)==pheno.t1d.2$Array)


#get only preSV times 
pheno.t1d.3 = pheno.t1d.2[which(pheno.t1d.2$casegroup2=="EV" | pheno.t1d.2$casegroup2=="PSV"),]
duplicatedIDs = names(table(pheno.t1d.3$ID))[which(table(pheno.t1d.3$ID)==2)]
p1 = pheno.t1d.3[which(pheno.t1d.3$ID %in% duplicatedIDs),]
p1 = p1[which(p1$casegroup2=="PSV"),]
pheno.t1d.presv = rbind(p1, pheno.t1d.3[-which(pheno.t1d.3$ID %in% duplicatedIDs),])


M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.presv$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d.presv$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d.presv$Array==colnames(M.norm.batch.2)[i]))
}
pheno.t1d.presv = pheno.t1d.presv[order,]
table(colnames(M.norm.batch.2)==pheno.t1d.presv$Array)

##remove any subjects run on the EPIC 
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")

pheno.t1d.presv =pheno.t1d.presv[-which(pheno.t1d.presv$IVYOmicsMethylationMaster_ID %in% longEPIC.subjects),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.presv$Array)]
table(colnames(M.norm.batch.2)==pheno.t1d.presv$Array)


preSV.450.arrays = data.frame(pheno.t1d.presv[,c(1,2)], "yes")
colnames(preSV.450.arrays) = c("key", "preSV.450.array", "preSV.450.analysis")
save(preSV.450.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/preSV.450K.arrays.Rdata")



### Statistical Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.presv$T1Dgroup)
t = as.factor(pheno.t1d.presv$t.y)
ID = as.factor(pheno.t1d.presv$ID)
gender=as.factor(pheno.t1d.presv$SEX.x)
clinage=pheno.t1d.presv$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender + clinage)
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
#test = M.norm.batch.2[1:10,]
#preSV.450.results = data.frame(t(apply(test, 1, getCrossSection)))
####YES!!!!! It works!!!!! Do on full dataset; 

preSV.450.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

preSV.450.results$FDR_group = p.adjust(preSV.450.results$groupT1D.control_pvalue, method="BH")
preSV.450.results$FDR_gender = p.adjust(preSV.450.results$genderMale_pvalue, method="BH")
preSV.450.results$FDR_clinage = p.adjust(preSV.450.results$clinage_pvalue, method="BH")

preSV.450.results$Mcase = preSV.450.results$X.Intercept.
preSV.450.results$Mcontrol = preSV.450.results$X.Intercept. + preSV.450.results$groupT1D.control

preSV.450.results$Bcontrol = getBeta(preSV.450.results$Mcontrol)
preSV.450.results$Bcase = getBeta(preSV.450.results$Mcase)

preSV.450.results$dfBeta = preSV.450.results$Bcase - preSV.450.results$Bcontrol

save(preSV.450.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/preSV.450K.Rdata")

################
# postSV: 450K #
################
rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#Did this for the CB cross-section
#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_crossSectionalCellTypeAdjModels.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/pheno450.long.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesame450K.batchAdj.Mmatrix.Rdata")

#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$T1Dgroup==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$Array==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$Array)

#remove technical replicates;
pheno.t1d.2 = pheno.t1d[!duplicated(pheno.t1d$key),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.2$Array)]
table(colnames(M.norm.batch.2)==pheno.t1d.2$Array)

#get only postSV times 
pheno.t1d.3 = pheno.t1d.2[which(pheno.t1d.2$casegroup2=="SV" | pheno.t1d.2$casegroup2=="PreT1D"),]
duplicatedIDs = names(table(pheno.t1d.3$ID))[which(table(pheno.t1d.3$ID)==2)]
p1 = pheno.t1d.3[which(pheno.t1d.3$ID %in% duplicatedIDs),]
p1 = p1[-which(p1$casegroup2=="PreT1D"),]
pheno.t1d.postsv = rbind(p1, pheno.t1d.3[-which(pheno.t1d.3$ID %in% duplicatedIDs),])


M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.postsv$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d.postsv$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d.postsv$Array==colnames(M.norm.batch.2)[i]))
}
pheno.t1d.postsv = pheno.t1d.postsv[order,]
table(colnames(M.norm.batch.2)==pheno.t1d.postsv$Array)

##remove the subject run on the EPIC;
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")
pheno.t1d.postsv =pheno.t1d.postsv[-which(pheno.t1d.postsv$IVYOmicsMethylationMaster_ID %in% longEPIC.subjects),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.postsv$Array)]

table(colnames(M.norm.batch.2)==pheno.t1d.postsv$Array)


postSV.450.arrays = data.frame(pheno.t1d.postsv[,c(1,2)], "yes")
colnames(postSV.450.arrays) = c("key", "postSV.450.array", "postSV.450.analysis")
save(postSV.450.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/postSV.450K.arrays.Rdata")

### Statistical Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.postsv$T1Dgroup)
t = as.factor(pheno.t1d.postsv$t.y)
ID = as.factor(pheno.t1d.postsv$ID)
gender=as.factor(pheno.t1d.postsv$SEX.x)
clinage=pheno.t1d.postsv$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender + clinage)
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
#test = M.norm.batch.2[1:10,]
#MMtest = data.frame(t(apply(test, 1, getCrossSection)))
####YES!!!!! It works!!!!! Do on full dataset; 

postSV.450.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

postSV.450.results$FDR_group = p.adjust(postSV.450.results$groupT1D.control_pvalue, method="BH")
postSV.450.results$FDR_gender = p.adjust(postSV.450.results$genderMale_pvalue, method="BH")
postSV.450.results$FDR_clinage = p.adjust(postSV.450.results$clinage_pvalue, method="BH")

postSV.450.results$Mcase = postSV.450.results$X.Intercept.
postSV.450.results$Mcontrol = postSV.450.results$X.Intercept. + postSV.450.results$groupT1D.control

postSV.450.results$Bcontrol = getBeta(postSV.450.results$Mcontrol)
postSV.450.results$Bcase = getBeta(postSV.450.results$Mcase)

postSV.450.results$dfBeta = postSV.450.results$Bcase - postSV.450.results$Bcontrol


save(postSV.450.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/postSV.450K.Rdata")


############
# CB: EPIC #
############

rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#saved this from the 450K cross section models ran just prior to this.
#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_cellTypeAdjModels.Rdata")

### Load 850K (individually normalized) Data ###

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesameEPIC.batchAdj.Mmatrix.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/EPICpheno.long.Rdata")

#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$group3==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$rgName)]

table(colnames(M.norm.batch.2)==pheno.t1d$rgName)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$rgName==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$rgName)

#keep only CWB times 
pheno.t1d.CB = pheno.t1d[which(pheno.t1d$casegroup2=="CWB"),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.CB$rgName)]

table(colnames(M.norm.batch.2)==pheno.t1d.CB$rgName)
pheno.t1d.CB$group3 = gsub("Case", "case", pheno.t1d.CB$group3)

#no technical replicates within the 850K to worry about;
#remove technical replicates between platforms;
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/CB.450K.arrays.Rdata")
CB450.ID = sapply(strsplit(CB450.arrays$key, split="_", fixed=TRUE), "[[", 1)

pheno.t1d.CB2 = pheno.t1d.CB[-which(pheno.t1d.CB$ID.x %in% CB450.ID),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.CB2$rgName)]

CB.EPIC.arrays = data.frame(pheno.t1d.CB2[,c(1,2)], "yes")
colnames(CB.EPIC.arrays) = c("key", "CB.EPIC.array", "CB.EPIC.analysis")
save(CB.EPIC.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/CB.EPIC.arrays.Rdata")

### Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.CB2$group3)
ID = as.factor(pheno.t1d.CB2$ID.x)
gender=as.factor(pheno.t1d.CB2$SEX.x)
clinage=pheno.t1d.CB2$clinage


getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender )
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
 
#test = M.norm.batch.2[356180:356190,]
#CB.850.results = data.frame(t(apply(test, 1, getCrossSection)))
####YES!!!!! It works!!!!! Do on full dataset; 

CB.850.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

CB.850.results$FDR_group = p.adjust(CB.850.results$groupT1D.control_pvalue, method="BH")
CB.850.results$FDR_gender = p.adjust(CB.850.results$genderMale_pvalue, method="BH")

CB.850.results$Mcase = CB.850.results$X.Intercept.
CB.850.results$Mcontrol = CB.850.results$X.Intercept. + CB.850.results$groupT1D.control
 
CB.850.results$Bcontrol = getBeta(CB.850.results$Mcontrol)
CB.850.results$Bcase = getBeta(CB.850.results$Mcase)

CB.850.results$dfBeta = CB.850.results$Bcase - CB.850.results$Bcontrol

save(CB.850.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/CB.EPIC.Rdata")

###############
# preSV: EPIC #
###############


rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#saved this from the 450K cross section models ran just prior to this.
#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_cellTypeAdjModels.Rdata")

### Load 850K (individually normalized) Data ###
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesameEPIC.batchAdj.Mmatrix.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/EPICpheno.long.Rdata")

#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$group3==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$rgName)]

table(colnames(M.norm.batch.2)==pheno.t1d$rgName)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$rgName==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$rgName)

#get only preSV times 
pheno.t1d.3 = pheno.t1d[which(pheno.t1d$casegroup2=="EV" | pheno.t1d$casegroup2=="PSV"),]
duplicatedIDs = names(table(pheno.t1d.3$ID.x))[which(table(pheno.t1d.3$ID.x)==2)]
p1 = pheno.t1d.3[which(pheno.t1d.3$ID.x %in% duplicatedIDs),]
p1 = p1[which(p1$casegroup2=="PSV"),]
pheno.t1d.presv = rbind(p1, pheno.t1d.3[-which(pheno.t1d.3$ID.x %in% duplicatedIDs),])


M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.presv$rgName)]
table(colnames(M.norm.batch.2)==pheno.t1d.presv$rgName)
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d.presv$rgName==colnames(M.norm.batch.2)[i]))
}
pheno.t1d.presv = pheno.t1d.presv[order,]
table(colnames(M.norm.batch.2)==pheno.t1d.presv$rgName)


##remove subjects from the 450K panel;
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")

length(which(pheno.t1d.presv$ID.x %in% long450.subjects))
#since this is 0, we are good to good

preSV.EPIC.arrays = data.frame(pheno.t1d.presv[,c(1,2)], "yes")
colnames(preSV.EPIC.arrays) = c("key", "preSV.EPIC.array", "preSV.EPIC.analysis")
save(preSV.EPIC.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/preSV.EPIC.arrays.Rdata")



### Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.presv$group3)
group = gsub("T1D Case", "T1D case", group)

ID = as.factor(pheno.t1d.presv$ID.x)
gender=as.factor(pheno.t1d.presv$SEX.x)
clinage=pheno.t1d.presv$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender + clinage)
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
#test = M.norm.batch.2[356180:356190,]
#preSV.850.results = data.frame(t(apply(test, 1, getCrossSection)))
####YES!!!!! It works!!!!! Do on full dataset; 

preSV.850.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

preSV.850.results$FDR_group = p.adjust(preSV.850.results$groupT1D.control_pvalue, method="BH")
preSV.850.results$FDR_clinage = p.adjust(preSV.850.results$clinage_pvalue, method="BH")
preSV.850.results$FDR_gender = p.adjust(preSV.850.results$genderMale_pvalue, method="BH")

preSV.850.results$Mcase = preSV.850.results$X.Intercept.
preSV.850.results$Mcontrol = preSV.850.results$X.Intercept. + preSV.850.results$groupT1D.control
 
preSV.850.results$Bcontrol = getBeta(preSV.850.results$Mcontrol)
preSV.850.results$Bcase = getBeta(preSV.850.results$Mcase)

preSV.850.results$dfBeta = preSV.850.results$Bcase - preSV.850.results$Bcontrol

save(preSV.850.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/preSV.EPIC.Rdata")

################
# postSV: EPIC #
################

rm(list=ls())
options(stringsAsFactors=FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(rtracklayer)

#saved this from the 450K cross section models ran just prior to this.
#session =sessionInfo()
#save(session, file="/home/vanderll/Norris/T1Dpaper/sessionInfo_cellTypeAdjModels.Rdata")

### Load 850K (individually normalized) Data ###
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesameEPIC.batchAdj.Mmatrix.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/EPICpheno.long.Rdata")

#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$group3==""),]
M.norm.batch.2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$rgName)]

table(colnames(M.norm.batch.2)==pheno.t1d$rgName)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d$rgName==colnames(M.norm.batch.2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm.batch.2)==pheno.t1d$rgName)


#get only postSV times 
pheno.t1d.3 = pheno.t1d[which(pheno.t1d$casegroup2=="SV" | pheno.t1d$casegroup2=="PreT1D"),]
duplicatedIDs = names(table(pheno.t1d.3$ID.x))[which(table(pheno.t1d.3$ID.x)==2)]
p1 = pheno.t1d.3[which(pheno.t1d.3$ID.x %in% duplicatedIDs),]
p1 = p1[-which(p1$casegroup2=="PreT1D"),]
pheno.t1d.postsv = rbind(p1, pheno.t1d.3[-which(pheno.t1d.3$ID.x %in% duplicatedIDs),])


M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.postsv$rgName)]
table(colnames(M.norm.batch.2)==pheno.t1d.postsv$rgName)
order = c()
for(i in 1:ncol(M.norm.batch.2)){
	order = c(order, which(pheno.t1d.postsv$rgName==colnames(M.norm.batch.2)[i]))
}
pheno.t1d.postsv = pheno.t1d.postsv[order,]
table(colnames(M.norm.batch.2)==pheno.t1d.postsv$rgName)


##remove subjects from the 450K panel;
load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")
pheno.t1d.postsv =pheno.t1d.postsv[-which(pheno.t1d.postsv$ID.x %in% long450.subjects),]
M.norm.batch.2 = M.norm.batch.2[,which(colnames(M.norm.batch.2) %in% pheno.t1d.postsv$rgName)]
table(colnames(M.norm.batch.2)==pheno.t1d.postsv$rgName)


postSV.EPIC.arrays = data.frame(pheno.t1d.postsv[,c(1,2)], "yes")
colnames(postSV.EPIC.arrays) = c("key", "postSV.EPIC.array", "postSV.EPIC.analysis")
save(postSV.EPIC.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/sampleInfo/postSV.EPIC.arrays.Rdata")


### Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm.batch.2)

group = as.factor(pheno.t1d.postsv$group3)
group = gsub("T1D Case", "T1D case", group)
ID = as.factor(pheno.t1d.postsv$ID.x)
gender=as.factor(pheno.t1d.postsv$SEX.x)
clinage=pheno.t1d.postsv$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}

getCrossSection = function(a){
	fit = lm(a ~ group + gender + clinage)
	
	estimate = summary(fit)$coefficients[,"Estimate"]
	stdErr = summary(fit)$coefficients[,"Std. Error"]
	names(stdErr) = paste(names(stdErr), "_stdError", sep="")
	tStat = summary(fit)$coefficients[,"t value"]
	names(tStat) = paste(names(tStat), "_tStatistic", sep="")
	pVal = summary(fit)$coefficients[,"Pr(>|t|)"]
	names(pVal) = paste(names(pVal), "_pvalue", sep="")
	
	want = c(estimate, pVal, stdErr, tStat)
	return(want)
}
 
#test = M.norm.batch.2[356180:356190,]
#postSV.850.results = data.frame(t(apply(test, 1, getCrossSection)))
####YES!!!!! It works!!!!! Do on full dataset; 

postSV.850.results = data.frame(t(apply(M.norm.batch.2, 1, getCrossSection)))

postSV.850.results$FDR_group = p.adjust(postSV.850.results$groupT1D.control_pvalue, method="BH")
postSV.850.results$FDR_clinage = p.adjust(postSV.850.results$clinage_pvalue, method="BH")
postSV.850.results$FDR_gender = p.adjust(postSV.850.results$genderMale_pvalue, method="BH")

postSV.850.results$Mcase = postSV.850.results$X.Intercept.
postSV.850.results$Mcontrol = postSV.850.results$X.Intercept. + postSV.850.results$groupT1D.control
 
postSV.850.results$Bcontrol = getBeta(postSV.850.results$Mcontrol)
postSV.850.results$Bcase = getBeta(postSV.850.results$Mcase)

postSV.850.results$dfBeta = postSV.850.results$Bcase - postSV.850.results$Bcontrol

save(postSV.850.results, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/crossSections/postSV.EPIC.Rdata")



