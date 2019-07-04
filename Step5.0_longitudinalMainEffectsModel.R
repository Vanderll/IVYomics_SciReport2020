#########################################
# 450K: Longitudinal Main Effects Model #
#########################################

rm(list=ls())
options(stringsAsFactors=FALSE)

library(plyr)
#library(bumphunter)
#library(minfi)
#library(ChAMP)
#library(waterRmelon)
library(doParallel)
#library(DMRcate)
#library(sva)
#library(rtracklayer)
library(nlme)

session =sessionInfo()
#save(session, file="/home/linuxshared/fromHome.v2/vanderll/Norris/T1Dpaper/sessionInfo_cellTypeAdjModels.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/pheno450.long.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesame450K.batchAdj.Mmatrix.Rdata")


#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$T1Dgroup==""),]
M.norm2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$Array)]

table(colnames(M.norm2)==pheno.t1d$Array)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm2)){
	order = c(order, which(pheno.t1d$Array==colnames(M.norm2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm2)==pheno.t1d$Array)

#remove any CWB times 
pheno.t1d.noCWB = pheno.t1d[-which(pheno.t1d$casegroup2=="CWB"),]
M.norm2 = M.norm2[,which(colnames(M.norm2) %in% pheno.t1d.noCWB$Array)]

table(colnames(M.norm2)==pheno.t1d.noCWB$Array)

#remove technical replicates;
pheno.t1d.noCWB2 = pheno.t1d.noCWB[!duplicated(pheno.t1d.noCWB$key),]
M.norm2 = M.norm2[,which(colnames(M.norm2) %in% pheno.t1d.noCWB2$Array)]
table(colnames(M.norm2)==pheno.t1d.noCWB2$Array)

#remove subject 00174-0 because more time points on EPIC
pheno.t1d.final = pheno.t1d.noCWB2[-which(pheno.t1d.noCWB2$IVYOmicsMethylationMaster_ID=="00174-0"),]

M.norm2 = M.norm2[,which(colnames(M.norm2) %in% pheno.t1d.final$Array)]
table(colnames(M.norm2)==pheno.t1d.final$Array)

if(FALSE){
long450.subjects = unique(pheno.t1d.final$IVYOmicsMethylationMaster_ID)
long450.arrays = data.frame(pheno.t1d.final[,c(1, 15,24,2)], "yes")
colnames(long450.arrays) = c("key", "ID", "DOV", "Array", "T1Dlong450Analysis")

save(long450.subjects, long450.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")
}

### Longitudinal Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm2)

group = as.factor(pheno.t1d.final$T1Dgroup)
t = as.factor(pheno.t1d.final$t.y)
ID = as.factor(pheno.t1d.final$ID)
gender=as.factor(pheno.t1d.final$SEX.x)
clinage=pheno.t1d.final$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}
#a = M.want[1,]

getLongMM = function(a){
	fit.mv2 <- try(lme (a ~  clinage+group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit))
	want = rep(NA, length=36)
	if("try-error" %in% class(fit.mv2)){want = want}else{
		
	#get ANOVA values;
	p2<- anova(fit.mv2)
	numDF = p2$numDF 
	names(numDF) = paste(rownames(p2), "numDF", sep="_")
	denDF = p2$denDF
	names(denDF) = paste(rownames(p2), "denDF", sep="_")
	Fvalue = p2$`F-value`	
	names(Fvalue) = paste(rownames(p2), "Fvalue", sep="_")
	pvalue = p2$`p-value`	
	names(pvalue) = paste(rownames(p2), "pvalue", sep="_")

	#get t-test values;
	t_pvalue = summary(fit.mv2)$tTable[,"p-value"]
	names(t_pvalue) = paste(rownames(summary(fit.mv2)$tTable), "t_pvalue", sep="_")
	t_stat = summary(fit.mv2)$tTable[,"t-value"]
	names(t_stat) = paste(rownames(summary(fit.mv2)$tTable), "t_stat", sep="_")
	t_stdError = summary(fit.mv2)$tTable[,"Std.Error"]
	names(t_stdError) = paste(rownames(summary(fit.mv2)$tTable), "t_stdError", sep="_")
	t_df = summary(fit.mv2)$tTable[,"DF"]
	names(t_df) = paste(rownames(summary(fit.mv2)$tTable), "t_df", sep="_")

	want = c(summary(fit.mv2)$coefficients$fixed, pvalue, numDF, denDF, Fvalue, t_pvalue, t_stat, t_stdError, t_df)
	
	}
	
	names(want) = c("Intercept","clinage","groupT1D control","genderMale","Intercept_pvalue","clinage_pvalue","group_pvalue","gender_pvalue","Intercept_numDF","clinage_numDF","group_numDF","gender_numDF", "Intercept_denDF","clinage_denDF", "group_denDF","gender_denDF", "Intercept_Fvalue","clinage_Fvalue", "group_Fvalue","gender_Fvalue", "Intercept_t_pvalue","clinage_t_pvalue", "groupT1D control_t_pvalue","genderMale_t_pvalue", "Intercept_t_stat", "clinage_t_stat", "groupT1D control_t_stat","genderMale_t_stat", "Intercept_t_stdError","clinage_t_stdError", "groupT1D control_t_stdError", "genderMale_t_stdError", "Intercept_t_df" ,"clinage_t_df", "groupT1D control_t_df","genderMale_t_df")

	return(want)
}
 
 #test = M.norm2[356180:356190,]
#MMtest = data.frame(t(apply(test, 1, getLongMM)))
####YES!!!!! It works!!!!! Do on full dataset; 

long450.sesame.combatAdj = data.frame(t(apply(M.norm2, 1, getLongMM)))

long450.sesame.combatAdj$FDR_group = p.adjust(long450.sesame.combatAdj$group_pvalue, method="BH")
long450.sesame.combatAdj$FDR_clinage = p.adjust(long450.sesame.combatAdj$clinage_pvalue, method="BH")
long450.sesame.combatAdj$FDR_gender = p.adjust(long450.sesame.combatAdj$gender_pvalue, method="BH")

#long450.sesame.combatAdj$FDR_t_group = p.adjust(long450.sesame.combatAdj$`groupT1D control_t_pvalue`, method="BH")

#long450.sesame.combatAdj$Mcase = long450.sesame.combatAdj$`(Intercept)`
#long450.sesame.combatAdj$Mcontrol = long450.sesame.combatAdj$`(Intercept)` + long450.sesame.combatAdj$`groupT1D control` 

#long450.sesame.combatAdj$Bcontrol = getBeta(long450.sesame.combatAdj$Mcontrol)
#long450.sesame.combatAdj$Bcase = getBeta(long450.sesame.combatAdj$Mcase)

#long450.sesame.combatAdj$dfBeta = long450.sesame.combatAdj$Bcase - long450.sesame.combatAdj$Bcontrol

save(long450.sesame.combatAdj, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/long450.sesame.combatAdj.Rdata")

#########################################
# EPIC: Longitudinal Main Effects Model #
#########################################

rm(list=ls())
options(stringsAsFactors=FALSE)

library(plyr)
#library(bumphunter)
#library(minfi)
#library(ChAMP)
#library(waterRmelon)
library(doParallel)
#library(DMRcate)
#library(sva)
#library(rtracklayer)
library(nlme)

session =sessionInfo()
#save(session, file="/home/linuxshared/fromHome.v2/vanderll/Norris/T1Dpaper/sessionInfo_cellTypeAdjModels.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/phenotypeForLong/EPICpheno.long.Rdata")

load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameProcessing/sesameEPIC.batchAdj.Mmatrix.Rdata")


#subset for t1d 
pheno.t1d = pheno.forAnalysis[-which(pheno.forAnalysis$T1Dgroup==""),]
M.norm2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% pheno.t1d$rgName)]

table(colnames(M.norm2)==pheno.t1d$rgName)
#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm2)){
	order = c(order, which(pheno.t1d$rgName==colnames(M.norm2)[i]))
}
pheno.t1d = pheno.t1d[order,]
table(colnames(M.norm2)==pheno.t1d$rgName)

#remove any CWB times 
pheno.t1d.noCWB = pheno.t1d[-which(pheno.t1d$casegroup2=="CWB"),]
M.norm2 = M.norm2[,which(colnames(M.norm2) %in% pheno.t1d.noCWB$rgName)]

table(colnames(M.norm2)==pheno.t1d.noCWB$rgName)


#remove TR from 450K;
load(file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")
#actual sample TR
pheno.t1d.noTR = pheno.t1d.noCWB[-which(pheno.t1d.noCWB$key %in% long450.arrays$key),]
#subject level TR
pheno.t1d.final = pheno.t1d.noTR[-which(pheno.t1d.noTR$ID.x %in% as.matrix(long450.subjects)),]


M.norm2 = M.norm2[,which(colnames(M.norm2) %in% pheno.t1d.final$rgName)]

table(colnames(M.norm2)==pheno.t1d.final$rgName)

if(FALSE){
longEPIC.subjects = unique(pheno.t1d.final$ID.x)
longEPIC.arrays = data.frame(pheno.t1d.final[,c(1, 4,39,2)], "yes")
colnames(longEPIC.arrays) = c("key", "ID", "DOV", "Array", "T1DlongEPICanalysis")

save(longEPIC.subjects, longEPIC.arrays, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")

##get Fran's dataframe (columns wanted: "ID", "DOVisit", "methylation450Kdata", "methylation850K")
longEPIC = data.frame(longEPIC.arrays[,c(2,3)], "", longEPIC.arrays[,5])
long450 = data.frame(long450.arrays[,c(2,3,5)], "")

colnames(longEPIC) = c("ID", "DOVisit", "metnylation450Kdata", "methylation850Kdata")
colnames(long450) = c("ID", "DOVisit", "metnylation450Kdata", "methylation850Kdata")

forFran = rbind(long450, longEPIC)
write.csv(forFran, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/forFran_T1DlongAnalysisPhenotype.csv", row.names=FALSE)

write.csv(long450.subjects, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/T1Dlong.450K.subjects.csv", row.names=FALSE)
write.csv(longEPIC.subjects, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/T1Dlong.EPIC.subjects.csv", row.names=FALSE)
}

### Longitudinal Model ###

### T1D Case vs T1D Control ###
nCpG <- nrow(M.norm2)

group = as.factor(pheno.t1d.final$T1Dgroup)
t = as.factor(pheno.t1d.final$t)
ID = as.factor(pheno.t1d.final$ID.x)
gender=as.factor(pheno.t1d.final$SEX.x)
clinage=pheno.t1d.final$clinage

getBeta = function(a){
  b = 2^a/(1 + 2^a)  
  return(b)
}
#a = M.want[1,]

getLongMM = function(a){
	fit.mv2 <- try(lme (a ~  clinage+group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit))
	want = rep(NA, length=36)
	if("try-error" %in% class(fit.mv2)){want = want}else{
		
	#get ANOVA values;
	p2<- anova(fit.mv2)
	numDF = p2$numDF 
	names(numDF) = paste(rownames(p2), "numDF", sep="_")
	denDF = p2$denDF
	names(denDF) = paste(rownames(p2), "denDF", sep="_")
	Fvalue = p2$`F-value`	
	names(Fvalue) = paste(rownames(p2), "Fvalue", sep="_")
	pvalue = p2$`p-value`	
	names(pvalue) = paste(rownames(p2), "pvalue", sep="_")

	#get t-test values;
	t_pvalue = summary(fit.mv2)$tTable[,"p-value"]
	names(t_pvalue) = paste(rownames(summary(fit.mv2)$tTable), "t_pvalue", sep="_")
	t_stat = summary(fit.mv2)$tTable[,"t-value"]
	names(t_stat) = paste(rownames(summary(fit.mv2)$tTable), "t_stat", sep="_")
	t_stdError = summary(fit.mv2)$tTable[,"Std.Error"]
	names(t_stdError) = paste(rownames(summary(fit.mv2)$tTable), "t_stdError", sep="_")
	t_df = summary(fit.mv2)$tTable[,"DF"]
	names(t_df) = paste(rownames(summary(fit.mv2)$tTable), "t_df", sep="_")

	want = c(summary(fit.mv2)$coefficients$fixed, pvalue, numDF, denDF, Fvalue, t_pvalue, t_stat, t_stdError, t_df)
	
	}
	
	names(want) =  c("Intercept","clinage","groupT1D control","genderMale","Intercept_pvalue","clinage_pvalue","group_pvalue","gender_pvalue","Intercept_numDF","clinage_numDF","group_numDF","gender_numDF", "Intercept_denDF","clinage_denDF", "group_denDF","gender_denDF", "Intercept_Fvalue","clinage_Fvalue", "group_Fvalue","gender_Fvalue", "Intercept_t_pvalue","clinage_t_pvalue", "groupT1D control_t_pvalue","genderMale_t_pvalue", "Intercept_t_stat", "clinage_t_stat", "groupT1D control_t_stat","genderMale_t_stat", "Intercept_t_stdError","clinage_t_stdError", "groupT1D control_t_stdError", "genderMale_t_stdError", "Intercept_t_df" ,"clinage_t_df", "groupT1D control_t_df","genderMale_t_df")

	return(want)
}
 
 #test = M.norm2[356180:356190,]
#MMtest = data.frame(t(apply(test, 1, getLongMM)))
####YES!!!!! It works!!!!! Do on full dataset; 

longEPIC.sesame.combatAdj = data.frame(t(apply(M.norm2, 1, getLongMM)))

longEPIC.sesame.combatAdj$FDR_group = p.adjust(longEPIC.sesame.combatAdj$group_pvalue, method="BH")
longEPIC.sesame.combatAdj$FDR_clinage = p.adjust(longEPIC.sesame.combatAdj$clinage_pvalue, method="BH")
longEPIC.sesame.combatAdj$FDR_gender = p.adjust(longEPIC.sesame.combatAdj$gender_pvalue, method="BH")

#longEPIC.sesame.combatAdj$FDR_t_group = p.adjust(longEPIC.sesame.combatAdj$`groupT1D control_t_pvalue`, method="BH")

#longEPIC.sesame.combatAdj$Mcase = longEPIC.sesame.combatAdj$`(Intercept)`
#longEPIC.sesame.combatAdj$Mcontrol = longEPIC.sesame.combatAdj$`(Intercept)` + longEPIC.sesame.combatAdj$`groupT1D control` 

#longEPIC.sesame.combatAdj$Bcontrol = getBeta(longEPIC.sesame.combatAdj$Mcontrol)
#longEPIC.sesame.combatAdj$Bcase = getBeta(longEPIC.sesame.combatAdj$Mcase)

#longEPIC.sesame.combatAdj$dfBeta = longEPIC.sesame.combatAdj$Bcase - longEPIC.sesame.combatAdj$Bcontrol

save(longEPIC.sesame.combatAdj, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/longEPIC.sesame.combatAdj.Rdata")

