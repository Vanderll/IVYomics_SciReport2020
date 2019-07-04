########################################
# 450K: Longitudinal Interaction Model #
########################################


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

load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/final450Klong.arraysSubjects.Rdata")

## Get same longitudinal subjects ###

pheno.t1d.final = pheno.forAnalysis[which(pheno.forAnalysis$Array %in% long450.arrays$Array),]
M.norm2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% long450.arrays$Array)]
table(colnames(M.norm2)==pheno.t1d.final$Array)

#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm2)){
	order = c(order, which(pheno.t1d.final$Array==colnames(M.norm2)[i]))
}
pheno.t1d.final = pheno.t1d.final[order,]
table(colnames(M.norm2)==pheno.t1d.final$Array)


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
a = M.norm2[1,]
test = lme (a ~  clinage*group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit)
getLongMM = function(a){
	fit.mv2 <- try(lme (a ~  clinage*group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit))
	want = rep(NA, length=45)
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
	
	names(want) = c("Intercept", "clinage", "groupT1D_control", "genderMale", "clinage:groupT1D_control", "Intercept_pvalue", "clinage_pvalue", "group_pvalue", "gender_pvalue", "clinage:group_pvalue", "Intercept_numDF", "clinage_numDF", "group_numDF", "gender_numDF", "clinage:group_numDF", "Intercept_denDF", "clinage_denDF", "group_denDF", "gender_denDF", "clinage:group_denDF", "Intercept_Fvalue", "clinage_Fvalue", "group_Fvalue", "gender_Fvalue", "clinage:group_Fvalue", "Intercept_t_pvalue", "clinage_t_pvalue", "groupT1D_control_t_pvalue", "genderMale_t_pvalue", "clinage:groupT1D_control_t_pvalue", "Intercept_t_stat", "clinage_t_stat", "groupT1D_control_t_stat", "genderMale_t_stat", "clinage:groupT1D_control_t_stat", "Intercept_t_stdError", "clinage_t_stdError", "groupT1D_control_t_stdError", "genderMale_t_stdError", "clinage:groupT1D_control_t_stdError", "Intercept_t_df", "clinage_t_df", "groupT1D_control_t_df", "genderMale_t_df", "clinage:groupT1D_control_t_df")

	return(want)
}
 
#MMtest = M.norm2[356180:356190,]
#MMtest = data.frame(t(apply(MMtest, 1, getLongMM)))
###YES!!!!! It works!!!!! Do on full dataset; 

long450.sesame.interaction = data.frame(t(apply(M.norm2, 1, getLongMM)))

save(long450.sesame.interaction, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/interaction/long450.sesame.interaction.Rdata")

###########################
# EPIC: Interaction Model #
###########################

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

load("/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/finalEPIClong.arraysSubjects.Rdata")

## Get same longitudinal subjects ###

pheno.t1d.final = pheno.forAnalysis[which(pheno.forAnalysis$rgName %in% longEPIC.arrays$Array),]
M.norm2 = M.sesame.batch[,which(colnames(M.sesame.batch) %in% longEPIC.arrays$Array)]
table(colnames(M.norm2)==pheno.t1d.final$rgName)

#reorder pheno file so arrays match up;
order = c()
for(i in 1:ncol(M.norm2)){
	order = c(order, which(pheno.t1d.final$rgName==colnames(M.norm2)[i]))
}
pheno.t1d.final = pheno.t1d.final[order,]
table(colnames(M.norm2)==pheno.t1d.final$rgName)


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
a = M.norm2[1,]
test = lme (a ~  clinage*group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit)
getLongMM = function(a){
	fit.mv2 <- try(lme (a ~  clinage*group+gender, random=~1|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML",  na.action=na.omit))
	want = rep(NA, length=45)
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
	
	names(want) = c("Intercept", "clinage", "groupT1D_control", "genderMale", "clinage:groupT1D_control", "Intercept_pvalue", "clinage_pvalue", "group_pvalue", "gender_pvalue", "clinage:group_pvalue", "Intercept_numDF", "clinage_numDF", "group_numDF", "gender_numDF", "clinage:group_numDF", "Intercept_denDF", "clinage_denDF", "group_denDF", "gender_denDF", "clinage:group_denDF", "Intercept_Fvalue", "clinage_Fvalue", "group_Fvalue", "gender_Fvalue", "clinage:group_Fvalue", "Intercept_t_pvalue", "clinage_t_pvalue", "groupT1D_control_t_pvalue", "genderMale_t_pvalue", "clinage:groupT1D_control_t_pvalue", "Intercept_t_stat", "clinage_t_stat", "groupT1D_control_t_stat", "genderMale_t_stat", "clinage:groupT1D_control_t_stat", "Intercept_t_stdError", "clinage_t_stdError", "groupT1D_control_t_stdError", "genderMale_t_stdError", "clinage:groupT1D_control_t_stdError", "Intercept_t_df", "clinage_t_df", "groupT1D_control_t_df", "genderMale_t_df", "clinage:groupT1D_control_t_df")

	return(want)
}
 
#test = M.norm2[356180:356190,]
#MMtest = data.frame(t(apply(test, 1, getLongMM)))
####YES!!!!! It works!!!!! Do on full dataset; 

longEPIC.sesame.interaction = data.frame(t(apply(M.norm2, 1, getLongMM)))

save(longEPIC.sesame.interaction, file="/home/linuxshared/vanderll/fromHome.v2/Norris/T1Dpaper/data/sesameResults_v2/interaction/longEPIC.sesame.interaction.Rdata")



