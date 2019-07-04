###########################################################
# DVP Analysis NEW                                        #
#                                                         #
# This is for new cross-sections (CWB, preSV and postSV)  #
# Date: 10/19/2017                                        #
# Author: Lauren Vanderlinden                             #
###########################################################

#Set up workspace 
rm(list=ls())
options(stringsAsFactors = FALSE)


library(ggplot2)
library(RColorBrewer)
library(plyr)
library(limma)
library(bumphunter)
library(doParallel)
library(sva)
library(nlme)
library(rtracklayer)
library(car)
library(knitr)

#load methylation data;
load("/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/forAnalysis/Mixedmodel.SWANnorm.Rdata")

#load in new phenotype classification from Fran 
newGrouping = read.csv(file="/home/vanderll/Norris/prepost group.LVdate.csv")

#so the date in Fran's file doesn't force the month and day to be 2 digits, wheras the date in the master list does.  This is causing issues merging.  
newGrouping$year = as.numeric(sapply(strsplit(newGrouping$DOVISIT, split="/", fixed=TRUE), "[[", 3))
newGrouping$month = sapply(strsplit(newGrouping$DOVISIT, split="/", fixed=TRUE), "[[", 1)
newGrouping$day = sapply(strsplit(newGrouping$DOVISIT, split="/", fixed=TRUE), "[[", 2)


DOVISITnew = c()
for(i in 1:nrow(newGrouping)){
  if(newGrouping[i, "year"]<10){DOVISITnew[i]=paste(newGrouping[i,"month"], newGrouping[i,"day"], paste("200", newGrouping[i,"year"], sep=""), sep="/")}
  if(newGrouping[i, "year"]>9&newGrouping[i, "year"]<50){DOVISITnew[i]=paste(newGrouping[i,"month"], newGrouping[i,"day"], paste("20", newGrouping[i,"year"], sep=""), sep="/")}
  if(newGrouping[i, "year"]>50){DOVISITnew[i]=paste(newGrouping[i,"month"], newGrouping[i,"day"], paste("19", newGrouping[i,"year"], sep=""), sep="/")}
}

newGrouping$DOVISITnew = DOVISITnew

newGrouping.2merge = newGrouping[,c("ID", "DOVISITnew", "casegroup4")]

pheno.want = merge(pheno.v1, newGrouping.2merge, by.x=c("IVYOmicsMethylationMaster_ID", "DOVISIT"), by.y=c("ID", "DOVISITnew"))


### DVP: Longitudinal Model ###

#get specific arrays for T1D paper 
M.mm = M.norm.batch.mm[,which(colnames(M.norm.batch.mm) %in% as.matrix(pheno.want$Array))]

table(colnames(M.mm) == pheno.want$Array)

order=c()
for(i in 1:nrow(pheno.want)){
  order = c(order, which(colnames(M.mm)==pheno.want[i,"Array"]))
}
M.mm = M.mm[,order]
table(colnames(M.mm) == pheno.want$Array)

ID = as.matrix(pheno.want$ID)
clinage = as.matrix(pheno.want$clinage)
group4 = as.matrix(pheno.want$group4)
sex = as.matrix(pheno.want$SEX)

getLongitudinalDVPresults = function(a){
  #model = lme (a ~  clinage+ sex, random=~0+clinage|ID, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8), method = "REML")
  model <- gls(a ~ clinage + sex, method="REML", correlation=corAR1(,form=~1|ID))
  
  #get residuals;
  resids = residuals(model)
  #perform the mean-trimmed Levene's test for homogeneous variances
  dvp.test = leveneTest(resids, as.factor(group4), center=mean, trim=0.1)
  dvp.Fvalue = dvp.test[1,2]
  dvp.pvalue = dvp.test[1,3]
  
  #put the statistics we want in 1 object
  results = cbind(dvp.Fvalue, dvp.pvalue)
  return(results)
}

#perform the DVP on all probes
DVP.longitudinal = as.data.frame(t(apply(M.mm, 1, getLongitudinalDVPresults)))
colnames(DVP.longitudinal) = c("dvp.Fvalue", "dvp.pvalue")
#adjust for multiple comparisons
DVP.longitudinal$FDR = p.adjust(DVP.longitudinal[,"dvp.pvalue"], method="BH")

### DVP: CWB Model ###

M.CWB = M.mm[,which(pheno.want$casegroup4=="CWB")]
pheno.CWB = pheno.want[which(pheno.want$casegroup4=="CWB"),]

clinage.cwb = as.matrix(pheno.CWB$clinage)
group4.cwb = as.matrix(pheno.CWB$group4)
sex.cwb = as.matrix(pheno.CWB$SEX)

getCWB.DVPresults = function(a){
  #perform the original test (in this case cross-sectional at timepoint T1D)
  model = lm(a~ sex.cwb + clinage.cwb)
  #get residuals;
  resids = residuals(model)
  #perform the mean-trimmed Levene's test for homogeneous variances
  dvp.test = leveneTest(resids, as.factor(group4.cwb), center=mean, trim=0.1)
  dvp.Fvalue = dvp.test[1,2]
  dvp.pvalue = dvp.test[1,3]
  
  #put the statistics we want in 1 object
  results = cbind(dvp.Fvalue, dvp.pvalue)
  return(results)
}

DVP.CWB.crossSectional = as.data.frame(t(apply(M.CWB, 1, getCWB.DVPresults)))
colnames(DVP.CWB.crossSectional) = c("dvp.Fvalue", "dvp.pvalue")
#adjust for multiple comparisons
DVP.CWB.crossSectional$FDR = p.adjust(DVP.CWB.crossSectional[,"dvp.pvalue"], method="BH")

### DVP: preSV Model ###

M.preSV = M.mm[,which(pheno.want$casegroup4=="Pre SV")]
pheno.preSV = pheno.want[which(pheno.want$casegroup4=="Pre SV"),]

clinage.preSV = as.matrix(pheno.preSV$clinage)
group4.preSV = as.matrix(pheno.preSV$group4)
sex.preSV = as.matrix(pheno.preSV$SEX)

getpreSV.DVPresults = function(a){
  #perform the original test (in this case cross-sectional at timepoint T1D)
  model = lm(a~ sex.preSV + clinage.preSV)
  #get residuals;
  resids = residuals(model)
  #perform the mean-trimmed Levene's test for homogeneous variances
  dvp.test = leveneTest(resids, as.factor(group4.preSV), center=mean, trim=0.1)
  dvp.Fvalue = dvp.test[1,2]
  dvp.pvalue = dvp.test[1,3]
  
  #put the statistics we want in 1 object
  results = cbind(dvp.Fvalue, dvp.pvalue)
  return(results)
}

DVP.preSV.crossSectional = as.data.frame(t(apply(M.preSV, 1, getpreSV.DVPresults)))
colnames(DVP.preSV.crossSectional) = c("dvp.Fvalue", "dvp.pvalue")
#adjust for multiple comparisons
DVP.preSV.crossSectional$FDR = p.adjust(DVP.preSV.crossSectional[,"dvp.pvalue"], method="BH")



### DVP: postSV Model ###

M.postSV = M.mm[,which(pheno.want$casegroup4=="Post SV")]
pheno.postSV = pheno.want[which(pheno.want$casegroup4=="Post SV"),]

clinage.postSV = as.matrix(pheno.postSV$clinage)
group4.postSV = as.matrix(pheno.postSV$group4)
sex.postSV = as.matrix(pheno.postSV$SEX)

getpostSV.DVPresults = function(a){
  #perform the original test (in this case cross-sectional at timepoint T1D)
  model = lm(a~ sex.postSV + clinage.postSV)
  #get residuals;
  resids = residuals(model)
  #perform the mean-trimmed Levene's test for homogeneous variances
  dvp.test = leveneTest(resids, as.factor(group4.postSV), center=mean, trim=0.1)
  dvp.Fvalue = dvp.test[1,2]
  dvp.pvalue = dvp.test[1,3]
  
  #put the statistics we want in 1 object
  results = cbind(dvp.Fvalue, dvp.pvalue)
  return(results)
}

DVP.postSV.crossSectional = as.data.frame(t(apply(M.postSV, 1, getpostSV.DVPresults)))
colnames(DVP.postSV.crossSectional) = c("dvp.Fvalue", "dvp.pvalue")
#adjust for multiple comparisons
DVP.postSV.crossSectional$FDR = p.adjust(DVP.postSV.crossSectional[,"dvp.pvalue"], method="BH")


### save results;
save(DVP.longitudinal, DVP.CWB.crossSectional,DVP.postSV.crossSectional,DVP.preSV.crossSectional, file="/home/vanderll/Norris/DVP/DVPresults.newTimeGroups.Rdata")

