#####################################
# 450K DATA 						#
#####################################


## Setup Workspace
rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("/home/linuxshared/vanderll/Norris")

library(minfi)
library(wateRmelon)
library(ChAMP)
library(RColorBrewer)
library(sva)
library(plyr)
library(ChAMP, lib.loc="/home/vanderll/R/x86_64-redhat-linux-gnu-library/3.2")
library("FlowSorted.Blood.450k")

session = sessionInfo()
save(session, file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/RsessionInfo.Rdata")

#########################
# Step 1: Read in Data 	#
#########################

### CLINICAL PHENOTYPE DATA ###
pheno = read.csv(file="/home/linuxshared/vanderll/Norris/Data/BasicFormattingUpdates/IVYOmicsMethylationDemo.savedAsCSV.csv")
dim(pheno)
head(pheno)
length(unique(pheno$NID))

table(pheno$group2)
table(pheno$SEX)

missingGroup.sampleIDs = pheno[which(pheno$group2==""), "NID"]

pheno[which(pheno$NID %in% missingGroup.sampleIDs),]


### METHYLATION (ARRAY) DATA ###
##look for the methylation data

#First look at the pilot array
baseDirPilot = c("/home/linuxshared/vanderll/Norris/DataRaw/Pilot4_09222015")
targetPilot = read.450k.sheet(baseDirPilot)
targetPilot.forUse = targetPilot[-which(targetPilot$Basename=="character(0)"),]
targetPilot.forUse$Sample_Plate = gsub("Plate 1", "Plate Extra4", targetPilot.forUse$Sample_Plate)
 

test.baseDir = c("/home/linuxshared/vanderll/Norris/DataRaw/pilotPlate")
targetPilot.24.forUse = read.450k.sheet(test.baseDir )
targetPilot.24.forUse$Sample_Plate = gsub("Plate 1", "Plate Pilot", targetPilot.24.forUse$Sample_Plate)

test.baseDir.confirm = c("/home/linuxshared/vanderll/Norris/DataRaw/confirmationArrays")
targetConfirm.forUse = read.450k.sheet(test.baseDir.confirm)
targetConfirm.forUse$Sample_Plate = rep("Confirmation", times=nrow(targetConfirm.forUse))


 #Now look at the standard 7 arrays
baseDir = c("/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate1", "/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate2_09302015", "/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_3_10162015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_4_10162015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_5_10292015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_6_12022015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_7_12022015")

targets <- data.frame()
for(i in 1:length(baseDir)){
	targets = rbind(targets, read.450k.sheet(baseDir[i]))
}
nrow(targets)
length(unique(targets$Sample_Name))

#include the pilot arrays in target list
targets.wPilot = rbind(targets, targetPilot.forUse, targetPilot.24.forUse, targetConfirm.forUse)

nrow(targets.wPilot)
length(unique(targets.wPilot$Sample_Name))

table(targets.wPilot$Sample_Plate)

#actually read in the data
##RUN THIS: rgSet <- read.450k.exp(targets=targets.wPilot, extended=T)

##Need to have extended=FALSE in order for estimating proportions of cells
#rgSet <- read.450k.exp(targets=targets.wPilot, extended=FALSE)
#start time 9:39 am...
#end time 9:58 am.   

##RUN THIS: save(rgSet, file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet.Rdata")
#took about 10 minutes to save



### COMBINE CLINICAL PHENOTYPES & METHYLATION ARRAY ###
#get the same ID in both the clinical data and the methylation data
load("/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet.Rdata")

colnames(rgSet)

#rgSet sample names are slide_sample
targets.wPilot$rgName = paste(targets.wPilot$Slide, targets.wPilot$Array, sep="_")
table(targets.wPilot$rgName==colnames(rgSet))

toOrder = cbind(targets.wPilot, c(1:nrow(targets.wPilot)))
colnames(toOrder)[10] = "ArrayOrder"

pheno.v2 = merge(toOrder, pheno, by.x="Sample_Name", by.y="NID", all=TRUE)

#sort clinic dataset into the order of methylation column order. 
#pheno.v3 = pheno.v2[!is.na(pheno.v2$rgName),]
pheno.v3 = pheno.v2[order(pheno.v2$ArrayOrder),]
#double check that the order is correct;
table(pheno.v3$rgName==colnames(rgSet))

stopifnot(all(pheno.v3$Sample_Name==rgSet$Sample_Name))
pData(rgSet)$Sex <- pheno.v3$SEX 
pData(rgSet)$group2  <- pheno.v3$group2 
pData(rgSet)$Sample_Plate <- pheno.v3$Sample_Plate
pData(rgSet)$Sample_Well <- pheno.v3$Sample_Well
pData(rgSet)$Slide <- pheno.v3$Slide
pData(rgSet)$rgName <- pheno.v3$rgName
pData(rgSet)$ID <- pheno.v3$Sample_Name
pData(rgSet)$IVYomicID <- pheno.v3$IVYOmicsMethylationMaster_ID
pData(rgSet)$labID <- pheno.v3$Labid

pd <- pData(rgSet) 
dim(pd)
head(pd)

save(pd, rgSet, file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet_and_phenoData.Rdata")
#Jeremy must have changed my memory allotment on yampa becuase this took almost 30 min.



#################################
# EPIC ARRAYS					#
#################################


## Setup Workspace
rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("/home/linuxshared/vanderll/Norris/Reports/EPIC805K.processingAndQC/")

library(minfi)
library(wateRmelon)
#library(ChAMP)
library(RColorBrewer)
library(sva)
library(plyr)
library(DMRcate)
#library("FlowSorted.Blood.450k")

session = sessionInfo()
save(session, file="/home/linuxshared/vanderll/Norris/Reports/EPIC805K.processingAndQC/RsessionInfo.Rdata")

#########################
# Step 1: Read in Data 	#
#########################

#### SKIP THE PHENOTYPE DATA TILL I HEAR BACK FROM RANDI ####
### CLINICAL PHENOTYPE DATA ###
pheno = read.csv(file="/home/linuxshared/vanderll/Norris/Data/BasicFormattingUpdates/IVYOmicsMethylationDemo.savedAsCSV.csv")
dim(pheno)
head(pheno)
length(unique(pheno$NID))

table(pheno$group2)
table(pheno$SEX)

missingGroup.sampleIDs = pheno[which(pheno$group2==""), "NID"]

pheno[which(pheno$NID %in% missingGroup.sampleIDs),]

### METHYLATION (ARRAY) DATA ###
##look for the methylation data

#First look at the pilot array
baseDirPilot = c("/home/linuxshared/vanderll/Norris/DataRaw/850Karrays")
targetPilot = read.metharray.sheet(baseDirPilot)

dim(targetPilot)
head(targetPilot)
table(targetPilot$Sample_Plate)
length(unique(targetPilot$Basename))

#need to change plate name in the pilot sample sheet;
targetPilot.forUse = targetPilot
targetPilot.forUse$Sample_Plate = targetPilot$Sample_Plate = gsub("Plate 1", "Plate 8", targetPilot.forUse$Sample_Plate)

##### EXPLORE THIS ####
###Something looks weird with names;
table(targetPilot.forUse$Slide)
#there should only be 8/slide

arrayNames = paste(targetPilot$Slide, targetPilot$Array, sep="_")
length(arrayNames)
length(unique(arrayNames))


###If you have the other arrays and want to loop them in (and they already have the plate named as you would like...);

#Now look at the standard 7 arrays
baseDir = c("/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate1", "/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate2_09302015", "/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_3_10162015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_4_10162015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_5_10292015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_6_12022015",
"/home/linuxshared/vanderll/Norris/DataRaw/JNorris_plate_7_12022015")

targets <- data.frame()
for(i in 1:length(baseDir)){
	targets = rbind(targets, read.metharray.sheet(baseDir[i]))
}
nrow(targets)
length(unique(targets$Sample_Name))

#include the pilot arrays in target list
targets.wPilot = rbind(targets, targetPilot.forUse)

nrow(targets.wPilot)
length(unique(targets.wPilot$Sample_Name))

table(targets.wPilot$Sample_Plate)

###actually read in the idat data files

####As of May 9, 2017: We just want to process the pilot data
targets.wPilot = targetPilot

##Note: Need to have extended=FALSE in order for estimating proportions of cells using Houseman
rgSet <- read.metharray.exp(targets=targets.wPilot, extended=TRUE)
#about 30 minutes

#issues with read.metharray.exp identifying what type of array we have.  Go ahead and force the array type annotation:
rgSet@annotation=c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19")

save(rgSet, file="rgSet.Rdata")
#took about 10 minutes to save.


### COMBINE CLINICAL PHENOTYPES & METHYLATION ARRAY ###
#get the same ID in both the clinical data and the methylation data
load("/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet.Rdata")

colnames(rgSet)

#rgSet sample names are slide_sample
targets.wPilot$rgName = paste(targets.wPilot$Slide, targets.wPilot$Array, sep="_")
table(targets.wPilot$rgName==colnames(rgSet))

toOrder = cbind(targets.wPilot, c(1:nrow(targets.wPilot)))
colnames(toOrder)[10] = "ArrayOrder"

pheno.v2 = merge(toOrder, pheno, by.x="Sample_Name", by.y="NID", all=TRUE)

#sort clinic dataset into the order of methylation column order. 
#pheno.v3 = pheno.v2[!is.na(pheno.v2$rgName),]
pheno.v3 = pheno.v2[order(pheno.v2$ArrayOrder),]
#double check that the order is correct;
table(pheno.v3$rgName==colnames(rgSet))

stopifnot(all(pheno.v3$Sample_Name==rgSet$Sample_Name))
pData(rgSet)$Sex <- pheno.v3$SEX 
pData(rgSet)$group2  <- pheno.v3$group2 
pData(rgSet)$Sample_Plate <- pheno.v3$Sample_Plate
pData(rgSet)$Sample_Well <- pheno.v3$Sample_Well
pData(rgSet)$Slide <- pheno.v3$Slide
pData(rgSet)$rgName <- pheno.v3$rgName
pData(rgSet)$ID <- pheno.v3$Sample_Name
pData(rgSet)$IVYomicID <- pheno.v3$IVYOmicsMethylationMaster_ID
pData(rgSet)$labID <- pheno.v3$Labid

pd <- pData(rgSet) 
dim(pd)
head(pd)

save(pd, rgSet, file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet_and_phenoData.Rdata")
#This took almost 30 min!
