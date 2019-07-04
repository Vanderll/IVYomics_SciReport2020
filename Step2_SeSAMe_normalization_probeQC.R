#############
# 450K DATA #
#############

##install packages
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sesame")

##load packages
library(sesameData)
library(sesame)

### Load minfi object data

#do this part in console on yampa...
if(FALSE){
load(file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet_and_phenoData.Rdata")
rgSet
rgSet = updateObject(rgSet)
save(rgSet, pd, file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet_and_phenoData_updatedObject.Rdata")
}

load(file="/home/linuxshared/vanderll/Norris/Reports/v2.processingAndQC/rgSet_and_phenoData_updatedObject.Rdata")

### Make into sesame object
rgSet.sesame = sesamize(rgSet, parallel=TRUE)
save(rgSet.sesame, file="/home/linuxshared/vanderll/Norris/Reports/v3.sesame/rgSet_450K.Rdata")

#############
# EPIC DATA #
#############

##install packages
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sesame")

##load packages
library(sesameData)
library(sesame)

### Load minfi object data 

load(file="/home/linuxshared/vanderll/Norris/Reports/EPIC805K.processingAndQC/all850Kplates/rgSet_and_phenoData.Rdata")

### Make into sesame object
rgSet.sesame = sesamize(rgSet, parallel=TRUE)
save(rgSet.sesame, file="/home/linuxshared/vanderll/Norris/Reports/v3.sesame/rgSet_EPIC.Rdata")

