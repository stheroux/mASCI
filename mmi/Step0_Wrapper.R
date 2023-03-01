# wrapper script for creating MMIs 

# to do 
# 1. clean up scripts 
# 2. check otu5.rbcL 
# 3. run performance scripts for output 
# 4. index scores when ind metrics are NA 



diatoms = T # one at a time
sba = F
hybrid = F
ASV = T

{ 
###Loading Needed Packages
library(vegan)
library(indicspecies)
library(reshape2)
library(tidyverse)
library(plyr) # WARNING: problems arise if you have dplyr already loaded, so maybe start from scratch
options(gsubfn.engine = "R")
library(sqldf)
library(ggplot2)
library(pscl)
library(caret)
library(randomForest)
library(ggplot2)

if (diatoms == T) { assemblage <- "diatoms" }
if (sba == T) { assemblage <- "sba" }
if (hybrid == T) { assemblage <- "hybrid" }

}

# get OTU ---------------

if((diatoms == T) & (ASV==F)) {
  otu=otu5.rbcL
  foo<-as.data.frame(do.call(rbind, strsplit(as.character(otu$ConsensusLineage), ";")))
  colnames(foo)<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus", "Species")
  foo$Species<-trimws(foo$Species, which = "both")
  otu<-cbind(otu, foo) }

if((diatoms == T) & (ASV==T)) {
  otu=otu5.rbcL; otu$Species <- otu$OTUID }


# Step 1 
source("R/mmi/Step1_EvalMetrics_4.R")

# Step 2a 
source("R/mmi/Step2a_ModelMetrics3.R")

# Step 2b 
source("R/mmi/Step2b_PredictValues.R")

# Step 3
source("R/mmi/Step3_ScoreMetrics.R")

# Step 4a 
source("R/mmi/Step4a_ScreenMetrics_6.R")

# Step 4b
source("R/mmi/Step4b_iterate.metrics_all.R") # requires entering "x" 

# Step 4c save winning metrics
source("R/mmi/Step4c_WinningMetrics.R")

# Step 4d
source("R/mmi/Step4d_Winning_PMMIs.R")

# following steps are done on all assemblages at the same time 

# Step 5 
#source("Step5_ASCI_Scores_3.R")
#setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 6
#source("Step6_PullRFmodels.R")
#setwd("~/Documents/R/ASCI/pMMI_newnew/")

# Step 7 
source("R/mmi/Step7_Compare.R")

################################
# 

# Performance 1 
#source("~/Documents/R/ASCI/PERFORMANCE_newnew/001_Performance14.R")
#source("~/Documents/R/ASCI/PERFORMANCE_newnew/002_Figures10.R")
#source("~/Documents/R/ASCI/PERFORMANCE_newnew/003_CompareCSCI.R")
#source("~/Documents/R/ASCI/PERFORMANCE_newnew/004_Fig5_LowE2.R")



