# mASCI


library(reshape)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(reshape2)
my_palette<-brewer.pal(11, "Spectral")
palette(my_palette)
# library(dplyr)
# library(plyr)
library(funrar)
# randomForest, speakeasy
# EnvironmentalDNA index == inverse of wisconsin index 
# proportions in sample, for each species, highest proportion, removes amplification efficiency

#  env -------------------------------------------
env<-read.csv("data/env/env5.csv", header=T, stringsAsFactors = F)
env[env==""] <- NA
env[env=="#N/A"] <- NA

foo<-which(env$Replicate %in% c("MB", "eDNA", "LB"))
env<-env[-foo,]

env$SampleID<-paste0(env$stationcode.correct, "_", env$Date)

env$AFDM<-as.numeric(as.character(env$AFDM_Algae..Particulate))
env$ChlA<-as.numeric(as.character(env$Chlorophyll.a..Particulate))
env$DOC<-as.numeric(as.character(env$Dissolved.Organic.Carbon..Dissolved))
env$SpCond<-as.numeric(as.character(env$SpecificConductivity..Total))

## new chem
# chem<-read.csv("~/Documents/R/SMC/SMC_data/unified_chem_filtered.csv")
# foo<-which(chem$stationcode %in% env$stationcode)
# chem<-chem[foo,]
# chem$SampleID<-paste0(chem$stationcode, "_", chem$sampledate)
# 
# sort(unique(chem$analytename))
# 
# keep<-c(#"Ash Free Dry Mass", 
#         "AFDM_Algae, Particulate",
#         #"Chlorophyll a", 
#         #"Chlorophyll a, Total", 
#         "Chlorophyll a, Particulate", 
#         #"Dissolved Organic Carbon",
#         "Dissolved Organic Carbon, Dissolved",
#         "SpecificConductivity, Total" 
#         )
# foo<-which(chem$analytename %in% keep)
# chem<-chem[foo,]
# 
# chem<-as.data.frame(acast(chem, SampleID~analytename,
#                             value.var="result",
#                             fun.aggregate=sum))
# write.csv(chem, "data/env/extrachem.csv")
#chem<-read.csv("data/env/extrachem.csv")

#foo<-which(chem$SampleID %in% env$SampleID)

#env<-merge(env, chem, by="SampleID", all.x=T, all.y=F, sort=F)


# chem<-read.csv("data/env/chem/unified_chem_filtered.csv")
# unique(sort(chem$analytename))
# keep<-c("AFDM_Algae, Total" , 
#         "Chlorophyll a, Total",
#         "Dissolved Organic Carbon, Total", 
#         "SpecificConductivity, Total")
# foo<-which(chem$analytename %in% keep)
# chem.sub<-chem[foo,]
# foo<-which(chem.sub$stationcode %in% env$masterid)
# chem.sub<-chem.sub[foo,]
# unique(chem.sub$analytename)
# chem.sub<-droplevels(chem.sub)
# write.csv(chem.sub, "data/env/chem/chem.sub.csv")


