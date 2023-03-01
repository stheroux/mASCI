# modeling metrics for reference sites 

# Import metric observed scores 

if (diatoms == T) { 
obs<-read.csv("R/mmi/diatoms/diatoms.all.metrics.csv", row.names=1)
obs<-obs[order(row.names(obs)),]
dir.create("R/mmi/diatoms/diatoms.rf.models/")
  }

if (sba == T) { 
  obs<-read.csv("R/mmi/sba/sba.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  dir.create("R/mmi/sba/sba.rf.models/")
  }

if (hybrid == T) { 
  obs<-read.csv("R/mmi/hybrid/hybrid.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  dir.create("R/mmi/hybrid/hybrid.rf.models/")
  }

# import stations data 
stations<-stations[order(stations$SeqID),]
drop<-"XEFE"
foo<-which(stations$SeqID %in% drop)
if(length(foo)>0) {stations<-stations[-foo,]}
stations$new_lat<-stations$Lat
stations$new_long<-stations$Long
stations_tl<-stations
names(stations_tl)<-tolower(names(stations_tl))
stations_tl$new_lat<-stations_tl$lat
stations_tl$new_long<-stations_tl$long
row.names(stations)<-stations$SeqID

source("bin/cand.var.bu.R")
cand.var<-tolower(cand.var)
cand.var %in% names(stations) #do you have them all?

cal<-subset(stations, stations$CalVal=="Cal")
val<-subset(stations, stations$CalVal=="Val")

stations.rc<-subset(stations, (stations$SeqID) %in% (cal$SeqID))
stations.rv<-subset(stations, (stations$SeqID) %in% val$SeqID)
obs.rc<-subset(obs, row.names(obs) %in% cal$SeqID)
obs.rv<-subset(obs, row.names(obs) %in% val$SeqID)

stations.rc$SeqID == row.names(obs.rc)
stations.rv$SeqID == row.names(obs.rv)

stations.rc<-droplevels(stations.rc)

#foo<-which(is.na(colSums(obs.rc)))
#if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }
#foo<-which(colSums(obs.rc)==0)
#if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }

metrics.list<-sort(colnames(obs.rc))
#metrics.list<-c("cnt.spp.BCG2", "cnt.spp.BCG3") 

##################################################
##### BEGIN LOOP #######
##### this loop will build models to predict metric distributions #######


source("R/mmi/1model.loop.R")
#source("R/mmi/1screen.metrics.R")
library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
base<-model.loop
wd<-"R/mmi/diatoms/"
clusterExport(cl, c("base",
                    "metrics.list", 
                    "obs.rc", "stations.rc", 
                    "ctrl", "assemblage", "cand.var", "wd" ))

parSapply(cl, metrics.list,
          function(x)
            fun=base(x))

stopCluster(cl)


##### END LOOP #######
######################### 




