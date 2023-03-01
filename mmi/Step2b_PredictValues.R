# modeling metrics for reference sites 

source("bin/cbind.na.R")

results<-read.table("refcal.rsq.txt")
foo<-as.matrix(results)
foo2<-do.call(rbind, strsplit(as.character(foo[,2]), " "))
foo2<-as.data.frame(foo2[,c(1,3)])
names(foo2)<-c('top.metrics', "rsq")

foo2$rsq<-as.numeric(as.character(foo2$rsq))

if (diatoms ==T) { 
write.csv(foo2, 'diatoms.refcal.rsq.csv') }

if (sba ==T) { 
  write.csv(foo2, 'sba.refcal.rsq.csv') }

if (hybrid ==T) { 
  write.csv(foo2, 'hybrid.refcal.rsq.csv') }

metrics<-foo2$top.metrics
metrics<-droplevels(metrics)
#results$rsq<-as.numeric(as.character(results$rsq))
tops <-droplevels(subset(foo2, rsq>0.2)) # or pr2.1>0.1
top.metrics<-tops$top.metrics


if (diatoms == T ) {write.csv(top.metrics, "diatoms.top.metrics.csv") }
if (sba == T ) {write.csv(top.metrics, "sba.top.metrics.csv") }
if (hybrid == T ) {write.csv(top.metrics, "hybrid.top.metrics.csv") }

##########################################################
#using RF models, predict values 
##########################################################

predicted<-data.frame(stations$SeqID, stations$SeqID)
row.names(predicted)<-predicted[,1] 
predicted.out<-predicted

for (i in unique(metrics)) {        # i<-"cnt.spp.BCG2"
 
  # import RF model 
  filename<-paste0("R/mmi/", assemblage,"/",assemblage,".rf.models/",assemblage,".",i,"RF.model.Rdata" )
  load(filename)
  pred.var<-row.names(rf.out$importance)
  
  foo<-complete.cases(stations[,pred.var])
  stations1<-stations[foo,]
  
  ### apply new RF model to ALL data 
  stations1<-stations1[,pred.var]
  stations2<-as.data.frame(lapply(stations1, as.numeric))
  row.names(stations2)<-row.names(stations1)
  set.seed(10)
  pred<-as.data.frame(predict(rf.out, stations2))
  names(pred)[1]<-i
  predicted<-cbind.na(predicted, pred) 
  predicted.out<-cbind.na(predicted.out, pred) 
  } 

predicted.out<-predicted.out[,-1]
colnames(predicted.out)<-unique(metrics)

if (diatoms == T ) {write.csv(predicted.out, "diatoms.predicted.metrics.csv") }
if (sba == T ) {write.csv(predicted.out, "sba.predicted.metrics.csv") }
if (hybrid == T ) {write.csv(predicted.out, "hybrid.predicted.metrics.csv") }

predicted.top<-predicted.out[,top.metrics]

if (diatoms == T ) {write.csv(predicted.top, "diatoms.predicted.top.metrics.csv") }
if (sba == T ) {write.csv(predicted.top, "sba.predicted.top.metrics.csv") }
if (hybrid == T ) {write.csv(predicted.top, "hybrid.predicted.top.metrics.csv") }






