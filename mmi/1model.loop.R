
model.loop<-function(i) { 
  library(pscl)
  library(caret)
  library(randomForest)
  source("bin/OE.load.and.source.R")
  source("bin/OE.caret.load.and.source.R")
  source("bin/cand.var.bu.R")
  
  cand.var<-tolower(cand.var)

### run CARET to find optimal vars 
#setwd(wd)
stations.rc.sub<-stations.rc[,cand.var]
foo<-complete.cases(stations.rc.sub)
stations.rc.sub<-stations.rc.sub[foo,]
obs.rc.sub<-subset(obs.rc, row.names(obs.rc) %in% row.names(stations.rc.sub))
foo<-complete.cases(obs.rc.sub)
obs.rc.sub<-obs.rc.sub[foo,]

if(length(row.names(obs.rc.sub))>0) {
stations.rc.sub<-subset(stations.rc.sub, row.names(stations.rc.sub) %in% row.names(obs.rc.sub))

stations.rc.sub<-as.data.frame(lapply(stations.rc.sub, as.numeric))

set.seed(10)
rfe.out<-rfe(y=obs.rc.sub[,i], x=stations.rc.sub[,cand.var], 
             sizes=2:10,  rfeControl=ctrl, maximize=T, na.action=na.omit) 

# record optvars
optvars.rfe <- rfe.out$optVariables
optvars.rfe

### build RF with optimal vars from caret 
set.seed(10)
rf.out<-randomForest(y=obs.rc.sub[,i], x=stations.rc.sub[,optvars.rfe], 
                     ntree=500,importance=TRUE, norm.votes=TRUE, 
                     keep.forest=TRUE,na.action=na.omit, proximity=T)
saveas<-paste0("R/mmi/",assemblage,"/",assemblage,".rf.models/",assemblage, ".", i, "RF.model.Rdata")
save(rf.out, file=saveas)

rsq<-as.data.frame(rf.out$rsq)
rsq<-(rsq[500,]) 

### apply new RF model to CAL data 
set.seed(10)
predicted<-predict(rf.out, stations.rc.sub[,optvars.rfe])
observed<-obs.rc.sub[,i]

set.seed(10)
fit1<-glm(predicted~observed)
pr2<- as.data.frame(pR2(fit1))
set.seed(10)
fit2<-lm(predicted~observed)
r2 <- format(summary(fit2)$adj.r.squared, digits=3)
#plot(predicted,observed)

#capture.output(paste(i, "rsq",rsq, "r2",r2, "pr2",pr2[4,], "pr2",pr2[6,]), file="refcal.lm.out.txt", append=T) 
capture.output(paste(i, "rsq",rsq), file="refcal.rsq.txt", append=T) 

}}
##### END LOOP #######
######################### convert txt to csv 