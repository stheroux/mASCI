# learning lapply and making a screen.metrics function

### trying in parallel 

#lapply(1:3, function(x) c(x, x^2, x^3))


# test 
#win.metrics<-c("EpiRho.richness","prop.spp.Trophic.I", 
#               "OrgN.NALON.richness_raw", "OxyReq.DO_100.richness_raw")


screen.metrics<-function(i) {
  keep<-i
  combined.win.metrics<- data.frame(combined[,keep]) # should this be scored, not scored? originally combined.sc
  stations<-subset(stations, stations$SeqID %in% row.names(combined.win.metrics))
  stations$SeqID == row.names(combined.win.metrics)
  scores<-data.frame(stations$SiteStatus, rowMeans(combined.win.metrics, na.rm = F))
  colnames(scores)<-c("Type", "Means")
  fit<-aov(scores$Means~scores$Type)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  saveas<-paste(colnames(combined.win.metrics), collapse = ",")
  Pv.list[saveas]<-Pv
  Fv.list[saveas]<-Fv
  Fv
  
  # sd refcal
  scores.rc<-subset(scores, row.names(scores) %in% rc$SeqID)
  sd1<-sd(scores.rc$Means, na.rm = T)
  sd1.list[saveas]<-sd1
  sd1
  
  # PSA anova 
  #stations.rc<-subset(stations, stations$SeqID %in% row.names(scores.rc))
  fit<-aov(scores.rc$Means~rc$psa6c)
  Fv2<-summary(fit)[[1]][["F value"]][[1]]
  psa.anova.list[saveas]<-Fv2
  Fv2
 
  paste(Fv,",",sd1,",",Fv2)
  #return(paste("out",",",toString(i),",",Fv,",",sd1,",",Fv2))
  #capture.output(paste("out",",",toString(i),",",Fv,",",sd1,",",Fv2),file="out.txt", append=T)

  return(paste(Fv,",",sd1,",",Fv2))
  #capture.output(paste("out",",",as.character(i),",",Fv,",",sd1,",",Fv2),file="out.txt", append=T)
  
}



#sapply(foo, function(x) fun=screen.metrics(x))
#lapply(foo, function(x) fun=screen.metrics(x))
