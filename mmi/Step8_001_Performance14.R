
# compare ASCI performance 

#################################
# prep ------------------------
#################################

# combine O/E and MMI scores 

oe.asv<-read.csv("OEout/diatoms.asv.best/Scores.all.csv")
oe.asv<-as.data.frame(oe.asv[c(1,4)])
names(oe.asv)<-c("SeqID","OE.asv")

oe<-read.csv("OEout/diatoms.genus.best/Scores.all.csv")
oe<-oe[,c(1,4)]
names(oe)<-c("SeqID", "OE")

mmi.asv<-read.csv("R/mmi/diatoms/diatoms.asv/diatoms.combined.win.metrics.csv")
mmi.asv<-mmi.asv[,c("X", "MMI_scaled")]
names(mmi.asv)<-c("SeqID", "MMI.asv")

mmi<-read.csv("R/mmi/diatoms/diatoms/diatoms.combined.win.metrics.csv")
mmi<-mmi[,c("X", "MMI_scaled")]
names(mmi)<-c("SeqID", "MMI")

scores<-merge(oe.asv,oe, by="SeqID", all.x=T, all.y=T, sort=F)
scores<-merge(scores,mmi.asv, by="SeqID", all.x=T, all.y=T, sort=F)
scores<-merge(scores,mmi, by="SeqID", all.x=T, all.y=T, sort=F)
scores<-merge(scores, env,by="SeqID", all.x=T, all.y=F, sort=F )

ind<-c("OE.asv","OE","MMI.asv", "MMI")

scores$new_lat<-scores$Lat
scores$new_long<-scores$Long

scores.cal<-subset(scores, scores$CalVal %in% "Cal")
scores.cal<-droplevels(scores.cal)
scores.val<-subset(scores, scores$CalVal %in% "Val")
scores.val<-droplevels(scores.val)

scores.str.cal<-subset(scores, scores$StrCalVal %in% "Cal")
scores.str.val<-subset(scores, scores$StrCalVal %in% "Val")

scores.str<-rbind(scores.str.cal, scores.str.val)

##############################################
# Means at ref sites -------------
##############################################
mean.out<-list()

for (i in ind) {     #   i<-"OE"
  mean.out[[paste(i, "RefCal")]]<-mean(scores.cal[,i], na.rm=T) 
  mean.out[[paste(i, "RefVal")]]<-mean(scores.val[,i], na.rm=T)
  }

mean.out<-as.data.frame(mean.out)
mean.out<-as.data.frame(t(mean.out))
write.csv(mean.out, "means.out.csv")

##############################################
# ANOVA PSA region  -------------
##############################################
anova.out<-list()

for (i in ind) {     #   i<-"OoverE.d"
  foo<-aov(scores.cal[,i]~scores.cal$psa6c)
  anova.out[[paste(i, "RefCal")]]<-summary(foo)[[1]][["F value"]][[1]]
  foo<-aov(scores.val[,i]~scores.val$psa6c)
  anova.out[[paste(i, "RefVal")]]<-summary(foo)[[1]][["F value"]][[1]]
}

anova.out<-as.matrix(t(anova.out))
anova.out<-as.matrix(t(anova.out))
write.csv(anova.out, "anova.psa.out.csv")


##############################################
# Variance explained by natural and str gradients   -------------
##############################################

library(randomForest)
source("bin/cand.var.bu.R")
envgrads<-tolower(cand.var)

# env gradients at calibration sites ------------
out1<-list()

# cal 
for (i in ind) { # i<-"MMI.d"
  sub1<-scores.cal
  foo<-which(is.na(sub1[[i]]))
  if (length(foo) > 0) {sub1<-sub1[-foo,]}
  
  set.seed(13)
  rf.out<-randomForest(y=sub1[[i]], x=sub1[,envgrads], 
                       ntree=500,importance=TRUE, norm.votes=TRUE, 
                       keep.forest=TRUE,na.action=na.omit, proximity=T)
  rsq<-as.data.frame(rf.out$rsq) ; rsq<-(rsq[500,]) 
  out1[paste(i, "RefCal")]<-rsq
}

# val 
for (i in ind) { # i<-"MMI.d"
  sub1<-scores.val
  foo<-which(is.na(sub1[[i]]))
  if (length(foo) > 0) {sub1<-sub1[-foo,]}
  
  set.seed(13)
  rf.out<-randomForest(y=sub1[[i]], x=sub1[,envgrads], 
                       ntree=500,importance=TRUE, norm.votes=TRUE, 
                       keep.forest=TRUE,na.action=na.omit, proximity=T)
  rsq<-as.data.frame(rf.out$rsq) ; rsq<-(rsq[500,]) 
  out1[paste(i, "RefVal")]<-rsq
}

# str gradients ------------

str.grads<-c(#"Ag_2006_1K",
  #"Ag_2006_5K",
  #"Ag_2006_WS",
  "Ag_2011_1K",
  "Ag_2011_5K",
  "Ag_2011_WS",
  #"AtmSO4",
  #"CODE_21_2006_1K",
  #"CODE_21_2006_5K",
  #"CODE_21_2006_WS",
  "CODE_21_2011_1K",
  "CODE_21_2011_5K",
  "CODE_21_2011_WS",
  #"CondQR01",
  #"CondQR05",
  #"CondQR25",
  #"CondQR75",
  #"CondQR95",
  #"CondQR99",
  #"CondQR99c",
  #"GravelMineDensL_R5K",
  "InvDamDist",
  #"Log_N_MEAN",
  #"Log_P_MEAN",
  "MINES_5K",
  #"N_MEAN",
  #"NHD_SO",
  #"P_MEAN",
  "PAVED_INT_1K",
  "PAVED_INT_5K",
  "PAVED_INT_WS",
  #"PerManMade_WS",
  #"Replicate",
  "RoadDens_1K",
  "RoadDens_5K",
  #"RoadDens_corrected"
  "RoadDens_WS",
  #"SumAve_P",
  #"URBAN_2006_1K",
  #"URBAN_2006_5K",
  #"URBAN_2006_WS",
  "URBAN_2011_1K",
  "URBAN_2011_5K",
  "URBAN_2011_WS",
  NULL)

str.grads<-tolower(str.grads)

out2<-list()

scores.int<-subset(scores, scores$SiteStatus %in% "Intermediate") 
scores.cal2<-rbind(scores.cal, scores.int, scores.str.cal)
scores.val2<-rbind(scores.val, scores.int, scores.str.val)

# cal 
for (i in ind) { # i<-"MMI.d"
  sub1<-scores.cal2
  foo<-which(is.na(sub1[[i]]))
  if (length(foo) > 0) {sub1<-sub1[-foo,]}
  
  foo<-complete.cases(sub1[,str.grads])
  sub1<-sub1[foo,]
  
  set.seed(11)
  rf.out<-randomForest(y=sub1[[i]], x=sub1[,str.grads], 
                       ntree=500,importance=TRUE, norm.votes=TRUE, 
                       keep.forest=TRUE,na.action=na.omit, proximity=T)
  rsq<-as.data.frame(rf.out$rsq) ; rsq<-(rsq[500,]) 
  out2[paste(i, "RefCal")]<-rsq
}

# val 
for (i in ind) { # i<-"MMI.d"
  sub1<-scores.val2
  foo<-which(is.na(sub1[[i]]))
  if (length(foo) > 0) {sub1<-sub1[-foo,]}
  
  foo<-complete.cases(sub1[,str.grads])
  sub1<-sub1[foo,]
  
  set.seed(11)
  rf.out<-randomForest(y=sub1[[i]], x=sub1[,str.grads], 
                       ntree=500,importance=TRUE, norm.votes=TRUE, 
                       keep.forest=TRUE,na.action=na.omit, proximity=T)
  rsq<-as.data.frame(rf.out$rsq) ; rsq<-(rsq[500,]) 
  out2[paste(i, "RefVal")]<-rsq
}


out3<-cbind(out1, out2)
out3<-as.data.frame(out3)
names(out3)<-c("envgrads", "strgrads")
write.csv(as.matrix(out3), "envgrads.strgrads.csv")


##############################################
# standard deviations ref sites -------------
##############################################
sd.out<-list()

for (i in ind) {     #   i<-"OoverE.d"
  sd.out[[paste(i, "RefCal")]]<-sd(scores.cal[,i], na.rm=T) 
  sd.out[[paste(i, "RefVal")]]<-sd(scores.val[,i], na.rm=T)
}

sd.out<-as.data.frame(sd.out)
sd.out<-as.data.frame(t(sd.out))
write.csv(sd.out, "sd.out.csv")


# ##############################################
# # get repeat stations  -------------
# ##############################################
#library(plyr)

foo<-plyr::count(scores$stationcode.correct)
foo<-subset(foo, foo$freq>1)
foo<-droplevels(foo)
repeat.stations.list<-foo$x

repeat.stations<-subset(scores, scores$stationcode.correct %in% repeat.stations.list)
repeat.stations.cal<-subset(repeat.stations, CalVal=="Cal")
repeat.stations.val<-subset(repeat.stations, CalVal=="Val")


##############################################
# sd Station Code within sites   -------------
##############################################
sd.out<-list()

for (i in ind) {     #   i<-"MMI.asv"
  for (j in (unique(repeat.stations.cal$stationcode.correct))) { # j<-"312NPPCBR"
    foo1<-subset(scores, scores$stationcode.correct %in% j)
    foo1<-droplevels(foo1)
    sd.out[[paste(i,j, "Cal")]]<-sd(foo1[[i]], na.rm = T)
  }}

for (i in ind) {     #   i<-"MMI.d"
  for (j in (unique(repeat.stations.val$stationcode.correct))) { # j<-"312NPPCBR"
    foo1<-subset(scores, scores$stationcode.correct %in% j)
    foo1<-droplevels(foo1)
    sd.out[[paste(i,j, "Val")]]<-sd(foo1[[i]], na.rm = T)
  }}

sd.out<-as.data.frame(t(t(sd.out)))
foo<-as.data.frame(do.call(rbind, strsplit(as.character(row.names(sd.out)), " ")))
sd.out<-cbind(sd.out, foo)

foo<-which(is.na(sd.out$V1))
sd.out<-as.data.frame(sd.out[-foo,])
names(sd.out)<-c("value", "ind", "stationcode", "CalVal")

sd.out2<-list()

for (i in unique(sd.out$ind)) { # i<-"MMI"
  for (j in c("Cal", "Val")) { # j<-"Cal"
  foo1<-subset(sd.out, ind %in% i) 
  foo1<-subset(sd.out, CalVal %in% j) 
  mn<-mean(as.numeric(foo1$value), na.rm = T)
  sd.out2[paste(i, j)]<-mn
}}

sd.out2[paste("OE", "Cal")] <- NA
sd.out2[paste("OE", "Val")] <- NA

write.csv(as.matrix(sd.out2), "sd.stationcode.csv")

########################
# ttest ref/str -------
########################

ttest.out<-list()

for (i in ind) {     #   i<-"OE"
  foo<-t.test(scores.cal[,i],scores.str.cal[,i])
  ttest.out[paste(i, "RefCal")]<-foo$statistic
  foo<-t.test(scores.val[,i],scores.str.val[,i])
  ttest.out[paste(i, "RefVal")]<-foo$statistic
  
}

ttest.out<-as.matrix(t(ttest.out))
ttest.out<-as.matrix(t(ttest.out))
write.csv(ttest.out, "ttest.out.csv")


##############################################
# Levene's test  -------------
##############################################

# library(car)
# levenetest.out<-list()
# 
# for (i in ind) {     #   i<-"OoverE.d"
#   foo<-leveneTest(scores.cal[,i]~scores.cal$psa6c)
#   levenetest.out[[paste(i, "RefCal")]]<-foo$`F value`[1]
#   foo<-leveneTest(scores.val[,i]~scores.val$psa6c)
#   levenetest.out[[paste(i, "RefVal")]]<-foo$`F value`[1]
# }
# levenetest.out<-as.matrix(t(levenetest.out))
# levenetest.out<-as.matrix(t(levenetest.out))
# write.csv(levenetest.out, "levenetest.out.csv")


######################
# Percentiles -----------
######################
first.perc.ref<-list()
tenth.perc.ref<-list()
thirt.perc.ref<-list()
perc.rc.above.tenth.perc<-list()
perc.rv.above.tenth.perc<-list()
perc.str.below.tenth<-list()

for (i in ind) { # i<-"MMI.d" 
  foo1<-qnorm(p = 0.1, mean = mean(scores.cal[,i], na.rm=T), sd = sd(scores.cal[,i], na.rm=T))
  tenth.perc.ref[paste(i)]<-foo1
  perc.rc.above.tenth.perc[paste(i)] <- (length(which(scores.cal[,i] > foo1)) / length(which(scores.cal[,i] > 0))) * 100
  perc.rv.above.tenth.perc[paste(i)] <- (length(which(scores.val[,i] > foo1)) / length(which(scores.val[,i] > 0))) * 100
  perc.str.below.tenth[paste(i)] <- (length(which(scores.str[,i] < foo1)) / length(which(scores.str[,i] > 0))) * 100
  
  foo1<-qnorm(p = 0.01, mean = mean(scores.cal[,i], na.rm=T), sd = sd(scores.cal[,i], na.rm=T))
  first.perc.ref[paste(i)]<-foo1
  
  foo1<-qnorm(p = 0.3, mean = mean(scores.cal[,i], na.rm=T), sd = sd(scores.cal[,i], na.rm=T))
  thirt.perc.ref[paste(i)]<-foo1
  
}

qnorm(p = 0.01, mean = mean(scores.cal$MMI.asv, na.rm=T), sd = sd(scores.cal$MMI.asv, na.rm=T))
qnorm(p = 0.1, mean = mean(scores.cal$MMI.asv, na.rm=T), sd = sd(scores.cal$MMI.asv, na.rm=T))
qnorm(p = 0.3, mean = mean(scores.cal$MMI.asv, na.rm=T), sd = sd(scores.cal$MMI.asv, na.rm=T))


out1<-cbind(first.perc.ref, tenth.perc.ref,thirt.perc.ref, 
            perc.rc.above.tenth.perc, perc.rv.above.tenth.perc, perc.str.below.tenth)
out1<-as.data.frame(out1)

out1.perc<-out1

write.csv(as.matrix(out1), "percentiles.csv")


#############################################
# percent classified "correctly" ------------
#############################################

perccor<-list()

for (i in ind) {     #   i<-"MMI.asv"
  # total length of refcal
  len<-length(scores.cal$SeqID)

  # get threshold
  sub1<-subset(out1.perc, row.names(out1.perc) %in% i)
  thr<-sub1$tenth.perc.ref[[1]]

  # get refcal for index
  sub2<-subset(scores.cal, scores.cal[i] > thr ) # refcal

  # length of refcal > thr
  len2<-length(sub2$SeqID)
  perc1<-len2/len*100
  perccor[i]<-perc1
}


write.csv(as.matrix(perccor), "percentiles.above.thr.csv")

####################################################
#### Spearman ref/int/str ####### -----------
####################################################

# fit<-lm(otu[[i]]~geochem$Sal)
# r2 <- format(summary(fit)$adj.r.squared, digits=3)
r2.out<-list()
lms<-c("TN", "TP", "SpCond" )

for (i in ind) { # i<-"MMI.hybrid"
  for (j in lms) { # j<-"Nitrogen_Total_mgPerL"
    
    sub1<-scores.cal2
    foo<-which(is.na(sub1[[i]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    foo<-which(is.na(sub1[[j]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    
    # spearman 
    fit<- cor.test(y=as.numeric(sub1[[i]]), x=as.numeric(sub1[[j]]), method = 'spearman', exact = F)
    r2 <- format((fit$estimate[1]), digits=3)
    r2.out[paste(i, j, "Cal")] <- r2
    
    sub1<-scores.val2
    foo<-which(is.na(sub1[[i]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    foo<-which(is.na(sub1[[j]]))
    if (length(foo) > 0) {sub2<-sub2[-foo,]}
    
    # spearman 
    fit<- cor.test(y=as.numeric(sub1[[i]]), x=as.numeric(sub1[[j]]), method = 'spearman', exact = F)
    r2 <- format((fit$estimate[1]), digits=3)
    r2.out[paste(i, j, "Val")] <- r2
    
  }}
    
spearman.out<-as.matrix(t(t(r2.out)))
write.csv(as.matrix(spearman.out), "spearman.out.csv")

# expanded 
r2.out<-list()
lms<-c(#"Chloride_mgPerL", 
       "DOC",
       "TN",
       #"PCT_SAFN",
       "pH", 
       "TP",
       #"Phosphorus_as_P_mgPerL", 
       "RoadDens_5K",
       "SpCond",
       "TEMP_00_09",
       "URBAN_2011_5K",
       "MaxW1_HALL",
       #"W1_HALL_SWAMP",
       #"DayOfYear",
       #"Month", 
       "AREA_SQKM", 
       "AtmCa", 
       "AtmMg",
       "BDH_AVE", 
       "CaO_Mean", 
       "CondQR50", 
       "ELEV_RANGE", 
       "KFCT_AVE",
       #"LogWSA",
       "LPREM_mean",
       "LST32AVE",
       "MAX_ELEV",
       "MAXWD_WS",
       "MEANP_WS",
       "MgO_Mean",
       "MINP_WS",
       "New_Lat",
       "New_Long",
       "PCT_CENOZ",
       "PCT_NOSED",
       "PCT_QUART",
       "PCT_SEDIM",
       "PCT_VOLCNC",
       "PPT_00_09",
       "PRMH_AVE",
       #"PSA6c",
       "S_Mean",
       "SITE_ELEV",
       "TEMP_00_09",
       "TMAX_WS",
       "UCS_Mean",
       "XWD_WS") 

lms<-tolower(lms)
scores$ph<-scores$pH
scores$doc<-scores$DOC
scores$tn<-scores$TN
scores$tp<-scores$TP
scores$spcond<-scores$SpCond
scores$maxw1_hall<-scores$MaxW1_HALL

scores.cal2$ph<-scores.cal2$pH
scores.cal2$doc<-scores.cal2$DOC
scores.cal2$tn<-scores.cal2$TN
scores.cal2$tp<-scores.cal2$TP
scores.cal2$spcond<-scores.cal2$SpCond
scores.cal2$maxw1_hall<-scores.cal2$MaxW1_HALL

scores.val2$ph<-scores.val2$pH
scores.val2$doc<-scores.val2$DOC
scores.val2$tn<-scores.val2$TN
scores.val2$tp<-scores.val2$TP
scores.val2$spcond<-scores.val2$SpCond
scores.val2$maxw1_hall<-scores.val2$MaxW1_HALL


for (i in ind) { # i<-"MMI.hybrid"
  for (j in lms) { # j<-"Nitrogen_Total_mgPerL"
    
    sub1<-scores.cal2
    foo<-which(is.na(sub1[[i]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    foo<-which(is.na(sub1[[j]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    
    # spearman 
    fit<- cor.test(y=as.numeric(sub1[[i]]), x=as.numeric(sub1[[j]]), method = 'spearman', exact = F)
    r2 <- format((fit$estimate[1]), digits=3)
    r2.out[paste(i, j, "Cal")] <- r2
    
    sub1<-scores.val2
    foo<-which(is.na(sub1[[i]]))
    if (length(foo) > 0) {sub1<-sub1[-foo,]}
    foo<-which(is.na(sub1[[j]]))
    if (length(foo) > 0) {sub2<-sub2[-foo,]}
    
    # spearman 
    fit<- cor.test(y=as.numeric(sub1[[i]]), x=as.numeric(sub1[[j]]), method = 'spearman', exact = F)
    r2 <- format((fit$estimate[1]), digits=3)
    r2.out[paste(i, j, "Val")] <- r2
    
  }}


spearman.out<-as.matrix(t(t(r2.out)))
write.csv(as.matrix(spearman.out), "spearman.expanded.out.csv")

########################################################
### make performance table ----------------------------
########################################################
v1=read.csv('means.out.csv')
names(v1)<-c("Ind","Means")
v1<-v1[order(v1$Ind),]

v2=read.csv('anova.psa.out.csv')
names(v2)<-c("Ind","anvPSA")
v2<-v2[order(v2$Ind),]

v3=read.csv('envgrads.strgrads.csv')
names(v3)<-c("Ind","env.grads", "str.grads")
v3<-v3[order(v3$Ind),]

v4=read.csv('sd.out.csv')
names(v4)<-c("Ind","SD1")
v4<-v4[order(v4$Ind),]

v5=read.csv('sd.stationcode.csv')
names(v5)<-c("Ind","SD2")
v5<-v5[order(v5$Ind),]

v6=read.csv('ttest.out.csv')
v6<-as.data.frame(v6)
names(v6)<-c("Ind","ttest")
v6<-v6[order(v6$Ind),]

v7=read.csv('envgrads.strgrads.csv')
names(v7)<-c("Ind","env.grads", "str.grads")
v7<-v7[order(v7$Ind),]

v8=read.csv('spearman.out.csv')
names(v8)<-c("Index","Rho")
foo<-do.call(rbind, strsplit(as.character(v8$Index), " "))
foo<-as.data.frame(foo)
names(foo)<-c("Ind", "Stressor")
v8<-cbind(v8, foo)
v8<-v8[order(v8$Ind),]

v9=read.csv('spearman.expanded.out.csv')
names(v9)<-c("Index","Rho")
foo<-do.call(rbind, strsplit(as.character(v9$Ind), " "))
foo<-as.data.frame(foo)
names(foo)<-c("Ind", "Stressor")
v9<-cbind(v9, foo)
v9<-v9[order(v9$Ind),]

v.all<-as.data.frame(cbind(v1$Means, v2$anvPSA, v3$env.grads, v4$SD1,v5$SD2, v6$ttest, v7$str.grads))
row.names(v.all)<-v1$Ind
names(v.all)<-c("means", "PSA", "Var1", "SD1","SD2", "ttest", "Var2")

tn<-as.data.frame(subset(v8, v8$Stressor=="TN"))
tp<-as.data.frame(subset(v8, v8$Stressor=="TP"))
sc<-as.data.frame(subset(v8, v8$Stressor=="SpCond"))

v.all$TN<-tn$Rho
v.all$TP<-tp$Rho
v.all$SpCond<-sc$Rho

write.csv(v.all, "table.out.csv")

# 
# #####################################
# # get randomforest importance -------
# 
# # OE 
# 
# setwd("~/Documents/R/ASCI/PERFORMANCE_newnew/")
# library(randomForest)
# library(reshape2)
# 
# rf1<-load("~/Documents/R/ASCI/OE/diatoms.best/diatom.RF.OE.Rdata")
# foo<-importance(diatom.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/diatoms.rf.oe.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE/sba.best/sba.RF.OE.Rdata")
# foo<-importance(sba.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/sba.rf.oe.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE/hybrid.best/hybrid.RF.OE.Rdata")
# foo<-importance(hybrid.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/hybrid.rf.oe.csv")
# 
# # ecoregion -----
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SN/diatoms.best/diatom.RF.OE.Rdata")
# foo<-importance(diatom.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/diatoms.rf.oe.sn.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SN/sba.best/sba.RF.OE.Rdata")
# foo<-importance(sba.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/sba.rf.oe.sn.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SN/hybrid.best/hybrid.RF.OE.Rdata")
# foo<-importance(hybrid.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/hybrid.rf.oe.sn.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SC/diatoms.best/diatom.RF.OE.Rdata")
# foo<-importance(diatom.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/diatoms.rf.oe.sc.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SC/sba.best/sba.RF.OE.Rdata")
# foo<-importance(sba.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/sba.rf.oe.sc.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/SC/hybrid.best/hybrid.RF.OE.Rdata")
# foo<-importance(hybrid.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/hybrid.rf.oe.sc.csv")
# 
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/CH/diatoms.best/diatom.RF.OE.Rdata")
# foo<-importance(diatom.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/diatoms.rf.oe.ch.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/CH/sba.best/sba.RF.OE.Rdata")
# foo<-importance(sba.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/sba.rf.oe.ch.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/CH/hybrid.best/hybrid.RF.OE.Rdata")
# foo<-importance(hybrid.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/hybrid.rf.oe.ch.csv")
# 
# 
# 
# # extravars
# rf1<-load("~/Documents/R/ASCI/OE_PSA/extravars/diatoms.best/diatom.RF.OE.Rdata")
# foo<-importance(diatom.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/diatoms.rf.oe.ev.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/extravars/sba.best/sba.RF.OE.Rdata")
# foo<-importance(sba.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/sba.rf.oe.ev.csv")
# 
# rf1<-load("~/Documents/R/ASCI/OE_PSA/extravars/hybrid.best/hybrid.RF.OE.Rdata")
# foo<-importance(hybrid.rf.oe)
# write.csv(as.matrix(foo), "randomforest_imps/hybrid.rf.oe.ev.csv")
# 
# 
# # MMI --------------------------------------------------------
# # diatoms
# win<-read.csv("~/Documents/R/ASCI/pMMI_newnew/diatoms/win.metrics.csv", row.names=1, header=T)
# foo<-grep("*_raw", win$x)
# if (length(foo) >0) {win<-win[-foo,]}
# win<-droplevels(win)
# win
# for (i in win$x) {  # i<-"Salinity.BF.richness"
#   load(paste0("~/Documents/R/ASCI/pMMI_newnew/diatoms/diatoms.rf.models/diatoms.",i,"RF.model.Rdata"))
#   print(i)
#   foo<-as.data.frame(importance(rf.out,type=1))
#   out2<-paste(i, as.vector(row.names(foo)), as.vector(foo$`%IncMSE`))   
#   capture.output(out2, file = "randomforest_imps/diatoms.mmi.rf.imp.csv", append = T)
# 
# }
# 
# # sba
# 
# # hybrid
# win<-read.csv("~/Documents/R/ASCI/pMMI_newnew/hybrid/win.metrics.csv", row.names=1, header=T)
# foo<-grep("*_raw", win$x)
# if (length(foo) >0) {win<-win[-foo,]}
# win<-droplevels(win)
# win
# for (i in win) {  # i<-"prop.spp.Trophic.E"
#   load(paste0("~/Documents/R/ASCI/pMMI_newnew/hybrid/hybrid.rf.models/hybrid.",i,"RF.model.Rdata"))
#   print(i)
#   foo<-as.data.frame(importance(rf.out,type=1))
#   out2<-paste(i, as.vector(row.names(foo)), as.vector(foo$`%IncMSE`))   
#   capture.output(out2, file = "randomforest_imps/hybrid.mmi.rf.imp.csv", append = T)
# 
# }
# 
# # combine all -----
# a<-read.csv("randomforest_imps/diatoms.rf.oe.csv", row.names=1)
# b<-read.csv("randomforest_imps/sba.rf.oe.csv", row.names=1)
# c<-read.csv("randomforest_imps/hybrid.rf.oe.csv", row.names=1)
# 
# d<-read.csv("randomforest_imps/diatoms.rf.oe.sn.csv", row.names=1)
# e<-read.csv("randomforest_imps/sba.rf.oe.sn.csv", row.names=1)
# f<-read.csv("randomforest_imps/hybrid.rf.oe.sn.csv", row.names=1)
# 
# g<-read.csv("randomforest_imps/diatoms.rf.oe.sc.csv", row.names=1)
# h<-read.csv("randomforest_imps/sba.rf.oe.sc.csv", row.names=1)
# i<-read.csv("randomforest_imps/hybrid.rf.oe.sc.csv", row.names=1)
# 
# j<-read.csv("randomforest_imps/diatoms.rf.oe.ev.csv", row.names=1)
# k<-read.csv("randomforest_imps/sba.rf.oe.ev.csv", row.names=1)
# l<-read.csv("randomforest_imps/hybrid.rf.oe.ev.csv", row.names=1)
# 
# source("~/Documents/R/ASCI/OE_PSA/extravars/cand.var.bu.morevars.R")
# vars<-as.data.frame(cand.var)
# 
# a$cand.var<-row.names(a)
# b$cand.var<-row.names(b)
# c$cand.var<-row.names(c)
# d$cand.var<-row.names(d)
# e$cand.var<-row.names(e)
# f$cand.var<-row.names(f)
# g$cand.var<-row.names(g)
# h$cand.var<-row.names(h)
# i$cand.var<-row.names(i)
# j$cand.var<-row.names(j)
# k$cand.var<-row.names(k)
# l$cand.var<-row.names(l)
# 
# vars1<-merge(vars, a[,c("cand.var","MeanDecreaseGini")], all.x = T, all.y = T,by = "cand.var")
# vars1<-merge(vars1, b[,c("cand.var","MeanDecreaseGini")], all.x = T, all.y = T,by = "cand.var")
# vars1<-merge(vars1, c[,c("cand.var","MeanDecreaseGini")], all.x = T, all.y = T,by = "cand.var")
# vars1<-merge(vars1, d[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, e[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, f[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, g[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, h[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, i[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, j[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, k[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# vars1<-merge(vars1, l[,c("cand.var","MeanDecreaseGini")], all.x = T,by = "cand.var")
# names(vars1)<-c("cand.var", "d.oe", "sba.oe", "h.oe", 
#                 "d.oe.sn", "sba.oe.sn", "h.oe.sn",
#                 "d.oe.sc", "sba.oe.sc", "h.oe.sc",
#                 "d.oe.ev", "sba.oe.ev", "h.oe.ev")
# 
# # 
# comb<-read.csv("randomforest_imps/mmis_predvars.csv")
# 
# comb.m<-acast(comb, comb$Predvars~comb$metric, 
#               value.var="MSE", 
#               fun.aggregate=sum) 
# comb.m<-as.data.frame(comb.m)
# comb.m$cand.var<-row.names(comb.m)
# 
# vars1<-merge(vars1, comb.m, all.x = T)
# 
# 
# write.csv(vars1, "randomforest_imps/all.imp.csv")

####################### 
# misc O/E stats -------------------

library(randomForest)

dat<-read.csv("OEout/diatom.best/Scores.all.csv", header=T, row.names=1)
fit<-lm(dat$OE.scores.O~dat$OE.scores.E)
coef(fit)["(Intercept)"]
coef(fit)["dat$OE.scores.E"]
ggplot(dat, aes(x=OE.scores.E, y=OE.scores.O)) + geom_point() + 
  geom_smooth(method="lm") + theme_bw() + 
  xlab(label = "Expected") + ylab(label = "Observed") +
  geom_abline(intercept=0, slope=1, color="red") + 
  ggsave("OEout/diatom.best/OEplot.pdf")



