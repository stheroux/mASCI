# screen metrics

#################################
# Import data ------------------
#################################

combined<-read.csv("combined.not.scaled.csv", row.names=1, header=T)

# Import stations data 
row.names(combined) == row.names(stations1)
#row.names(combined) == row.names(stations.num)

z<-length(row.names(combined))
foo<-which (colSums(!is.na(combined) == 0) > (100)) # find columns with 2000+ NAs 
if (length(foo) > 0 ) {combined<-combined[,-foo]}

metrics.list<-colnames(combined)

combined.str<-subset(combined, row.names(combined) %in% row.names(stations.str))
combined.rc<-subset(combined, row.names(combined) %in% row.names(stations.rc))

##########################################################
# Performance of individual metrics -----------------------
##########################################################
# ANOVA across PSA6c regions for ref stations ----------
# < 2 (or 3? )
Pv.list<-list()
Fv.list<-list()
metrics1<-list()

  for (i in metrics.list) {
    fit<-aov(combined.rc[[i]]~stations.rc$psa6c)
    Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
    Fv<-summary(fit)[[1]][["F value"]][[1]]
    Pv.list[[i]]<-Pv
    Fv.list[[i]]<-Fv
  } 

anova.psa<-cbind(Pv.list, Fv.list)
#write.csv(anova.psa, paste0(assemblage,".anova.PSA6c.ref.csv"))

i<-"cnt.spp.BCG2_raw"
ggplot(combined.rc, aes(x=stations.rc$psa6c, y=combined.rc[,i])) + 
  geom_boxplot() +
  geom_jitter(cex=0.4) + ylab(label = i) + xlab("PSA6C")


# ANOVA across ref/int/str ------------------

Pv.list<-list()
Fv.list<-list()

for (i in metrics.list) {
  foo<-as.matrix(combined[,i])
  fit<-aov(foo~stations1$SiteStatus)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  Pv.list[[i]]<-Pv
  Fv.list[[i]]<-Fv
}

anova.refintstr<-cbind(Pv.list, Fv.list)
#write.csv(anova.refintstr, file=paste0(assemblage,".anova.ref.stress.csv"))

# ttest across ref/str ------------------------------

tt<-list()
ttp<-list()

for (i in metrics.list) {
tt.out<-t.test(combined.rc[,i], combined.str[,i])
tt[i]<-tt.out$statistic
ttp[i]<-tt.out$p.value
}

tt.out2<-cbind(tt,ttp)

########################################################################################
# Frequence of one's and zero's in the ref dataset -------------------------
########################################################################################
# MaxFreqZero<-.33
# MaxCount<-5

library(plyr)

FreqZeroRefCal<-list()
for (i in metrics.list) { # i<-"prop.spp.Salinity.BF"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==0, na.rm = T)
  if (is.na(foo)) { FreqZeroRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqZeroRefCal[i]<-(foo/length1)*100 }
}

FreqOneRefCal<-list()
for (i in metrics.list) { # i<-"Salinity.BF.richness"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==1, na.rm = T)
  if (is.na(foo)) { FreqOneRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqOneRefCal[i]<-(foo/length1)*100 }
}

Freqs<-cbind(FreqZeroRefCal, FreqOneRefCal)

# Rafi SN --------------------------------------------------------------
# F value less than 3
# variance across all sites, variance at repeat samplings

DF<-combined
DF$SeqID<-row.names(DF)
StationCode.mean<-ddply(DF, .(SeqID), colwise(mean))
StationCode.sd<-ddply(DF, .(SeqID), colwise(sd))
StationCode.var<-ddply(DF, .(SeqID), colwise(var))
SN.rafi<-cbind(StationCode.mean, StationCode.sd, StationCode.var)
#write.csv(SN.rafi, paste0(assemblage,".Rafi.SN.csv"))

# Rafi ANOVA, variance within station code, want F < 3 -------------------------------
Pv.list<-list()
Fv.list<-list()

refz<-droplevels(subset(stations1, stations1$SiteStatus=="Reference"))
refz.stations.codes<-unique(refz$SeqID)
combined.rcz<-subset(combined, rownames(combined) %in% refz.stations.codes)
combined.rcz$SeqID<-row.names(combined.rcz)
combined.rcz<-merge(combined.rcz, stations1[,c("stationcode.correct", "SeqID")], by="SeqID", sort=F, all.x=T, all.y=F)
                        
for (i in metrics.list) { 
  fit<-aov(combined.rcz[[i]]~combined.rcz$stationcode.correct)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  Pv.list[[i]]<-Pv
  Fv.list[[i]]<-Fv
    }
rafi.anova.out<-cbind(Pv.list, Fv.list)
#write.csv(rafi.anova.out, file=paste0(assemblage,".anova.within.stationcode.csv"))


# Range test (Stevenson 2013) --------------------------------------------
# the median value of a candidate metric was > 0 in either reference or disturbed sites
Range.ref<-list()
Range.str<-list()
for (i in metrics.list) { 
  a <- median(combined.rc[[i]], na.rm=T) 
  b <- median(combined.str[[i]], na.rm=T) 
  Range.ref[i] <- a
  Range.str[i] <- b 
  }
Range<-cbind(Range.ref, Range.str)
#write.csv(Range, file = paste0(assemblage,".Range.Stevenson.csv"))


# signal to noise (Stoddard) ----------------------------------------
# variance across all sites /  variance at repeat samplings 
# > 2 is best , periphyton 1 or 1.5 okay 
SN<-list()
for (i in metrics.list) { 
SN[i] <- var(combined[[i]], na.rm=T) / mean(StationCode.var[[i]], na.rm=T) }
#write.csv(SN, file=paste0(assemblage,".SN.stoddard.csv"))

# Now calculate lm to disturbance gradients ----------------------
# from Cao et al 2007
#list3<- c("DayOfYear", "Elevation", "XSLOPE", "NHDSLOPE", "New_Lat", "New_Long", "LST32AVE", "MEANP_WS", "PPT_00_09", "CondQR50", "TMAX_WS", "XWD_WS", 
#  "KFCT_AVE", "BDH_AVE", "PRMH_AVE" )
list3<-c("maxw1_hall", "Ag_2011_5K", "URBAN_2011_5K", "CODE_21_2011_5K","RoadDens_5K",
               "PAVED_INT_1K")

list3<-tolower(list3)
r2.list<-list()
slope.list<-list()

colnames(stations1)<-tolower(colnames(stations1))

for (x in metrics.list) { 
  for (y in list3) {  #or list2 for all vars
    NAsum<-sum(is.na(stations1[[y]]))
    z<-length(row.names(stations1))
    if((NAsum/z) < .20) {  
      fit<-lm(stations1[[y]]~combined[[x]], na.action = "na.exclude")
      r2 <- format(summary(fit)$adj.r.squared, digits=3) 
      slope <-format(coef(fit)[2], digits=3)
      saveas<-paste(x,y)
      r2.list[[saveas]]<-r2
      slope.list[[saveas]]<-slope
    }}}

lm.out<-as.matrix(cbind(r2.list, slope.list))
#write.csv(lm.out, "lm.out.csv")

if (diatoms==T) {write.csv(lm.out, "diatoms.lm.out.gradients.csv")}
if (sba==T) {write.csv(lm.out, "sba.lm.out.gradients.csv")}
if (hybrid==T) {write.csv(lm.out, "hybrid.lm.out.gradients.csv")}


# SD
SD<-list()
#cal<-droplevels(subset(calval, calval$CalVal=="Cal"))
#combined.rc<-droplevels(subset(combined, row.names(combined) %in% row.names(cal)))

for (i in metrics.list) { 
  SD[i] <- sd(combined.rc[[i]], na.rm=T) }


# combine results ------------------------------------------------------------
results<-cbind(anova.psa,
               anova.refintstr,
               tt.out2,
               FreqZeroRefCal,
               FreqOneRefCal,
               #freq.max,
               #freq.mode,
               #SN.rafi,
               rafi.anova.out,
               Range,
               SN, 
               SD)
results<-as.data.frame(results)

colnames(results)<-c("anova.psa.pv", "anova.psa.fv",
                     "anova.refintstr.pv","anova.refintstr.fv",
                     "tt.out2.t", "tt.out2.p",
                     "FreqZero", "FreqOne",
                     "anova.rafi.pv", "anova.rafi.fv",
                     "Range.ref", "Range.str",
                     "SN", "SD")
results<-as.matrix(results)
write.csv(results, file=paste0(assemblage,".results.csv"))
results<-as.data.frame(results)

results2<-results
##################################
# add rsq info --------------
##################################

results2$metric<-row.names(results2)
rsq<-read.csv(file = paste0(assemblage, ".refcal.rsq.csv"), header=T)
rsq$metric<-rsq$top.metrics
results3<-merge(results2, rsq, by="metric", sort=F, all.x=T, all.y=F)

del1<-which(results3$rsq<0.2)
results3<-results3[-del1,]
write.csv(as.matrix(results3), file=paste0(assemblage, ".results3.csv") )



