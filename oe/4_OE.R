# Build an O/E model

#Set a working directory

#### Load packages 
library("sampling")
library("rms")
library("plyr")
library("ggplot2")
library("gridExtra")
library(cluster)
library(fpc)
library(caret)
source("bin/OE.load.and.source.R")
source("bin/OE.cand.vars.R")
source("bin/OE.caret.load.and.source.R")
source("bin/model.predict.RanFor.4.2_ed2.r") #overwrites earlier model.predict

dir.create("OEout")

diatoms = T  # pick 1
sba =    F #
combo =   F # 

#### KNOBS TO TURN 
prezabs = T
genuslev = T
plotz = T
exportfiles = T
remlowtaxa = F # default F

removeNA = T # deletes taxa with the name 'NA'
removeQual = T # remove qual sample
removeplanktonic = F # remove planktonic species

removeoutliersites=F # default F

if (diatoms==T) { 
  # predictornumz<-c(10) ; clusternumz<-c(4:12) ; Pcapz <- c(0.4,0.5) } 
 # predictornumz<-c(10) ; clusternumz<-c(11) ; Pcapz <- c(0.5) } # ASV
  predictornumz<-c(10) ; clusternumz<-c(5) ; Pcapz <- c(0.4) }  # by genus

if (sba==T) {
  #predictornumz<-c(10) ; clusternumz<-c(4:15) ; Pcapz <- c(0.4,0.5) }
  predictornumz<-c(10) ; clusternumz<-c(8) ; Pcapz <- c(0.4) }

if (combo==T) { 
  #predictornumz<-c(10) ; clusternumz<-c(4:15) ; Pcapz <- c(0.4,0.5) }
  predictornumz<-c(10) ; clusternumz<-c(15) ; Pcapz <- c(0.4) }


#export order 
#Scores.cal.rfe
#Scores.rv
#Scores.new.int
#Scores.new.str

#Import your bug data
if(diatoms==T) { 
tax<-otu.m.rbcL}
#tax.orig<-tax

#Import stations data, make sure #N/A are NA
stations<-env
source(file = "bin/ASCI_fonts.OE.R")
foo<-which(is.na(stations$SeqID))
stations<-stations[-foo,]
row.names(stations)<-stations$SeqID
stations.orig<-stations

stations<-subset(stations, psa6c!="NA", drop=T)
stations<-subset(stations, psa6c!="NV", drop=T)
stations<-droplevels(stations)
stations<-subset(stations, MostRecent=="y", drop=T)

###################################################################################
# ASSIGN CAL/VAL -------------------------------------------------
###################################################################################
stations.rc<-subset(stations,CalVal %in% "Cal", drop=T )
stations.rv<-subset(stations, CalVal=="Val",drop=T )
stations.r<-subset(stations, CalVal!="NA",)

###################################################################################
# Step 1 DATA PREP -------------------------------------------------
###################################################################################

# get taxonomy

foo<-which(tax$SeqID %in% stations$SeqID)
tax<-tax[foo,]

if (diatoms==T) { 
tax<-subset(tax, Phylum=="Bacillariophyta") }
if (sba==T) { 
  tax<-subset(tax, Phylum!="Bacillariophyta") }
tax<-droplevels(tax)

#remove low number taxa in each sample
if (remlowtaxa ==T) { tax<-subset(tax, tax$value>0.006) }
  
tax.sub<-tax
if (genuslev==T) {tax.m<-acast(tax.sub, SeqID~Genus, value.var="value", fun.aggregate=sum)}
if (genuslev==F) {tax.m<-acast(tax.sub, SeqID~OTUID, value.var="value", fun.aggregate=sum)}
tax.m.sp<-acast(tax.sub, SeqID~Species, value.var="value", fun.aggregate=sum)
write.csv(tax.m.sp, "OEout/species.csv")

#subset bug matrices to just the refcal sites
foo<-which(row.names(tax.m) %in% stations.rc$SeqID)
tax.rc<-as.data.frame(subset(tax.m, row.names(tax.m) %in% stations.rc$SeqID, drop = T))
foo<-which(colSums(tax.rc)==0) # no zero sum columns
if(length(foo)>0) {tax.rc<-tax.rc[,-foo]}

#subset stations.rc matrices to just the bug sites
stations.rc<-subset(stations.rc, row.names(stations.rc) %in% rownames(tax.rc))

#Check to make sure rows are aligned correctly
row.names(tax.rc)==row.names(stations.rc)

#If they aren't, you can do this:
tax.rc<-tax.rc[row.names(stations.rc),]
row.names(tax.rc)==row.names(stations.rc)

#Make a presence-absence version
if (prezabs==T) {tax.rc<-ifelse(tax.rc>0,1,0)}
if (prezabs==F) {tax.rc<-tax.rc}
tax.rc<-as.data.frame(tax.rc)

# export matrices 
if (diatoms==T) {write.csv(tax.m, "OEout/diatoms.tax.m.csv")}
if (sba==T) {write.csv(tax.m, "OEout/sba.tax.m.csv")}
if (combo==T) {write.csv(tax.m, "OEout/hybrid.tax.m.csv")}

tax.rc.norare<-tax.rc
tax.rc.norare<-ifelse(tax.rc.norare>0,1,0) #always cluster on a PA 

# lower threshold 
thresh<-(nrow(tax.rc)*0.025)
foo<-which(colSums(tax.rc.norare) < thresh)
if (length(foo)>0) {tax.rc.norare<-as.data.frame(tax.rc.norare[,-foo])}

# upper threshold 
thresh<-(nrow(tax.rc)*0.95)
foo<-which(colSums(tax.rc.norare) > thresh)
if (length(foo) > 0) { tax.rc.norare<-as.data.frame(tax.rc.norare[,-foo]) } 

foo<-which(rowSums(tax.rc.norare)==0)
if (length(foo)>0) {tax.rc.norare<-tax.rc.norare[-foo,]}

stations.rc<-stations.rc[which(row.names(stations.rc) %in% row.names(tax.rc.norare)),]
tax.rc<-tax.rc[which(row.names(tax.rc) %in% row.names(stations.rc)),]

foo<-paste("prezabs", prezabs, "genuslev", genuslev, "remlowtaxa", remlowtaxa, 
           "removeNA", removeNA, "removeQual", removeQual, "removeplanktonic", removeplanktonic)

write.csv(foo, "OEout/settings.txt")

#################################################
####### Start Loop ----- ***********
#################################################

for (x in predictornumz) { 
  for (y in clusternumz) { 
    for (p in Pcapz) {
      predictornum<-x
      clusternum<-y
      Pcap<-p
      
      ##################################################################;
      # CLUSTERING -----------------------------------
      ##################################################################;
      
      #Recommend doing clustering on data sets with rare taxa removed. 
      #But build models using the data set that includes rare taxa
      #https://www.r-statistics.com/2013/08/k-means-clustering-from-r-in-action/
      #fit.km <- kmeans(tax.rc,clusternum, nstart=25)     
      #stations.rc$BG<-as.factor(fit.km$cluster)
      
      set.seed(10)
      fit.km <- kmeans(tax.rc.norare,clusternum, nstart=25)     
      stations.rc$BG<-as.factor(fit.km$cluster)
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p, ".BG.csv")
      write.csv(data.frame(row.names(stations.rc), stations.rc$BG), saveas)
      
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p, "clusterplot.pdf")
      pdf(saveas,width=6,height=4) 
      plotcluster(tax.rc.norare, fit.km$cluster)
      dev.off()
      
      
      #clusplot(tax.rc.norare, fit.km$cluster, color=TRUE, shade=TRUE, 
      #         labels=2, lines=0)
      
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p, "_anosim")
      set.seed(10)
      anosim1<-anosim(tax.rc.norare, stations.rc$BG, permutations = 99, distance="bray")
      anosim.out<-paste(saveas, "R", anosim1$statistic, "Sig", anosim1$signif) 
      capture.output(anosim.out, file="OEout/anosim.out.txt", append=T)
      
      
      ###################################################################################
      # EVALUATE CANDIDATE VARS ----------
      ###################################################################################
      
      #import from file above
      source("bin/cand.var.bu.R")
      cand.var
      foo<-which(cand.var %in% names(stations.rc)) #do you have them all?
      cand.var<-cand.var[foo]
      #############################
      # CARET analysis to select env vars -------
      #############################
      
      # prep CARET dataframe 
      dat<-data.frame(stations.rc[,cand.var])
      foo<-complete.cases(dat)
      dat<-dat[foo,]
      dat.BG<-droplevels(stations.rc$BG)
      
      set.seed(10)
      rfe.out<-rfe(y=dat.BG, x=dat, sizes=2:predictornum,  rfeControl=ctrl, maximize=T)
      
      #record optvars
      optvars.rfe <- rfe.out$optVariables
      optvars.rfe
      
      #view results
      #ggplot(data=rfe.out$results, aes(x=Accuracy, y=Kappa))+
      #geom_text(aes(label=Variables))
      
      ##########################################################
      # RF with env vars from CARET -----------------
      ##########################################################
      
      #from RFE
      set.seed(10)
      rfe.out2<-randomForest(x = dat[,optvars.rfe], y = dat.BG, 
                             ntree=5000, importance=TRUE, norm.votes=TRUE, keep.forest=TRUE)
      print(rfe.out2)
      
      
      if (diatoms == T) { 
        diatom.rf.oe<-rfe.out2
        save(diatom.rf.oe, file="OEout/diatom.RF.OE.Rdata") 
        write.csv(row.names(diatom.rf.oe$importance), "OEout/diatoms.predvars.txt")
      }
      
      if (sba == T)  { 
        sba.rf.oe<-rfe.out2
        save(sba.rf.oe, file="OEout/sba.RF.OE.Rdata")
        write.csv(row.names(sba.rf.oe$importance), "sba.predvars.txt")}
      
      if (combo == T) { 
        hybrid.rf.oe<-rfe.out2
        save(hybrid.rf.oe, file="OEout/hybrid.RF.OE.Rdata") 
        write.csv(row.names(hybrid.rf.oe$importance), "hybrid.predvars.txt")}
      
      ############################################
      # Step 9 Calculate scores -----------------
      ############################################
      
      # need dataframes to be numeric
      if (prezabs==T) {tax.rc<-as.data.frame(ifelse(tax.rc>0,1,0))} 
      
      # save opt vars
      saveas<-paste0("OEout/","PredNums_",x,"_Clustnum_",y,"_Pcap_", p, "_RFEoptvars", ".txt")
      foo<-paste(saveas,unlist(optvars.rfe))
      capture.output(foo, file="OEout/optvars.rfe.all.csv")
      
      # save RFE scores
      set.seed(10)
      Scores.cal<-model.predict.RanFor.4.2(
        bugcal.pa = tax.rc,
        grps.final = dat.BG,
        preds.final = optvars.rfe, # for caret
        ranfor.mod = rfe.out2, # for caret
        prednew = stations.rc,
        bugnew = (tax.rc),
        Pc=Pcap,
        Cal.OOB=TRUE);
      
      
      # Accuracy: Intercept should be zero
      # Precision: Points should be tightly clustered around the regression line
      saveas2<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p,".pdf")
      ggplot(data=Scores.cal$OE.scores, aes(x=E, y=O), cex=2)+
        geom_point() + theme_bw() + xlab("Expected") + ylab("Observed") + 
        geom_smooth(method=lm, color="red") + xlim(0,20) + ylim(0,20) + geom_abline(intercept=0, slope=1) #Observed prediction
      if (plotz == T) { ggsave(file=saveas2, width = 6, height = 6) }
      
      #plot the O/E for ref sites by region 
      scores.out<-as.data.frame(cbind(Scores.cal$OE.scores, stations.rc$psa6c))
      saveas2<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"PSAregion.pdf")
      ggplot(scores.out, aes(y=as.numeric(scores.out$OoverE))) + 
        geom_boxplot(aes(x=factor(scores.out$`stations.rc$psa6c`))) +
        geom_jitter(aes(x=factor(scores.out$`stations.rc$psa6c`)), cex=0.5) +
        geom_hline(yintercept = 1) + theme_bw() + 
        xlab("") + ylab("O/E") 
      if (plotz == T) { ggsave(file=saveas2, width = 6, height = 6) }
      
      # Calculate replicate sampling SD of O/E
      #Then execute the function, using either in-bag or OOB predicted occurrence probs ('Capture probs') for the calibration data;
      saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_",p,"_repsamsd")
      repsamsd1<-rep.sam.sd(occprb.cal=Scores.cal$Capture.Probs,Pc=Pcap);
      repsamsd1.out<-paste(saveas2, repsamsd1)
      capture.output(print(repsamsd1.out), file="OEout/rep.sam.sd.txt", append=T)
      
      ######################################################
      # Apply model to reference validation  ---------
      ######################################################
      tax.rv<-as.data.frame(subset(tax.m, row.names(tax.m) %in% row.names(stations.rv), drop = T))
      stations.rv<-as.data.frame(subset(stations.rv, row.names(stations.rv) %in% row.names(tax.rv), drop=T))
      
      set.seed(10)
      Scores.rv<-model.predict.RanFor.4.2(
        bugcal=tax.rc,
        grps.final=dat.BG,
        preds.final=optvars.rfe, 
        ranfor.mod=rfe.out2,
        prednew=stations.rv,
        bugnew=(tax.rv),
        Pc=Pcap,
        Cal.OOB=F)
      
      saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_",p,"_repsamsd_rv")
      repsamsd2<-rep.sam.sd(occprb.cal=Scores.rv$Capture.Probs,Pc=Pcap);
      repsamsd2.out<-paste(saveas2, repsamsd2)
      capture.output(print(repsamsd2.out), file="OEout/rep.sam.sd.rv.txt", append=T)
      
      ####################################
      # Apply model to intermediate sites -------
      ####################################
      stations.int<-subset(stations, stations$SiteStatus=="Intermediate", drop = T)
      tax.int<-as.data.frame(subset(tax.m, row.names(tax.m) %in% row.names(stations.int), drop = T))
      foo<-which(stations.int$SeqID %in% row.names(tax.int))
      stations.int<-stations.int[foo,]
      foo<-complete.cases(stations.int[,optvars.rfe])
      stations.int<-stations.int[foo,]
      foo<-which(row.names(tax.int) %in% stations.int$SeqID)
      tax.int<-tax.int[foo,]
      
      set.seed(10)
      Scores.int<-model.predict.RanFor.4.2(
        bugcal=tax.rc,
        grps.final=dat.BG,
        preds.final=optvars.rfe, 
        ranfor.mod=rfe.out2,
        prednew=stations.int,
        bugnew=(tax.int),
        Pc=Pcap,
        Cal.OOB=F)
      
      ####################################
      # Apply model to stressed sites --------
      ####################################
      
      stations.str<-subset(stations, stations$SiteStatus=="Stressed", drop = T)
      stations.str<-stations.str[,optvars.rfe]
      foo<-which(is.na(rowSums(stations.str)))
      if (length(foo) > 0) { stations.str<-stations.str[-foo,] }
      
      tax.str<-as.data.frame(subset(tax.m, row.names(tax.m) %in% row.names(stations.str), drop = T))
      foo<-which(colSums(tax.str)==0)
      if (length(foo) > 0) {tax.str<-tax.str[,-foo] }
      
      stations.str<-subset(stations.str, row.names(stations.str) %in% row.names(tax.str), drop = T)
      
      set.seed(10)
      Scores.str<-model.predict.RanFor.4.2(
        bugcal=tax.rc,
        grps.final=dat.BG,
        preds.final=optvars.rfe, 
        ranfor.mod=rfe.out2,
        prednew=stations.str,
        bugnew=(tax.str),
        Pc=Pcap,
        Cal.OOB=F)
      
      ####################################
      # Apply model to ALL OTHERS sites --------
      ####################################
      
      foo<-which(stations.orig$MostRecent=="n")
      foo1<-stations.orig[foo,]
      
      foo<-which(is.na(stations.orig$SiteStatus))
      foo2<-stations.orig[-foo,]
      
      stations.all<-rbind(foo1, foo2)
      stations.all<-stations.all[,optvars.rfe]
      foo<-which(is.na(rowSums(stations.all)))
      stations.all<-stations.all[-foo,]
      
      tax.all<-tax.m
      tax.all<-ifelse(tax.all>0,1,0)
      foo<-which(colSums(tax.all)==0)
      if (length(foo) > 0) {tax.all<-tax.all[,-foo] }
      
      foo<-which(row.names(stations.all) %in% row.names(tax.all))
      if (length(foo)>0) {stations.all <-stations.all[foo,]}
      foo<-which(row.names(tax.all) %in% row.names(stations.all))
      if (length(foo)>0) {tax.all <-tax.all[foo,]}
      
      foo<-which(colSums(tax.all)==0)
      tax.all<-tax.all[,-foo]
      foo<-which(rowSums(tax.all)==0)
      
      tax.all<-as.data.frame(tax.all)
      tax.all<-tax.all[order(row.names(tax.all)),]
      stations.all<-stations.all[order(row.names(stations.all)),]
      
      row.names(tax.all) == row.names(stations.all)
      
      set.seed(10)
      Scores.all<-model.predict.RanFor.4.2(
        bugcal=tax.rc,
        grps.final=dat.BG,
        preds.final=optvars.rfe, 
        ranfor.mod=rfe.out2,
        prednew=stations.all,
        bugnew=(tax.all),
        Pc=Pcap,
        Cal.OOB=F)
      
      ####################################
      # Combined Scores --------
      ####################################
      
      Scores.cal<-as.data.frame(Scores.cal); Scores.cal$Type<-"Reference"
      Scores.rv<-as.data.frame(Scores.rv); Scores.rv$Type<-"Reference"
      Scores.int<-as.data.frame(Scores.int); Scores.int$Type<-"Intermediate"
      Scores.str<-as.data.frame(Scores.str); Scores.str$Type<-"Stressed"
      Scores.all<-as.data.frame(Scores.all); Scores.all$Type<-"Other"
      
      Scores.cat<-rbind(Scores.cal, Scores.rv, Scores.int, Scores.str, Scores.all)
      
      #type<-data.frame(row.names(stations.orig),stations.orig$SiteStatus, stations.orig$CalVal, stations.orig$StrCalVal )
      #row.names(type)<-type$row.names.stations.orig.
      
      #foo<-which(colnames(Scores.cat)=="Type")
      #Scores.cat<-Scores.cat[,-foo]
      
      #source("~/Documents/R/bin/cbind.na.R")
      #Scores.cat<-cbind.na(Scores.cat, type)
      
      ####################################
      # Plot ref/int/str --------
      ####################################
      
      # plot O/E by Ref/Int/Stressed 
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_OEbyType.pdf")
      ggplot(Scores.cat, aes(y=Scores.cat$OE.scores.OoverE, x=Scores.cat$Type)) +
        geom_boxplot(aes(fill=Scores.cat$Type)) + 
        geom_jitter(cex=0.02) +
        ylim(0,1.5) + scale_fill_manual(values=c("#FFCC66", "#339966", "#CC0000")) + 
        scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) +
        geom_hline(yintercept = 1) + theme_bw() + ylab("O/E") + xlab("")
      if (plotz == T) { ggsave(file=saveas, width = 6, height = 6) }
      
      
      ########################################
      # Stats -----------
      ########################################
      
      cat1<-Scores.cat
      cat1<-subset(cat1, Type!="Other")
      
      #calculating stats by ref/int/str
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_anova")
      set.seed(10)
      fit1<-aov(cat1$OE.scores.OoverE~cat1$Type); summary(fit1)
      anova.P<-summary(fit1)[[1]][["Pr(>F)"]][[1]]
      anova.F<-summary(fit1)[[1]][["F value"]][[1]]
      anova.out<-paste(saveas, "Pv", anova.P, "Fv", anova.F) 
      capture.output(anova.out, file="OEout/anova.out.txt", append=T)
      
      saveas<-paste0("OEout/","Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_ttest")
      set.seed(10)
      ttest1<-pairwise.t.test(cat1$OE.scores.OoverE, cat1$Type)
      ttest.RefInt<-ttest1$p.value[1]
      ttest.RefStr<-ttest1$p.value[2]
      ttest.IntStr<-ttest1$p.value[4]
      ttest.RefStr.t<-ttest1$p.value[2]
      ttest.out<-paste(saveas, "RefInt", ttest.RefInt, "RefStr", ttest.RefStr, "IntStr", ttest.IntStr) 
      capture.output(ttest.out, file="OEout/pw.ttest.out.txt", append=T)
      
      cat1.sub<-subset(cat1, Type!="Intermediate")
      set.seed(10)
      ttest1<-t.test(cat1.sub$OE.scores.OoverE~cat1.sub$Type)
      ttest.RefStr<-ttest1$statistic[1]
      ttest.out<-paste(saveas,"RefStr", ttest.RefStr) 
      capture.output(ttest.out, file="OEout/ttest.out.txt", append=T)
      
      
    }}}



########## Export files -----------

if (exportfiles==T) { 
  
  if(diatoms==T) {write.csv(stations.rc, "OEout/diatoms.stations.rc.csv")}
  if(diatoms==T) {write.csv(tax.rc, "OEout/diatoms.tax.rc.csv")}
  if(diatoms==T) {write.csv(stations, "OEout/diatoms.stations.csv")}
  if(diatoms==T) {write.csv(tax.m, "OEout/diatoms.tax.m.csv")}
  if (diatoms==T) { if (genuslev==F) { write.csv(tax.m, "OEout/diatoms.tax.m.species.csv")}}
  
  if(sba==T) {write.csv(stations.rc, "OEout/sba.stations.rc.csv")}
  if(sba==T) {write.csv(tax.rc, "OEout/sba.tax.rc.csv")}
  if(sba==T) {write.csv(stations, "OEout/sba.stations.csv")}
  if(sba==T) {write.csv(tax.m, "OEout/sba.tax.m.csv")}
  if (sba==T) { if (genuslev==F) { write.csv(tax.m, "OEout/sba.tax.m.species.csv")}}
  
  if(combo==T) {write.csv(stations.rc, "OEout/combo.stations.rc.csv")}
  if(combo==T) {write.csv(tax.rc, "OEout/combo.tax.rc.csv")}
  if(combo==T) {write.csv(stations, "OEout/combo.stations.csv")}
  if(combo==T) {write.csv(tax.m, "OEout/combo.tax.m.csv")}
  if(combo==T) { if (genuslev==F) { write.csv(tax.m, "OEout/combo.tax.m.species.csv")}}
  
  
  
  write.csv(Scores.cat, "OEout/Scores.all.csv")
  
  
}


# import the RFE scores into df # move RFEscores.txt into OE folder
rfe.final<-read.table("RFEscores.txt", sep = " ", stringsAsFactors = F)
rfe.final <- data.frame(do.call('rbind', strsplit(as.character(rfe.final$V2),' ',fixed=TRUE)))
colnames(rfe.final)<-c("combo", "model.mean", "model.stdev", "null.mean", "null.stdev", "dunno")
combo<-data.frame(do.call('rbind', strsplit(as.character(rfe.final$combo),'_',fixed=TRUE)))                               
#colnames(combo)<-c("a", "prednums", "b", "clusternums", "c", "pcap", "d")
#del<-c("a", "b", "c", "d")
#foo<-which(colnames(combo) %in% del)
#combo<-combo[,-foo]
rfe.final<-cbind(combo, rfe.final)                               
write.csv(rfe.final, "OEout/rfe.all.csv")

# import the ANOVA scores into df
anova.final<-read.table("OEout/anova.out.txt", sep = " ", stringsAsFactors = F)
anova.final <- data.frame(do.call('rbind', strsplit(as.character(anova.final$V2),' ',fixed=TRUE)))
colnames(anova.final)<-c("combo", "Pv", "pvalue", "Fv", "Fvalue")
combo<-data.frame(do.call('rbind', strsplit(as.character(anova.final$combo),'_',fixed=TRUE)))                               
#colnames(combo)<-c("a", "prednums", "b", "clusternums", "c", "pcap", "d")
#del<-c("a", "b", "c", "d")
#foo<-which(colnames(combo) %in% del)
#combo<-combo[,-foo]
anova.final<-cbind(combo, anova.final)                               
write.csv(anova.final, "OEout/anova.all.csv")

##### 
# make compare file 

#foo1<-read.csv("OEout/rfe.all.csv")
#foo2<-read.csv("OEout/anova.all.csv")

#foo1<-subset(foo1, foo1$null.mean==1)
#foo3<-cbind(foo1, foo2$Fvalue)
#foo3$diff<-(foo3$model.stdev - foo3$null.stdev)
#foo3<-subset(foo3, foo3$diff<0)

#write.csv(foo3, "OEout/compare.csv")
