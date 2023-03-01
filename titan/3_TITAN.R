# titan 
#install.packages("TITAN2")
library("TITAN2")
library("dplyr")
#glades.taxa<-glades.taxa
#glades.env<-glades.env

cyanos=F
sba=F
diatoms=T

# assemblage 
if (cyanos==T) { (assemblage="cyanos")} 
if (sba==T) { (assemblage="sba")} 
if (diatoms==T) { (assemblage="diatoms")} 

# cyanos 
if (cyanos==T) { 
  tax<-otu5.16S
  foo<-as.data.frame(do.call(rbind, strsplit(as.character(tax$ConsensusLineage), ";")))
  names(foo)<-c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
  foo$Phylum<-trimws(foo$Phylum, which="both")
  tax<-cbind(tax, foo)
  tax<-subset(tax,tax$Phylum =="Cyanobacteria") # only use diatom ASVs
  drop<-c("ConsensusLineage.x","EB_2_14_21","EB_7_26_21","EB_8_2_2021",
          "EB_9_1_20", "EB_9_9_21","ntc","zymodna","zymodna_logphase","ConsensusLineage", "Kingdom",                  
          "Phylum","Class","Order","Family","Genus","Species" )
  foo<-which(colnames(tax) %in% drop)
  tax<-tax[,-foo]
  
  write.csv(tax,"titan.asv/16S/tax.cyanos.csv")
  
  OTUID<-data.frame(tax$OTUID) 
  names(OTUID)<-"OTUID.orig"
  OTUID$num<-1:length(OTUID$OTUID.orig)
  OTUID$OTUID<-paste0("otu",OTUID$num)
  row.names(tax)<-OTUID$OTUID
  write.csv(OTUID, "titan.asv/16S/OTUID.cyanos.csv") 
  }

# diatoms
if (diatoms==T) { 
tax<-otu5.rbcL
foo<-as.data.frame(do.call(rbind, strsplit(as.character(tax$ConsensusLineage), ";")))
names(foo)<-c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
foo$Phylum<-trimws(foo$Phylum, which="both")
tax<-cbind(tax, foo)
tax<-subset(tax,tax$Phylum =="Bacillariophyta") # only use diatom ASVs
drop<-c("ConsensusLineage.x","EB_2_14_21","EB_7_26_21","EB_8_2_2021",
        "EB_9_1_20", "EB_9_9_21","ntc","zymodna","zymodna_logphase","ConsensusLineage", "Kingdom",                  
"Phylum","Class","Order","Family","Genus","Species" )
foo<-which(colnames(tax) %in% drop)
tax<-tax[,-foo]

write.csv(tax,"titan.asv/rbcL/tax.rbcL.csv")

OTUID<-data.frame(tax$OTUID) 
names(OTUID)<-"OTUID.orig"
OTUID$num<-1:length(OTUID$OTUID.orig)
OTUID$OTUID<-paste0("otu",OTUID$num)
row.names(tax)<-OTUID$OTUID
write.csv(OTUID, "titan.asv/rbcL/OTUID.diatoms.csv") 
}

# SBA
if (sba==T) { 
  tax<-otu5.18S
  foo<-as.data.frame(do.call(rbind, strsplit(as.character(tax$ConsensusLineage), ";")))
  names(foo)<-c("Kingdom", "Phylum", "Class","Order", "Family", "Genus", "Species")
  foo$Phylum<-trimws(foo$Phylum, which="both")
  tax<-cbind(tax, foo)
  #tax<-subset(tax,tax$Phylum =="Cyanobacteria") 
  drop<-c("ConsensusLineage.x","EB_2_14_21","EB_7_26_21","EB_8_2_2021",
          "EB_9_1_20", "EB_9_9_21","ntc","zymodna","zymodna_logphase","ConsensusLineage", "Kingdom",                  
          "Phylum","Class","Order","Family","Genus","Species" )
  foo<-which(colnames(tax) %in% drop)
  tax<-tax[,-foo]
  
  write.csv(tax,"titan.asv/18S/tax.18S.csv")
  
  OTUID<-data.frame(tax$OTUID) 
  names(OTUID)<-"OTUID.orig"
  OTUID$num<-1:length(OTUID$OTUID.orig)
  OTUID$OTUID<-paste0("otu",OTUID$num)
  row.names(tax)<-OTUID$OTUID
  write.csv(OTUID, "titan.asv/18S/OTUID.18S.csv") 
}

######################################3
# prep tax and env ------------

drop<-"OTUID"
foo<-which(colnames(tax) %in% drop)
tax<-tax[,-foo]

tax[is.na(tax)]<-0
tax<-as.data.frame(t(tax))
foo<-which(rownames(tax) %in% env$SeqID)
tax.sub<-tax[foo,]
foo<-which(env$SeqID %in% rownames(tax.sub))
env.sub<-env[foo,]

tax.sub<-tax.sub[order(row.names(tax.sub)),]
env.sub<-env.sub[order((env.sub$SeqID)),]

(rownames(tax.sub)) == (env.sub$SeqID)

# Run TITAN on all env vars ----------------------------------------------- 

params<-c("TN", "TP", 
          "Oxygen_Dissolved", 
          "AFDM",
          "DOC", 
          "SpCond",
          "ChlA",
          NULL)

ls.inc<-list()
ls.dec<-list()

for (i in params) { # i<-"AFDM"
foo<-complete.cases(env.sub[,i])
env.sub1<-env.sub[foo,]
tax.sub1<-tax.sub[foo,]

tax.sub1<-mutate_all(tax.sub1, function(x) as.numeric(as.character(x)))
foo<-which(rowSums(tax.sub1) ==0)
if(length(foo)>0){tax.sub1<-tax.sub1[-foo,]}
if(length(foo)>0){env.sub1<-env.sub1[-foo,]}

foo<-which(colSums(tax.sub1) ==0)
if(length(foo)>0){tax.sub1<-tax.sub1[,-foo]}
tax.sub1<-as.data.frame(tax.sub1)

tax.sub.pa<-ifelse(tax.sub1>0,1,0)
foo<-which(colSums(tax.sub.pa)<3)
if(length(foo)>0){tax.sub1<-tax.sub1[,-foo]}
foo<-which(rowSums(tax.sub1) ==0 )

titan.out<-titan(env = env.sub1[,i], txa = tax.sub1, 
                    minSplt = 5, numPerm = 25, boot = TRUE, nBoot = 20, imax = FALSE,
                    ivTot = FALSE, pur.cut = 0.95, rel.cut = 0.95, ncpus = 4, memory = FALSE)

saveas<-paste0("titan.asv/titan.out.sppmax.",i,".",assemblage,".csv")
write.csv(as.matrix(titan.out$sppmax), saveas)

saveas<-paste0("titan.out.",i,".",assemblage)
assign(saveas, titan.out)

titan.out$sumz.cp[3]
titan.out$sumz.cp[4]

ls.dec[i]<-titan.out$sumz.cp[3] # 0.12125 decreasers
ls.inc[i]<-titan.out$sumz.cp[4] # 2.285 increasers

# plots ----
saveas<-paste0("titan.asv/plots/titan.sumz.",i,".pdf")
pdf(file = saveas, width = 12, height = 10)
plot_sumz_density(titan.out, xlabel=i)
dev.off()

saveas<-paste0("titan.asv/plots/titan.taxa",i,".pdf")
pdf(file = saveas, width = 12, height = 10)
plot_taxa(titan.out, xlabel = i); 
dev.off()

saveas<-paste0("titan.asv/plots/titan.taxaridge.",i,".pdf")
pdf(file = saveas, width = 12, height = 10)
plot_taxa_ridges(titan.out, xlabel=i, cex=0.5)
dev.off()

}

foo1<-t(t(ls.dec))
foo2<-t(t(ls.inc))

out<-cbind(foo1, foo2)
names(out)<-c("Dec", "Inc")
if(diatoms==T) {write.csv(out, "titan.asv/rbcL/out.cp.csv")}
if(cyanos==T) {write.csv(out, "titan.asv/16S/out.cp.csv")}
if(sba==T) {write.csv(out, "titan.asv/18S/out.cp.csv")}


#########################################################
# Figures -------------------------

# # diatoms TN decreaser
# dat<-read.csv("titan.asv/rbcL/titan.out.sppmax.TN.diatoms.csv")
# colnames(dat)[colnames(dat) == "X"] ="OTUID"
# otuid<-read.csv("titan.asv/rbcL/OTUID.diatoms.csv")
# dat<-merge(dat, otuid, by="OTUID", sort=F, all.x=T, all.y=F)
# dat$OTUID<-dat$OTUID.orig
# dat<-subset(dat, dat$filter=='1')
# maxnum<-max(dat$IndVal)
# foo<-which(dat$IndVal == maxnum)
# foo2<-dat[foo,]
# keep<-as.character(droplevels(foo2$OTUID))
# 
# keep2<-which(otu5.rbcL$OTUID %in% keep)
# foo<-otu5.rbcL[keep2,]
# foo<-as.data.frame(t(foo))
# foo$SeqID<-row.names(foo)
# foo<-merge(foo, env[,c("SeqID", "TN")],by="SeqID", all.x=T, all.y=F, sort=F)
# plt<-foo
# foo<-complete.cases(plt)
# plt<-plt[foo,]
# plt[,2]<-as.numeric(as.character(plt[,2]))
# plt<-subset(plt, plt[,2]>0.000000001)
# 
# ggplot(plt, aes(x=log10(TN), y=as.numeric(plt[,2]))) + geom_point() + theme_bw() + 
#   geom_smooth(method=glm) + ylab(label = "Diatom ASV decreaser") +
#   ggsave("titan.asv/plots/diatoms.tn.dec.pdf")
# 
# # diatoms TN decreaser
# dat<-read.csv("titan.asv/rbcL/titan.out.sppmax.TN.diatoms.csv")
# colnames(dat)[colnames(dat) == "X"] ="OTUID"
# otuid<-read.csv("titan.asv/OTUID.diatoms.csv")
# dat<-merge(dat, otuid, by="OTUID", sort=F, all.x=T, all.y=F)
# dat$OTUID<-dat$OTUID.orig
# dat<-subset(dat, dat$filter=='2')
# maxnum<-max(dat$IndVal)
# foo<-which(dat$IndVal == maxnum)
# foo2<-dat[foo,]
# keep<-as.character(droplevels(foo2$OTUID))
# 
# keep2<-which(otu5.rbcL$OTUID %in% keep)
# foo<-otu5.rbcL[keep2,]
# foo<-as.data.frame(t(foo))
# foo$SeqID<-row.names(foo)
# foo<-merge(foo, env[,c("SeqID", "TN")],by="SeqID", all.x=T, all.y=F, sort=F)
# plt<-foo
# foo<-complete.cases(plt)
# plt<-plt[foo,]
# plt$`4`<-as.numeric(as.character(plt$`4`))
# plt<-subset(plt, plt$`4`>0.000000001)
# 
# ggplot(plt, aes(x=log10(TN), y=as.numeric(`4`))) + geom_point() + theme_bw() + 
#   geom_smooth(method=glm) + ylab(label = "Diatom ASV increaser") + 
#   ylim(0,0.2) + 
#   ggsave("titan.asv/plots/diatoms.tn.inc.pdf")



