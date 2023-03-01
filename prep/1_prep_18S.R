# mASCI

# otu
# things to do == trying the eDNA index normalization
otu1<-read.csv("data/otu/18S_plate1.csv", stringsAsFactors = F, row.names=1)
otu2<-read.csv("data/otu/18S_plate2.csv", stringsAsFactors = F, row.names=1)
otu3<-read.csv("data/otu/18S_2022.csv", stringsAsFactors = F, row.names=1)

# note singletons 
otu1.pa<-otu1
foo<-ifelse(otu1[,1:90]>0,1,0)
otu1.pa[,1:90]<-foo
sings1<-which(rowSums(otu1.pa[,1:90])==1)

otu2.pa<-otu2
foo<-ifelse(otu2[,1:85]>0,1,0)
otu2.pa[,1:85]<-foo
sings2<-which(rowSums(otu2.pa[,1:85])==1)

otu3.pa<-otu3
foo<-ifelse(otu3[,1:243]>0,1,0)
otu3.pa[,1:243]<-foo
sings3<-which(rowSums(otu3.pa[,1:243])==1)

# make rel abd
foo<-sweep(otu1[,1:90], 2, colSums(otu1[,1:90]), '/')
otu1[,1:90]<-foo
foo<-sweep(otu2[,1:85], 2, colSums(otu2[,1:85]), '/')
otu2[,1:85]<-foo
foo<-sweep(otu3[,1:243], 2, colSums(otu3[,1:243]), '/')
otu3[,1:243]<-foo



# remove singletons 
otu1<-otu1[-sings1,]
otu2<-otu2[-sings2,]
otu3<-otu3[-sings3,]

otu1$OTUID<-row.names(otu1)
otu2$OTUID<-row.names(otu2)
otu3$OTUID<-row.names(otu3)

otu4<-merge(otu1, otu2, by="OTUID", sort=F, all.x=T, all.y=T)
otu5<-merge(otu4, otu3, by="OTUID", sort=F, all.x=T, all.y=T)
#write.csv(otu5, "data/otu/otu.merged.csv")

foo<-which(otu4$OTUID %in% otu3$OTUID)
incommon<-otu4[foo,]
foo<-which(otu5$OTUID %in% incommon$OTUID)
incommon.merged<-otu5[foo,]
#write.csv(incommon.merged, "data/otu/otu.merged.incommon.csv")

otu5[is.na(otu5)]<-0

otu5.18S<-otu5

write.csv(otu5.18S, "data/otu/otu5.18S.csv")

# combine OTU tables 
otu.m<-otu5
otu.m<-melt(otu.m)
otu.m$SeqID<-otu.m$variable

foo<-which(is.na(otu.m$value))
#otu.m<-otu.m[-foo,]
foo<-which(otu.m$value==0.000000000000000000)
otu.m<-otu.m[-foo,]

drop<-c("ConsensusLineage.x" ,"ConsensusLineage.y", "ConsensusLineage")
foo<-which(colnames(otu.m) %in% drop)
otu.m<-otu.m[,-foo]

# consensus lineage
cl<-rbind(otu1[,c("OTUID","ConsensusLineage")], otu2[,c("OTUID","ConsensusLineage")])
cl<-rbind(cl, otu3[,c("OTUID","ConsensusLineage")])
cl<-unique(cl)

otu.m<-merge(otu.m, cl, by="OTUID", sort=F)
foo<-as.data.frame(do.call(rbind, strsplit(as.character(otu.m$ConsensusLineage), "; ")))
names(foo)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu.m<-cbind(otu.m, foo)

# rename 
foo<-env[,c("SeqID","stationcode.susie",
            "masterid","Date","SiteStatus", "Replicate")]
xwalk<-unique(foo)
otu.m<-merge(otu.m, xwalk, by="SeqID", all.x = T, all.y=T, sort=F)

otu.m.18S<-otu.m

# save 
write.csv(otu.m,"data/otu/otu.m.18S.csv")






