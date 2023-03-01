
# Take titan output and make traits table 

t.spcond<-read.csv("titan.asv/rbcL/titan.out.sppmax.SpCond.diatoms.csv")
t.doc<-read.csv("titan.asv/rbcL/titan.out.sppmax.DOC.diatoms.csv")
t.afdm<-read.csv("titan.asv/rbcL/titan.out.sppmax.AFDM.diatoms.csv")
t.chl<-read.csv("titan.asv/rbcL/titan.out.sppmax.ChlA.diatoms.csv")
t.do<-read.csv("titan.asv/rbcL/titan.out.sppmax.Oxygen_Dissolved.diatoms.csv")
t.tp<-read.csv("titan.asv/rbcL/titan.out.sppmax.TP.diatoms.csv")
t.tn<-read.csv("titan.asv/rbcL/titan.out.sppmax.TN.diatoms.csv")

otuid<-read.csv("titan.asv/rbcL/OTUID.diatoms.csv")

# renmae
colnames(t.spcond)[colnames(t.spcond) == "filter"] = "trait_spcond"
colnames(t.doc)[colnames(t.doc) == "filter"] = "trait_doc"
colnames(t.afdm)[colnames(t.afdm) == "filter"] = "trait_afdm"
colnames(t.chl)[colnames(t.chl) == "filter"] = "trait_chl"
colnames(t.do)[colnames(t.do) == "filter"] = "trait_do"
colnames(t.tp)[colnames(t.tp) == "filter"] = "trait_tp"
colnames(t.tn)[colnames(t.tn) == "filter"] = "trait_tn"

colnames(t.spcond)[colnames(t.spcond) == "zenv.cp"] = "zenv.cp_spcond"
colnames(t.doc)[colnames(t.doc) == "zenv.cp"] = "zenv.cp_doc"
colnames(t.afdm)[colnames(t.afdm) == "zenv.cp"] = "zenv.cp_afdm"
colnames(t.chl)[colnames(t.chl) == "zenv.cp"] = "zenv.cp_chl"
colnames(t.do)[colnames(t.do) == "zenv.cp"] = "zenv.cp_do"
colnames(t.tp)[colnames(t.tp) == "zenv.cp"] = "zenv.cp_tp"
colnames(t.tn)[colnames(t.tn) == "zenv.cp"] = "zenv.cp_tn"

# get thr
thr.dec<-titan.out.SpCond.diatoms$sumz.cp[15]
thr.inc<-titan.out.SpCond.diatoms$sumz.cp[16]
drop1<-which(t.spcond$zenv.cp_spcond>thr.dec & t.spcond$trait_spcond==1)
drop2<-which(t.spcond$zenv.cp_spcond<thr.inc & t.spcond$trait_spcond==2)
if(length(drop1)>0) {t.spcond<-t.spcond[-drop1,]}
if(length(drop2)>0) {t.spcond<-t.spcond[-drop2,]}
t.spcond<-subset(t.spcond, t.spcond$trait_spcond !=0, drop = T)

thr.dec<-titan.out.DOC.diatoms$sumz.cp[15]
thr.inc<-titan.out.DOC.diatoms$sumz.cp[16]
drop1<-which(t.doc$zenv.cp_doc>thr.dec & t.doc$trait_doc==1)
drop2<-which(t.doc$zenv.cp_doc<thr.inc & t.doc$trait_doc==2)
if(length(drop1)>0) {t.doc<-t.doc[-drop1,]}
if(length(drop2)>0) {t.doc<-t.doc[-drop2,]}
t.doc<-subset(t.doc, t.doc$trait_doc !=0, drop = T)

thr.dec<-titan.out.AFDM.diatoms$sumz.cp[15]
thr.inc<-titan.out.AFDM.diatoms$sumz.cp[16]
drop1<-which(t.afdm$zenv.cp_afdm>thr.dec & t.afdm$trait_afdm==1)
drop2<-which(t.afdm$zenv.cp_afdm<thr.inc & t.afdm$trait_afdm==2)
if(length(drop1)>0) {t.afdm<-t.afdm[-drop1,]}
if(length(drop2)>0) {t.afdm<-t.afdm[-drop2,]}
t.afdm<-subset(t.afdm, t.afdm$trait_afdm !=0, drop = T)

thr.dec<-titan.out.ChlA.diatoms$sumz.cp[15]
thr.inc<-titan.out.ChlA.diatoms$sumz.cp[16]
drop1<-which(t.chl$zenv.cp_chl>thr.dec & t.chl$trait_chl==1)
drop2<-which(t.chl$zenv.cp_chl<thr.inc & t.chl$trait_chl==2)
if(length(drop1)>0) {t.chl<-t.chl[-drop1,]}
if(length(drop2)>0) {t.chl<-t.chl[-drop2,]}
t.chl<-subset(t.chl, t.chl$trait_chl !=0, drop = T)

thr.dec<-titan.out.Oxygen_Dissolved.diatoms$sumz.cp[15]
thr.inc<-titan.out.Oxygen_Dissolved.diatoms$sumz.cp[16]
drop1<-which(t.do$zenv.cp_do>thr.dec & t.do$trait_do==1)
drop2<-which(t.do$zenv.cp_do<thr.inc & t.do$trait_do==2)
if(length(drop1)>0) {t.do<-t.do[-drop1,]}
if(length(drop2)>0) {t.do<-t.do[-drop2,]}
t.do<-subset(t.do, t.do$trait_do !=0, drop = T)

thr.dec<-titan.out.TN.diatoms$sumz.cp[15]
thr.inc<-titan.out.TN.diatoms$sumz.cp[16]
drop1<-which(t.tn$zenv.cp_tn>thr.dec & t.tn$trait_tn==1)
drop2<-which(t.tn$zenv.cp_tn<thr.inc & t.tn$trait_tn==2)
if(length(drop1)>0) {t.tn<-t.tn[-drop1,]}
if(length(drop2)>0) {t.tn<-t.tn[-drop2,]}
t.tn<-subset(t.tn, t.tn$trait_tn !=0, drop = T)

thr.dec<-titan.out.TP.diatoms$sumz.cp[15]
thr.inc<-titan.out.TP.diatoms$sumz.cp[16]
drop1<-which(t.tp$zenv.cp_tp>thr.dec & t.tp$trait_tp==1)
drop2<-which(t.tp$zenv.cp_tp<thr.inc & t.tp$trait_tp==2)
if(length(drop1)>0) {t.tp<-t.tp[-drop1,]}
if(length(drop2)>0) {t.tp<-t.tp[-drop2,]}
t.tp<-subset(t.tp, t.tp$trait_tp !=0, drop = T)

# merge
traits.asv<-merge(t.spcond[,c(1,3,17)], t.doc[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)
traits.asv<-merge(traits.asv, t.afdm[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)
traits.asv<-merge(traits.asv, t.chl[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)
traits.asv<-merge(traits.asv, t.do[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)
traits.asv<-merge(traits.asv, t.tp[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)
traits.asv<-merge(traits.asv, t.tn[,c(1,3,17)], by="X", sort=F, all.x=T, all.y = T)

colnames(traits.asv)[colnames(traits.asv) == "X"] = "OTUID"

traits.asv<-merge(traits.asv, otuid[c("OTUID", "OTUID.orig")], by="OTUID", sort=F)

write.csv(traits.asv, "traits/traits.asv.csv")

traits.asv.cl<-traits.asv
traits.asv.cl$OTUID<-traits.asv.cl$OTUID.orig
traits.asv.cl<-merge(traits.asv.cl, otu5.rbcL[,c("OTUID", "ConsensusLineage")], by="OTUID",all.x=T, all.y=F, sort=F)
write.csv(traits.asv.cl, "titan.asv/rbcL/traits.asv.cl.csv")


# ###############
# # with rafi 
# 
# titan.out.AFDM.diatoms$sppmax
# titan.out.AFDM.diatoms$sumz.cp # f is filtered; look for dec decrease < 3.6, inc that increase > 107
# 
# foo<-
#   titan.out.AFDM.diatoms$sppmax %>% 
#   as_tibble() %>% 
#   mutate(taxon=row.names(titan.out.AFDM.diatoms$sppmax)) %>% 
#   inner_join(
#     tibble(maxgrp=c(1,2,1,2),
#            filtered=c(F,F,T,T), 
#            meancp=titan.out.AFDM.diatoms$sumz.cp[,4]) %>%
#       filter(filtered) %>%
#       select(-filtered)
#   ) %>%
#   mutate(trait.afdm=case_when(
#     #  maxgrp==1 & filter==1 & zenv.cp<meancp~"sensitive.dec",
#     #  maxgrp==2 & filter==2 & zenv.cp>meancp~"super.inc",
#     maxgrp==1 & zenv.cp<meancp~"sensitive.dec",
#     maxgrp==2 & zenv.cp>meancp~"super.inc",
#     T~"other"
#   ))
# 
# foo %>% group_by(trait.afdm) %>% tally()
# foo


