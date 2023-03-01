# calculating metrics for MMI 

library(sqldf)
library(dplyr)
library(plyr)

dir.create("R/mmi/diatoms/")
dir.create("R/mmi/diatoms/diatoms.plots")
dir.create("R/mmi/sba/")
dir.create("R/mmi/sba/sba.plots")
dir.create("R/mmi/hybrid/")
dir.create("R/mmi/hybrid/hybrid.plots")

###Loading Data
stations=read.csv('data/env/env5.csv',header=TRUE,strip.white=TRUE,check.names=FALSE) #this is most recent + the not most recent sites
stations<-stations[order(stations$SeqID),]
stations<-subset(stations, stations$SeqID %in% colnames(otu))
stations$SalinityPPT<-stations$Salinity
foo<-which(colnames(stations) %in% "Salinity")
stations<-stations[,-foo]

if(ASV==F){
traits=read.csv('traits/traits.new.csv',header=TRUE,strip.white=TRUE,check.names=FALSE, na.strings=c("","NA"))
traits$ConsensusLineage<-traits$FinalIDassigned
traits$Species<-traits$FinalIDassigned}

if(ASV==T){
  traits=read.csv('traits/traits.ASV.csv',header=TRUE,strip.white=TRUE,check.names=FALSE, na.strings=c("","NA"))
  traits$Species<-traits$OTUID.orig}

##################################################################################################################################################################
# Indicators
##################################################################################################################################################################

###Finding Indicators for Stressed/Reference Sites; If already know load data#####################################################################################
indicators.count=read.csv('R/mmi/indic.with.count.csv',header=TRUE,strip.white=TRUE,check.names=FALSE)
indicators=indicators.count
indicators$Species<-indicators$FinalIDassigned
##################################################################################################################################################################
# data prep 
##################################################################################################################################################################
###General Data Prep For Evaluating Metrics

taxonomy_pa_melt=melt(otu,id=c('OTUID', "ConsensusLineage", "Species")) #melting data into long format
taxonomy_pa_melt$value<-ifelse(taxonomy_pa_melt$value>0,1,0)
names(taxonomy_pa_melt)<-c("OTUID", "ConsensusLineage","Species","SeqID", "value")
taxonomy_pa_melt<-subset(taxonomy_pa_melt, value>0)

###Other Data Prep
#other_stations=sqldf("select * from stations where SiteSetSample2 = 'NA'") #selecting only NA stations
other_stations<-stations
#the next set of lines is the combining of tables to create one giant table that is to be used in the metrics calculations below

# other_combined1=sqldf("select * from (other_stations left join taxonomy_pa_melt using(SeqID)) left join traits using(Species)")
other_combined1=merge(other_stations, taxonomy_pa_melt, by="SeqID", sort=F, all.x=T, all.y=F)
other_combined1=merge(other_combined1, traits, by="Species", sort=F, all.x=T, all.y=F)
other_combined2=merge(other_combined1, indicators, by="Species", sort=F, all.x=T, all.y=F)

# other_combined3=sqldf("select * from other_combined2
#                          left join
#                          taxonomy_count_melt
#                          on other_combined2.SeqID=taxonomy_count_melt.SeqID
#                          and other_combined2.FinalIDassigned=taxonomy_count_melt.variable")
other_combined3<-other_combined2

#the next two lines create a genus variable by splitting the full taxa name
other_list=strsplit(other_combined3$Species," ")
other_combined3$Genus=lapply(other_list,FUN=function(x) x[1])

#duplicated(other_combined3)
other_combined3<-other_combined3[!duplicated(other_combined3), ] 

if(ASV==F) {
other_combined3.d<-subset(other_combined3, Phylum=="Bacillariophyta")
other_combined3.sba<-subset(other_combined3, Phylum!="Bacillariophyta")}

if(ASV==T){ 
  other_combined3.d<-other_combined3
  other_combined3.sba<-other_combined3}

################################################################################
# Metric calculations --------
################################################################################

# ASV 

# "trait_spcond" "trait_doc"    "trait_afdm"   "trait_chl"   
# "trait_do"     "trait_tp"     "trait_tn"

if (ASV==T) {
  
  if(diatoms ==T ){ 
  metrics.d=plyr::ddply(.data = other_combined3.d, .variables = ~SeqID,.fun = summarise,
   
    trait_tn.1.richness=sum(na.omit(trait_tn=='1')), 
    trait_tn.2.richness=sum(na.omit(trait_tn=='2')),
    trait_tn.1.cnt.spp=sum(na.omit(trait_tn=='1')),
    trait_tn.2.cnt.spp=sum(na.omit(trait_tn=='2')),
    trait_tn.1.prop.spp=sum(na.omit(trait_tn=='1'))/length(na.omit(trait_tn)),
    trait_tn.2.prop.spp=sum(na.omit(trait_tn=='2'))/length(na.omit(trait_tn)),
    
    trait_tp.1.richness=sum(na.omit(trait_tp=='1')), 
    trait_tp.2.richness=sum(na.omit(trait_tp=='2')),
    trait_tp.1.cnt.spp=sum(na.omit(trait_tp=='1')),
    trait_tp.2.cnt.spp=sum(na.omit(trait_tp=='2')),
    trait_tp.1.prop.spp=sum(na.omit(trait_tp=='1'))/length(na.omit(trait_tp)),
    trait_tp.2.prop.spp=sum(na.omit(trait_tp=='2'))/length(na.omit(trait_tp)),
    
    trait_do.1.richness=sum(na.omit(trait_do=='1')), 
    trait_do.2.richness=sum(na.omit(trait_do=='2')),
    trait_do.1.cnt.spp=sum(na.omit(trait_do=='1')),
    trait_do.2.cnt.spp=sum(na.omit(trait_do=='2')),
    trait_do.1.prop.spp=sum(na.omit(trait_do=='1'))/length(na.omit(trait_do)),
    trait_do.2.prop.spp=sum(na.omit(trait_do=='2'))/length(na.omit(trait_do)),
    
    trait_doc.1.richness=sum(na.omit(trait_doc=='1')), 
    trait_doc.2.richness=sum(na.omit(trait_doc=='2')),
    trait_doc.1.cnt.spp=sum(na.omit(trait_doc=='1')),
    trait_doc.2.cnt.spp=sum(na.omit(trait_doc=='2')),
    trait_doc.1.prop.spp=sum(na.omit(trait_doc=='1'))/length(na.omit(trait_doc)),
    trait_doc.2.prop.spp=sum(na.omit(trait_doc=='2'))/length(na.omit(trait_doc)),
    
    trait_spcond.1.richness=sum(na.omit(trait_spcond=='1')), 
    trait_spcond.2.richness=sum(na.omit(trait_spcond=='2')),
    trait_spcond.1.cnt.spp=sum(na.omit(trait_spcond=='1')),
    trait_spcond.2.cnt.spp=sum(na.omit(trait_spcond=='2')),
    trait_spcond.1.prop.spp=sum(na.omit(trait_spcond=='1'))/length(na.omit(trait_spcond)),
    trait_spcond.2.prop.spp=sum(na.omit(trait_spcond=='2'))/length(na.omit(trait_spcond)),
   
    trait_afdm.1.richness=sum(na.omit(trait_afdm=='1')), 
    trait_afdm.2.richness=sum(na.omit(trait_afdm=='2')),
    trait_afdm.1.cnt.spp=sum(na.omit(trait_afdm=='1')),
    trait_afdm.2.cnt.spp=sum(na.omit(trait_afdm=='2')),
    trait_afdm.1.prop.spp=sum(na.omit(trait_afdm=='1'))/length(na.omit(trait_afdm)),
    trait_afdm.2.prop.spp=sum(na.omit(trait_afdm=='2'))/length(na.omit(trait_afdm)),
    
    trait_chl.1.richness=sum(na.omit(trait_chl=='1')), 
    trait_chl.2.richness=sum(na.omit(trait_chl=='2')),
    trait_chl.1.cnt.spp=sum(na.omit(trait_chl=='1')),
    trait_chl.2.cnt.spp=sum(na.omit(trait_chl=='2')),
    trait_chl.1.prop.spp=sum(na.omit(trait_chl=='1'))/length(na.omit(trait_chl)),
    trait_chl.2.prop.spp=sum(na.omit(trait_chl=='2'))/length(na.omit(trait_chl))
    ) 
    
}}


### other Metric Calculations ------------------------------------------------------------------------------------------------------------

if (ASV ==F ) { 
if (diatoms == T) {
metrics.d=plyr::ddply(.data = other_combined3.d, .variables = ~SeqID,.fun = summarise,
              
              # richness 
              OrgN.NAHON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NAHON')), #count of NAHON - species
              OrgN.NALON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NALON')), #count of NALON - species
              OrgN.NHHONF.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF')), #count of NHHONF - species
              OrgN.NHHONForNHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONF and NHHONO - species
              OrgN.NHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONO - species
              OxyReq.DO_100.richness=sum(na.omit(OxygenRequirements=='DO_100')), #count of DO_100 - species
              OxyReq.DO_100orDO_75.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75')), #count of DO_100 and DO_75 - species
              OxyReq.DO_100orDO_75orDO_50.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'|OxygenRequirements=='DO_50')), #count of DO_100 and DO_75 and DO_50 - species
              OxyReq.DO_75.richness=sum(na.omit(OxygenRequirements=='DO_75')), #count of DO_75 - species
              OxyReq.DO_50.richness=sum(na.omit(OxygenRequirements=='DO_50')), #count of DO_50 - species
              OxyRed.DO_30.richness=sum(na.omit(OxygenRequirements=='DO_30')), #count of DO_30 - species
              OxyReq.DO_30andDO_10.richness=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10')), #count of DO_30 and DO_10 - species
              OxyReq.DO_10.richness=sum(na.omit(OxygenRequirements=='DO_10')), #count of DO_10 - species
              Salinity.B.richness=sum(na.omit(Salinity=='B')), #count of B - species
              Salinity.BF.richness=sum(na.omit(Salinity=='BF')), #count of BF - species
              Salinity.F.richness=sum(na.omit(Salinity=='F')), #count of F - species
              Salinity.FB.richness=sum(na.omit(Salinity=='FB')), #count of FB - species
              Saprobic.AM.richness=sum(na.omit(Saprobity=='AM')), #count of AM - species
              Saprobic.AMorAMPS.richness=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS')), #count of AM and AMPS - species
              Saprobic.AMPS.richness=sum(na.omit(Saprobity=='AMPS')), #count of AMPS - species
              Saprobic.BM.richness=sum(na.omit(Saprobity=='BM')), #count of BM - species
              Saprobic.OS.richness=sum(na.omit(Saprobity=='OS')), #count of OS - species
              Saprobic.OSandPS.richness=sum(na.omit(Saprobity=='OS'|Saprobity=='PS')), #count of OS and PS - species
              Saprobic.PS.richness=sum(na.omit(Saprobity=='PS')), #count of PS - species
              Trophic.E.richness=sum(na.omit(TrophicState=='E')), #count of E - species
              Trophic.EorI.richness=sum(na.omit(TrophicState=='E'|TrophicState=='I')), #count of E and I - species
              Trophic.I.richness=sum(na.omit(TrophicState=='I')), #count of I - species
              Trophic.M.richness=sum(na.omit(TrophicState=='M')), #count of M - species
              Trophic.ME.richness=sum(na.omit(TrophicState=='ME')), #count of ME - species
              Trophic.O.richness=sum(na.omit(TrophicState=='O')), #count of O - species
              Trophic.OorOMorPH.richness=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH')), #count of O, OM, and PH - species
              Trophic.OM.richness=sum(na.omit(TrophicState=='OM')), #count of OM
              Trophic.PH.richness=sum(na.omit(TrophicState=='PH')), #count of PH
              Achnanthes.richness=sum(na.omit(Genus=='Achnanthes')), #count of genus Achnanthes
              Amphora.richness=sum(na.omit(Genus=='Amphora')), #count of genus Amphora
              Cocconeis.richness=sum(na.omit(Genus=='Cocconeis')), #count of genus Cocconeis
              Cyclotella.richness=sum(na.omit(Genus=='Cyclotella')), #count of genus Cyclotella
              Cymbella.richness=sum(na.omit(Genus=='Cymbella')), #count of genus Cymbella
              Epithemia.richness=sum(na.omit(Genus=='Epithemia')), #count of genus Epithemia
              Eunotia.richness=sum(na.omit(Genus=='Eunotia')), #count of genus Eunotia
              Fragilaria.richness=sum(na.omit(Genus=='Fragilaria')), #count of genus Fragilaria
              Frustulia.richness=sum(na.omit(Genus=='Frustulia')), #count of genus Frustulia
              Gomphonema.richness=sum(na.omit(Genus=='Gomphonema')), #count of genus Gomphonema
              Navicula.richness=sum(na.omit(Genus=='Navicula')), #count of genus Navicula
              Nitzschia.richness=sum(na.omit(Genus=='Nitzschia')), #count of genus Nitzschia
              Rhoicosphenia.richness=sum(na.omit(Genus=='Rhoicosphenia')), #count of genus Rhoicosphenia
              Rhopalodia.richness=sum(na.omit(Genus=='Rhopalodia')), #count of genus Rhopalodia
              Surirella.richness=sum(na.omit(Genus=='Surirella')), #count of genus Surirella
              Synedra.richness=sum(na.omit(Genus=='Synedra')), #count of genus Synedra
              EpiRho.richness=sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #count of genus Epithemia and Rhopalodia
              
              # count ----------------------------------------
              cnt.spp.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella')), #count of NNS
              cnt.spp.HighMotility=sum(na.omit(Motility=='HM')), #count of high motility - species
              cnt.spp.ModMotility=sum(na.omit(Motility=='MM')), #count of moderate motility - species
              cnt.spp.MinMotility=sum(na.omit(Motility=='LM')), #count of low motility - species
              cnt.spp.NonMotile=sum(na.omit(Motility=='NM')), #count of no motility - species
              cnt.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high')), #from betty
              cnt.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low')), #from betty
              cnt.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low')), #from betty
              #cnt.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high')), 
              #cnt.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high')), 
              #cnt.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF')), 
              #cnt.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF')), 
              cnt.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high')), 
              #cnt.spp.Heterocy=sum(na.omit(Heterocystous=='yes')), 
              cnt.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph')), 
              cnt.spp.Halo=sum(na.omit(Salinity2=='Halo')), 
              #cnt.spp.ZHR=sum(na.omit(ZHR=='yes')), 
              #cnt.spp.CRUS=sum(na.omit(CRUS=='yes')), 
              #cnt.spp.Green=sum(na.omit(Green=='yes')), 
              cnt.spp.Planktonic=sum(na.omit(Habitat=='P')), 
              cnt.spp.least.tol=sum(na.omit(designation=='Reference')), #cnt of species that are least tolerant
              cnt.spp.most.tol=sum(na.omit(designation=='Stressed')), #cnt of species that are most tolerant
              cnt.spp.BCG1=sum(na.omit(BCG=='1')), #count of BCG 1 - species
              cnt.spp.BCG2=sum(na.omit(BCG=='2')), #count of BCG 2 - species
              cnt.spp.BCG3=sum(na.omit(BCG=='3')), #count of BCG 3 - species
              cnt.spp.BCG4=sum(na.omit(BCG=='4')), #count of BCG 4 - species
              cnt.spp.BCG5=sum(na.omit(BCG=='5')), #count of BCG 5 - species
              cnt.spp.BCG6=sum(na.omit(BCG=='6')), #count of BCG 6 - species
              cnt.spp.BCG12=sum(na.omit(BCG %in% c('1', "2"))), 
              cnt.spp.BCG45=sum(na.omit(BCG %in% c('4', "5"))), 
              
              # proportion -----------------------------------
              prop.spp.OrgN.NAHON=sum(na.omit(NitrogenUptakeMetabolism=='NAHON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NAHON - species
              prop.spp.OrgN.NALON=sum(na.omit(NitrogenUptakeMetabolism=='NALON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NALON - species
              prop.spp.OrgN.NHHONF=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF - species
              prop.spp.OrgN.NHHONForNHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF and NHHONO - species
              prop.spp.OrgN.NHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONO - species
              prop.spp.OxyReq.DO_100=sum(na.omit(OxygenRequirements=='DO_100'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 - species
              prop.spp.OxyReq.DO_100orDO_75=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 and DO_75 - species
              prop.spp.OxyReq.DO_75=sum(na.omit(OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_75 - species
              prop.spp.OxyReq.DO_50=sum(na.omit(OxygenRequirements=='DO_50'))/length(na.omit(OxygenRequirements)), #proportion of DO_50 - species
              prop.spp.OxyReq.DO_atleast50=sum(na.omit(OxygenRequirements %in% c('DO_50','DO_75','DO_100')))/length(na.omit(OxygenRequirements)), #proportion of at least DO_50 - species
              prop.spp.OxyRed.DO_30=sum(na.omit(OxygenRequirements=='DO_30'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 - species
              prop.spp.OxyReq.DO_30orDO_10=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 and DO_10 - species
              prop.spp.OxyReq.DO_10=sum(na.omit(OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_10 - species
              prop.spp.Salinity.B=sum(na.omit(Salinity=='B'))/length(na.omit(Salinity)), #proportion of B - species
              prop.spp.Salinity.BF=sum(na.omit(Salinity=='BF'))/length(na.omit(Salinity)), #proportion of BF - species
              prop.spp.Salinity.F=sum(na.omit(Salinity=='F'))/length(na.omit(Salinity)), #proportion of F - species
              prop.spp.Salinity.FB=sum(na.omit(Salinity=='FB'))/length(na.omit(Salinity)), #proportion of FB - species
              prop.spp.Saprobic.AM=sum(na.omit(Saprobity=='AM'))/length(na.omit(Saprobity)), #proportion of AM - species
              prop.spp.Saprobic.AMorAMPS=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AM and AMPS - species
              prop.spp.Saprobic.AMPS=sum(na.omit(Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AMPS - species
              prop.spp.Saprobic.BM=sum(na.omit(Saprobity=='BM'))/length(na.omit(Saprobity)), #proportion of BM - species
              prop.spp.Saprobic.OS=sum(na.omit(Saprobity=='OS'))/length(na.omit(Saprobity)), #proportion of OS - species
              prop.spp.Saprobic.OSorPS=sum(na.omit(Saprobity=='OS'|Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of OS and PS - species
              prop.spp.Saprobic.PS=sum(na.omit(Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of PS - species
              prop.spp.Trophic.E=sum(na.omit(TrophicState=='E'))/length(na.omit(TrophicState)), #proportion of E - species
              prop.spp.Trophic.EorI=sum(na.omit(TrophicState=='E'|TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of E and I - species
              prop.spp.Trophic.I=sum(na.omit(TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of I - species
              prop.spp.Trophic.M=sum(na.omit(TrophicState=='M'))/length(na.omit(TrophicState)), #proportion of M - species
              prop.spp.Trophic.ME=sum(na.omit(TrophicState=='ME'))/length(na.omit(TrophicState)), #proportion of ME - species
              prop.spp.Trophic.O=sum(na.omit(TrophicState=='O'))/length(na.omit(TrophicState)), #proportion of O - species
              prop.spp.Trophic.OorOMorPH=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of O, OM, and PH - species
              prop.spp.Trophic.OM=sum(na.omit(TrophicState=='OM'))/length(na.omit(TrophicState)), #proportion of OM
              prop.spp.Trophic.PH=sum(na.omit(TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of PH
              
              prop.Achnanthes=sum(na.omit(Genus=='Achnanthes'))/length(na.omit(Genus)), #proportion of genus Achnanthes
              prop.Amphora=sum(na.omit(Genus=='Amphora'))/length(na.omit(Genus)), #proportion of genus Amphora
              prop.Cocconeis=sum(na.omit(Genus=='Cocconeis'))/length(na.omit(Genus)), #proportion of genus Cocconeis
              prop.Cyclotella=sum(na.omit(Genus=='Cyclotella'))/length(na.omit(Genus)), #proportion of genus Cyclotella
              prop.Cymbella=sum(na.omit(Genus=='Cymbella'))/length(na.omit(Genus)), #proportion of genus Cymbella
              prop.Epithemia=sum(na.omit(Genus=='Epithemia'))/length(na.omit(Genus)), #proportion of genus Epithemia
              prop.Eunotia=sum(na.omit(Genus=='Eunotia'))/length(na.omit(Genus)), #proportion of genus Eunotia
              prop.Fragilaria=sum(na.omit(Genus=='Fragilaria'))/length(na.omit(Genus)), #proportion of genus Fragilaria
              prop.Frustulia=sum(na.omit(Genus=='Frustulia'))/length(na.omit(Genus)), #proportion of genus Frustulia
              prop.Gomphonema=sum(na.omit(Genus=='Gomphonema'))/length(na.omit(Genus)), #proportion of genus Gomphonema
              prop.Navicula=sum(na.omit(Genus=='Navicula'))/length(na.omit(Genus)), #proportion of genus Navicula
              prop.Nitzschia=sum(na.omit(Genus=='Nitzschia'))/length(na.omit(Genus)), #proportion of genus Nitzschia
              prop.Rhoicosphenia=sum(na.omit(Genus=='Rhoicosphenia'))/length(na.omit(Genus)), #proportion of genus Rhoicosphenia
              prop.Rhopalodia=sum(na.omit(Genus=='Rhopalodia'))/length(na.omit(Genus)), #proportion of genus Rhopalodia
              prop.Surirella=sum(na.omit(Genus=='Surirella'))/length(na.omit(Genus)), #proportion of genus Surirella
              prop.Synedra=sum(na.omit(Genus=='Synedra'))/length(na.omit(Genus)), #proportion of genus Synedra
              prop.AchOverAchPlusNav=sum(na.omit(Genus=='Achnanthes'))/sum(na.omit(Genus=='Achnanthes'|Genus=='Navicula')), #proportion of genus Achnanthes divided by genus Achnanthes or Navicula
              prop.CymOverCymPlusNav=sum(na.omit(Genus=='Cymbella'))/sum(na.omit(Genus=='Cymbella'|Genus=='Navicula')), #proportion of genus Cymbella divided by genus Cymbella or Navicula
              prop.EpiOverEpiPlusRho=sum(na.omit(Genus=='Epithemia'))/sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #proportion of genus Epithemia divided by genus Epithemia or Rhopalodia
              
              prop.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella'))/length(na.omit(Genus)), #proportion of NNS
              prop.spp.HighMotility=sum(na.omit(Motility=='HM'))/length(na.omit(Motility)), #proportion of high motility - species
              prop.spp.ModMotility=sum(na.omit(Motility=='MM'))/length(na.omit(Motility)), #proportion of moderate motility - species
              prop.spp.MinMotility=sum(na.omit(Motility=='LM'))/length(na.omit(Motility)), #proportion of low motility - species
              prop.spp.NonMotile=sum(na.omit(Motility=='NM'))/length(na.omit(Motility)), #proportion of no motility - species
              
              prop.spp.least.tol=sum(na.omit(designation=='Reference'))/length(na.omit(designation)), #proportion of species that are least tolerant
              prop.spp.most.tol=sum(na.omit(designation=='Stressed'))/length(na.omit(designation)), #proportion of species that are most tolerant
              
              prop.spp.BCG1=sum(na.omit(BCG=='1'))/length(na.omit(BCG)), #proportion of BCG 1 - species
              prop.spp.BCG2=sum(na.omit(BCG=='2'))/length(na.omit(BCG)), #proportion of BCG 2 - species
              prop.spp.BCG3=sum(na.omit(BCG=='3'))/length(na.omit(BCG)), #proportion of BCG 3 - species
              prop.spp.BCG4=sum(na.omit(BCG=='4'))/length(na.omit(BCG)), #proportion of BCG 4 - species
              prop.spp.BCG5=sum(na.omit(BCG=='5'))/length(na.omit(BCG)), #proportion of BCG 5 - species
              prop.spp.BCG6=sum(na.omit(BCG=='6'))/length(na.omit(BCG)), #proportion of BCG 6 - species
              prop.spp.BCG12=sum(na.omit(BCG %in% c('1', "2")))/length(na.omit(BCG)), #proportion of BCG 6 - species
              prop.spp.BCG45=sum(na.omit(BCG %in% c('4', "5")))/length(na.omit(BCG)), #proportion of BCG 6 - species
               
              prop.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high'))/length(na.omit(IndicatorClass_TP)), 
              prop.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low'))/length(na.omit(IndicatorClass_TP)), 
              prop.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low'))/length(na.omit(IndicatorClass_TN)), 
              #prop.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high'))/length(na.omit(IndicatorClass_Cu)), 
              #prop.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high'))/length(na.omit(IndicatorClass_DOC)), 
              #prop.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF'))/length(na.omit(IndicatorClass_Ref)), 
              #prop.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF'))/length(na.omit(IndicatorClass_Ref)), 
              prop.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high'))/length(na.omit(IndicatorClass_TN)), 
              #prop.spp.Heterocy=sum(na.omit(Heterocystous=='yes'))/length(na.omit(Heterocystous)),
              prop.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph'))/length(na.omit(NitrogenUptakeMetabolism2)), 
              prop.spp.Halo=sum(na.omit(Salinity2=='Halo'))/length(na.omit(Salinity2)),
              #prop.spp.ZHR=sum(na.omit(ZHR=='yes'))/length(na.omit(ZHR)), 
              #prop.spp.CRUS=sum(na.omit(CRUS=='yes'))/length(na.omit(CRUS)), 
              #prop.spp.Green=sum(na.omit(Green=='yes'))/length(na.omit(Green)), 
              prop.spp.Planktonic=sum(na.omit(Habitat=='P'))/length(na.omit(Habitat)) 
              
) }


if (sba == T ) {
metrics.sba=ddply(.data = other_combined3.sba, .variables = ~SeqID,.fun = summarise,
                
                # richness 
                #OrgN.NAHON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NAHON')), #count of NAHON - species
                # OrgN.NALON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NALON')), #count of NALON - species
                # OrgN.NHHONF.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF')), #count of NHHONF - species
                # OrgN.NHHONForNHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONF and NHHONO - species
                # OrgN.NHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONO - species
                # OxyReq.DO_100.richness=sum(na.omit(OxygenRequirements=='DO_100')), #count of DO_100 - species
                # OxyReq.DO_100orDO_75.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75')), #count of DO_100 and DO_75 - species
                # OxyReq.DO_100orDO_75orDO_50.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'|OxygenRequirements=='DO_50')), #count of DO_100 and DO_75 and DO_50 - species
                # OxyReq.DO_75.richness=sum(na.omit(OxygenRequirements=='DO_75')), #count of DO_75 - species
                # OxyReq.DO_50.richness=sum(na.omit(OxygenRequirements=='DO_50')), #count of DO_50 - species
                # OxyRed.DO_30.richness=sum(na.omit(OxygenRequirements=='DO_30')), #count of DO_30 - species
                # OxyReq.DO_30andDO_10.richness=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10')), #count of DO_30 and DO_10 - species
                # OxyReq.DO_10.richness=sum(na.omit(OxygenRequirements=='DO_10')), #count of DO_10 - species
                # Salinity.B.richness=sum(na.omit(Salinity=='B')), #count of B - species
                # Salinity.BF.richness=sum(na.omit(Salinity=='BF')), #count of BF - species
                # Salinity.F.richness=sum(na.omit(Salinity=='F')), #count of F - species
                # Salinity.FB.richness=sum(na.omit(Salinity=='FB')), #count of FB - species
                # Saprobic.AM.richness=sum(na.omit(Saprobity=='AM')), #count of AM - species
                # Saprobic.AMorAMPS.richness=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS')), #count of AM and AMPS - species
                # Saprobic.AMPS.richness=sum(na.omit(Saprobity=='AMPS')), #count of AMPS - species
                # Saprobic.BM.richness=sum(na.omit(Saprobity=='BM')), #count of BM - species
                # Saprobic.OS.richness=sum(na.omit(Saprobity=='OS')), #count of OS - species
                # Saprobic.OSandPS.richness=sum(na.omit(Saprobity=='OS'|Saprobity=='PS')), #count of OS and PS - species
                # Saprobic.PS.richness=sum(na.omit(Saprobity=='PS')), #count of PS - species
                # Trophic.E.richness=sum(na.omit(TrophicState=='E')), #count of E - species
                # Trophic.EorI.richness=sum(na.omit(TrophicState=='E'|TrophicState=='I')), #count of E and I - species
                # Trophic.I.richness=sum(na.omit(TrophicState=='I')), #count of I - species
                # Trophic.M.richness=sum(na.omit(TrophicState=='M')), #count of M - species
                # Trophic.ME.richness=sum(na.omit(TrophicState=='ME')), #count of ME - species
                # Trophic.O.richness=sum(na.omit(TrophicState=='O')), #count of O - species
                # Trophic.OorOMorPH.richness=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH')), #count of O, OM, and PH - species
                # Trophic.OM.richness=sum(na.omit(TrophicState=='OM')), #count of OM
                # Trophic.PH.richness=sum(na.omit(TrophicState=='PH')), #count of PH
                # Achnanthes.richness=sum(na.omit(Genus=='Achnanthes')), #count of genus Achnanthes
                # Amphora.richness=sum(na.omit(Genus=='Amphora')), #count of genus Amphora
                # Cocconeis.richness=sum(na.omit(Genus=='Cocconeis')), #count of genus Cocconeis
                # Cyclotella.richness=sum(na.omit(Genus=='Cyclotella')), #count of genus Cyclotella
                # Cymbella.richness=sum(na.omit(Genus=='Cymbella')), #count of genus Cymbella
                # Epithemia.richness=sum(na.omit(Genus=='Epithemia')), #count of genus Epithemia
                # Eunotia.richness=sum(na.omit(Genus=='Eunotia')), #count of genus Eunotia
                # Fragilaria.richness=sum(na.omit(Genus=='Fragilaria')), #count of genus Fragilaria
                # Frustulia.richness=sum(na.omit(Genus=='Frustulia')), #count of genus Frustulia
                # Gomphonema.richness=sum(na.omit(Genus=='Gomphonema')), #count of genus Gomphonema
                # Navicula.richness=sum(na.omit(Genus=='Navicula')), #count of genus Navicula
                # Nitzschia.richness=sum(na.omit(Genus=='Nitzschia')), #count of genus Nitzschia
                # Rhoicosphenia.richness=sum(na.omit(Genus=='Rhoicosphenia')), #count of genus Rhoicosphenia
                # Rhopalodia.richness=sum(na.omit(Genus=='Rhopalodia')), #count of genus Rhopalodia
                # Surirella.richness=sum(na.omit(Genus=='Surirella')), #count of genus Surirella
                # Synedra.richness=sum(na.omit(Genus=='Synedra')), #count of genus Synedra
                # EpiRho.richness=sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #count of genus Epithemia and Rhopalodia
                # 
                # count ----------------------------------------
                #cnt.spp.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella')), #count of NNS
                #cnt.spp.HighMotility=sum(na.omit(Motility=='HM')), #count of high motility - species
                #cnt.spp.ModMotility=sum(na.omit(Motility=='MM')), #count of moderate motility - species
                #cnt.spp.MinMotility=sum(na.omit(Motility=='LM')), #count of low motility - species
                #cnt.spp.NonMotile=sum(na.omit(Motility=='NM')), #count of no motility - species
                cnt.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high')), #from betty
                cnt.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low')), #from betty
                cnt.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low')), #from betty
                cnt.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high')), 
                cnt.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high')), 
                cnt.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF')), 
                cnt.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF')), 
                cnt.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high')), 
                cnt.spp.Heterocy=sum(na.omit(Heterocystous=='yes')), 
                cnt.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph')), 
                cnt.spp.Halo=sum(na.omit(Salinity2=='Halo')), 
                cnt.spp.ZHR=sum(na.omit(ZHR=='yes')), 
                cnt.spp.CRUS=sum(na.omit(CRUS=='yes')), 
                cnt.spp.Green=sum(na.omit(Green=='yes')), 
                cnt.spp.Planktonic=sum(na.omit(Habitat=='P')), 
                cnt.spp.least.tol=sum(na.omit(designation=='Reference')), #cnt of species that are least tolerant
                cnt.spp.most.tol=sum(na.omit(designation=='Stressed')), #cnt of species that are most tolerant
                cnt.spp.BCG1=sum(na.omit(BCG=='1')), #count of BCG 1 - species
                cnt.spp.BCG2=sum(na.omit(BCG=='2')), #count of BCG 2 - species
                cnt.spp.BCG3=sum(na.omit(BCG=='3')), #count of BCG 3 - species
                cnt.spp.BCG4=sum(na.omit(BCG=='4')), #count of BCG 4 - species
                cnt.spp.BCG5=sum(na.omit(BCG=='5')), #count of BCG 5 - species
                cnt.spp.BCG6=sum(na.omit(BCG=='6')), #count of BCG 6 - species
                cnt.spp.BCG12=sum(na.omit(BCG %in% c('1', "2"))), 
                cnt.spp.BCG45=sum(na.omit(BCG %in% c('4', "5"))), 
                
                # proportion -----------------------------------
                # prop.spp.OrgN.NAHON=sum(na.omit(NitrogenUptakeMetabolism=='NAHON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NAHON - species
                # prop.spp.OrgN.NALON=sum(na.omit(NitrogenUptakeMetabolism=='NALON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NALON - species
                # prop.spp.OrgN.NHHONF=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF - species
                # prop.spp.OrgN.NHHONForNHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF and NHHONO - species
                # prop.spp.OrgN.NHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONO - species
                # prop.spp.OxyReq.DO_100=sum(na.omit(OxygenRequirements=='DO_100'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 - species
                # prop.spp.OxyReq.DO_100orDO_75=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 and DO_75 - species
                # prop.spp.OxyReq.DO_75=sum(na.omit(OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_75 - species
                # prop.spp.OxyReq.DO_50=sum(na.omit(OxygenRequirements=='DO_50'))/length(na.omit(OxygenRequirements)), #proportion of DO_50 - species
                # prop.spp.OxyReq.DO_atleast50=sum(na.omit(OxygenRequirements  %in% c('DO_50','DO_75','DO_100')))/length(na.omit(OxygenRequirements)), #proportion of at least DO_50 - species
                # prop.spp.OxyRed.DO_30=sum(na.omit(OxygenRequirements=='DO_30'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 - species
                # prop.spp.OxyReq.DO_30orDO_10=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 and DO_10 - species
                # prop.spp.OxyReq.DO_10=sum(na.omit(OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_10 - species
                # prop.spp.Salinity.B=sum(na.omit(Salinity=='B'))/length(na.omit(Salinity)), #proportion of B - species
                # prop.spp.Salinity.BF=sum(na.omit(Salinity=='BF'))/length(na.omit(Salinity)), #proportion of BF - species
                # prop.spp.Salinity.F=sum(na.omit(Salinity=='F'))/length(na.omit(Salinity)), #proportion of F - species
                # prop.spp.Salinity.FB=sum(na.omit(Salinity=='FB'))/length(na.omit(Salinity)), #proportion of FB - species
                # prop.spp.Saprobic.AM=sum(na.omit(Saprobity=='AM'))/length(na.omit(Saprobity)), #proportion of AM - species
                # prop.spp.Saprobic.AMorAMPS=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AM and AMPS - species
                # prop.spp.Saprobic.AMPS=sum(na.omit(Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AMPS - species
                # prop.spp.Saprobic.BM=sum(na.omit(Saprobity=='BM'))/length(na.omit(Saprobity)), #proportion of BM - species
                # prop.spp.Saprobic.OS=sum(na.omit(Saprobity=='OS'))/length(na.omit(Saprobity)), #proportion of OS - species
                # prop.spp.Saprobic.OSorPS=sum(na.omit(Saprobity=='OS'|Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of OS and PS - species
                # prop.spp.Saprobic.PS=sum(na.omit(Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of PS - species
                # prop.spp.Trophic.E=sum(na.omit(TrophicState=='E'))/length(na.omit(TrophicState)), #proportion of E - species
                # prop.spp.Trophic.EorI=sum(na.omit(TrophicState=='E'|TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of E and I - species
                # prop.spp.Trophic.I=sum(na.omit(TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of I - species
                # prop.spp.Trophic.M=sum(na.omit(TrophicState=='M'))/length(na.omit(TrophicState)), #proportion of M - species
                # prop.spp.Trophic.ME=sum(na.omit(TrophicState=='ME'))/length(na.omit(TrophicState)), #proportion of ME - species
                # prop.spp.Trophic.O=sum(na.omit(TrophicState=='O'))/length(na.omit(TrophicState)), #proportion of O - species
                # prop.spp.Trophic.OorOMorPH=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of O, OM, and PH - species
                # prop.spp.Trophic.OM=sum(na.omit(TrophicState=='OM'))/length(na.omit(TrophicState)), #proportion of OM
                # prop.spp.Trophic.PH=sum(na.omit(TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of PH
                # 
                # prop.Achnanthes=sum(na.omit(Genus=='Achnanthes'))/length(na.omit(Genus)), #proportion of genus Achnanthes
                # prop.Amphora=sum(na.omit(Genus=='Amphora'))/length(na.omit(Genus)), #proportion of genus Amphora
                # prop.Cocconeis=sum(na.omit(Genus=='Cocconeis'))/length(na.omit(Genus)), #proportion of genus Cocconeis
                # prop.Cyclotella=sum(na.omit(Genus=='Cyclotella'))/length(na.omit(Genus)), #proportion of genus Cyclotella
                # prop.Cymbella=sum(na.omit(Genus=='Cymbella'))/length(na.omit(Genus)), #proportion of genus Cymbella
                # prop.Epithemia=sum(na.omit(Genus=='Epithemia'))/length(na.omit(Genus)), #proportion of genus Epithemia
                # prop.Eunotia=sum(na.omit(Genus=='Eunotia'))/length(na.omit(Genus)), #proportion of genus Eunotia
                # prop.Fragilaria=sum(na.omit(Genus=='Fragilaria'))/length(na.omit(Genus)), #proportion of genus Fragilaria
                # prop.Frustulia=sum(na.omit(Genus=='Frustulia'))/length(na.omit(Genus)), #proportion of genus Frustulia
                # prop.Gomphonema=sum(na.omit(Genus=='Gomphonema'))/length(na.omit(Genus)), #proportion of genus Gomphonema
                # prop.Navicula=sum(na.omit(Genus=='Navicula'))/length(na.omit(Genus)), #proportion of genus Navicula
                # prop.Nitzschia=sum(na.omit(Genus=='Nitzschia'))/length(na.omit(Genus)), #proportion of genus Nitzschia
                # prop.Rhoicosphenia=sum(na.omit(Genus=='Rhoicosphenia'))/length(na.omit(Genus)), #proportion of genus Rhoicosphenia
                # prop.Rhopalodia=sum(na.omit(Genus=='Rhopalodia'))/length(na.omit(Genus)), #proportion of genus Rhopalodia
                # prop.Surirella=sum(na.omit(Genus=='Surirella'))/length(na.omit(Genus)), #proportion of genus Surirella
                # prop.Synedra=sum(na.omit(Genus=='Synedra'))/length(na.omit(Genus)), #proportion of genus Synedra
                # prop.AchOverAchPlusNav=sum(na.omit(Genus=='Achnanthes'))/sum(na.omit(Genus=='Achnanthes'|Genus=='Navicula')), #proportion of genus Achnanthes divided by genus Achnanthes or Navicula
                # prop.CymOverCymPlusNav=sum(na.omit(Genus=='Cymbella'))/sum(na.omit(Genus=='Cymbella'|Genus=='Navicula')), #proportion of genus Cymbella divided by genus Cymbella or Navicula
                # prop.EpiOverEpiPlusRho=sum(na.omit(Genus=='Epithemia'))/sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #proportion of genus Epithemia divided by genus Epithemia or Rhopalodia
                # 
                # prop.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella'))/length(na.omit(Genus)), #proportion of NNS
                # prop.spp.HighMotility=sum(na.omit(Motility=='HM'))/length(na.omit(Motility)), #proportion of high motility - species
                # prop.spp.ModMotility=sum(na.omit(Motility=='MM'))/length(na.omit(Motility)), #proportion of moderate motility - species
                # prop.spp.MinMotility=sum(na.omit(Motility=='LM'))/length(na.omit(Motility)), #proportion of low motility - species
                # prop.spp.NonMotile=sum(na.omit(Motility=='NM'))/length(na.omit(Motility)), #proportion of no motility - species
                # 
                prop.spp.least.tol=sum(na.omit(designation=='Reference'))/length(na.omit(designation)), #proportion of species that are least tolerant
                prop.spp.most.tol=sum(na.omit(designation=='Stressed'))/length(na.omit(designation)), #proportion of species that are most tolerant
                
                prop.spp.BCG1=sum(na.omit(BCG=='1'))/length(na.omit(BCG)), #proportion of BCG 1 - species
                prop.spp.BCG2=sum(na.omit(BCG=='2'))/length(na.omit(BCG)), #proportion of BCG 2 - species
                prop.spp.BCG3=sum(na.omit(BCG=='3'))/length(na.omit(BCG)), #proportion of BCG 3 - species
                prop.spp.BCG4=sum(na.omit(BCG=='4'))/length(na.omit(BCG)), #proportion of BCG 4 - species
                prop.spp.BCG5=sum(na.omit(BCG=='5'))/length(na.omit(BCG)), #proportion of BCG 5 - species
                prop.spp.BCG6=sum(na.omit(BCG=='6'))/length(na.omit(BCG)), #proportion of BCG 6 - species
                prop.spp.BCG12=sum(na.omit(BCG %in% c('1', "2")))/length(na.omit(BCG)), #proportion of BCG 6 - species
                prop.spp.BCG45=sum(na.omit(BCG %in% c('4', "5")))/length(na.omit(BCG)), #proportion of BCG 6 - species
                
                prop.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high'))/length(na.omit(IndicatorClass_TP)), 
                prop.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low'))/length(na.omit(IndicatorClass_TP)), 
                prop.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low'))/length(na.omit(IndicatorClass_TN)), 
                prop.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high'))/length(na.omit(IndicatorClass_Cu)), 
                prop.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high'))/length(na.omit(IndicatorClass_DOC)), 
                prop.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF'))/length(na.omit(IndicatorClass_Ref)), 
                prop.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF'))/length(na.omit(IndicatorClass_Ref)), 
                prop.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high'))/length(na.omit(IndicatorClass_TN)), 
                prop.spp.Heterocy=sum(na.omit(Heterocystous=='yes'))/length(na.omit(Heterocystous)),
                #prop.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph'))/length(na.omit(NitrogenUptakeMetabolism2)), 
                #prop.spp.Halo=sum(na.omit(Salinity2=='Halo'))/length(na.omit(Salinity2)),
                prop.spp.ZHR=sum(na.omit(ZHR=='yes'))/length(na.omit(ZHR)), 
                prop.spp.CRUS=sum(na.omit(CRUS=='yes'))/length(na.omit(CRUS)), 
                #prop.spp.Planktonic=sum(na.omit(Habitat=='P'))/length(na.omit(Habitat)) 
                prop.spp.Green=sum(na.omit(Green=='yes'))/length(na.omit(Green)) 
                
) }

if (hybrid == T) { 
metrics.h=ddply(.data = other_combined3, .variables = ~SeqID,.fun = summarise,
                  
                  # richness 
                  OrgN.NAHON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NAHON')), #count of NAHON - species
                  OrgN.NALON.richness=sum(na.omit(NitrogenUptakeMetabolism=='NALON')), #count of NALON - species
                  OrgN.NHHONF.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF')), #count of NHHONF - species
                  OrgN.NHHONForNHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONF and NHHONO - species
                  OrgN.NHHONO.richness=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO')), #count of NHHONO - species
                  OxyReq.DO_100.richness=sum(na.omit(OxygenRequirements=='DO_100')), #count of DO_100 - species
                  OxyReq.DO_100orDO_75.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75')), #count of DO_100 and DO_75 - species
                  OxyReq.DO_100orDO_75orDO_50.richness=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'|OxygenRequirements=='DO_50')), #count of DO_100 and DO_75 and DO_50 - species
                  OxyReq.DO_75.richness=sum(na.omit(OxygenRequirements=='DO_75')), #count of DO_75 - species
                  OxyReq.DO_50.richness=sum(na.omit(OxygenRequirements=='DO_50')), #count of DO_50 - species
                  OxyRed.DO_30.richness=sum(na.omit(OxygenRequirements=='DO_30')), #count of DO_30 - species
                  OxyReq.DO_30andDO_10.richness=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10')), #count of DO_30 and DO_10 - species
                  OxyReq.DO_10.richness=sum(na.omit(OxygenRequirements=='DO_10')), #count of DO_10 - species
                  Salinity.B.richness=sum(na.omit(Salinity=='B')), #count of B - species
                  Salinity.BF.richness=sum(na.omit(Salinity=='BF')), #count of BF - species
                  Salinity.F.richness=sum(na.omit(Salinity=='F')), #count of F - species
                  Salinity.FB.richness=sum(na.omit(Salinity=='FB')), #count of FB - species
                  Saprobic.AM.richness=sum(na.omit(Saprobity=='AM')), #count of AM - species
                  Saprobic.AMorAMPS.richness=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS')), #count of AM and AMPS - species
                  Saprobic.AMPS.richness=sum(na.omit(Saprobity=='AMPS')), #count of AMPS - species
                  Saprobic.BM.richness=sum(na.omit(Saprobity=='BM')), #count of BM - species
                  Saprobic.OS.richness=sum(na.omit(Saprobity=='OS')), #count of OS - species
                  Saprobic.OSandPS.richness=sum(na.omit(Saprobity=='OS'|Saprobity=='PS')), #count of OS and PS - species
                  Saprobic.PS.richness=sum(na.omit(Saprobity=='PS')), #count of PS - species
                  Trophic.E.richness=sum(na.omit(TrophicState=='E')), #count of E - species
                  Trophic.EorI.richness=sum(na.omit(TrophicState=='E'|TrophicState=='I')), #count of E and I - species
                  Trophic.I.richness=sum(na.omit(TrophicState=='I')), #count of I - species
                  Trophic.M.richness=sum(na.omit(TrophicState=='M')), #count of M - species
                  Trophic.ME.richness=sum(na.omit(TrophicState=='ME')), #count of ME - species
                  Trophic.O.richness=sum(na.omit(TrophicState=='O')), #count of O - species
                  Trophic.OorOMorPH.richness=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH')), #count of O, OM, and PH - species
                  Trophic.OM.richness=sum(na.omit(TrophicState=='OM')), #count of OM
                  Trophic.PH.richness=sum(na.omit(TrophicState=='PH')), #count of PH
                  Achnanthes.richness=sum(na.omit(Genus=='Achnanthes')), #count of genus Achnanthes
                  Amphora.richness=sum(na.omit(Genus=='Amphora')), #count of genus Amphora
                  Cocconeis.richness=sum(na.omit(Genus=='Cocconeis')), #count of genus Cocconeis
                  Cyclotella.richness=sum(na.omit(Genus=='Cyclotella')), #count of genus Cyclotella
                  Cymbella.richness=sum(na.omit(Genus=='Cymbella')), #count of genus Cymbella
                  Epithemia.richness=sum(na.omit(Genus=='Epithemia')), #count of genus Epithemia
                  Eunotia.richness=sum(na.omit(Genus=='Eunotia')), #count of genus Eunotia
                  Fragilaria.richness=sum(na.omit(Genus=='Fragilaria')), #count of genus Fragilaria
                  Frustulia.richness=sum(na.omit(Genus=='Frustulia')), #count of genus Frustulia
                  Gomphonema.richness=sum(na.omit(Genus=='Gomphonema')), #count of genus Gomphonema
                  Navicula.richness=sum(na.omit(Genus=='Navicula')), #count of genus Navicula
                  Nitzschia.richness=sum(na.omit(Genus=='Nitzschia')), #count of genus Nitzschia
                  Rhoicosphenia.richness=sum(na.omit(Genus=='Rhoicosphenia')), #count of genus Rhoicosphenia
                  Rhopalodia.richness=sum(na.omit(Genus=='Rhopalodia')), #count of genus Rhopalodia
                  Surirella.richness=sum(na.omit(Genus=='Surirella')), #count of genus Surirella
                  Synedra.richness=sum(na.omit(Genus=='Synedra')), #count of genus Synedra
                  EpiRho.richness=sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #count of genus Epithemia and Rhopalodia
                  
                  # count ----------------------------------------
                  cnt.spp.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella')), #count of NNS
                  cnt.spp.HighMotility=sum(na.omit(Motility=='HM')), #count of high motility - species
                  cnt.spp.ModMotility=sum(na.omit(Motility=='MM')), #count of moderate motility - species
                  cnt.spp.MinMotility=sum(na.omit(Motility=='LM')), #count of low motility - species
                  cnt.spp.NonMotile=sum(na.omit(Motility=='NM')), #count of no motility - species
                  cnt.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high')), #from betty
                  cnt.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low')), #from betty
                  cnt.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low')), #from betty
                  cnt.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high')), 
                  cnt.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high')), 
                  cnt.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF')), 
                  cnt.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF')), 
                  cnt.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high')), 
                  cnt.spp.Heterocy=sum(na.omit(Heterocystous=='yes')), 
                  cnt.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph')), 
                  cnt.spp.Halo=sum(na.omit(Salinity2=='Halo')), 
                  cnt.spp.ZHR=sum(na.omit(ZHR=='yes')), 
                  cnt.spp.CRUS=sum(na.omit(CRUS=='yes')), 
                  cnt.spp.Green=sum(na.omit(Green=='yes')), 
                  cnt.spp.Planktonic=sum(na.omit(Habitat=='P')), 
                  cnt.spp.least.tol=sum(na.omit(designation=='Reference')), #cnt of species that are least tolerant
                  cnt.spp.most.tol=sum(na.omit(designation=='Stressed')), #cnt of species that are most tolerant
                  cnt.spp.BCG1=sum(na.omit(BCG=='1')), #count of BCG 1 - species
                  cnt.spp.BCG2=sum(na.omit(BCG=='2')), #count of BCG 2 - species
                  cnt.spp.BCG3=sum(na.omit(BCG=='3')), #count of BCG 3 - species
                  cnt.spp.BCG4=sum(na.omit(BCG=='4')), #count of BCG 4 - species
                  cnt.spp.BCG5=sum(na.omit(BCG=='5')), #count of BCG 5 - species
                  cnt.spp.BCG6=sum(na.omit(BCG=='6')), #count of BCG 6 - species
                  cnt.spp.BCG12=sum(na.omit(BCG %in% c('1', "2"))), 
                  cnt.spp.BCG45=sum(na.omit(BCG %in% c('4', "5"))), 
                  
                  # proportion -----------------------------------
                  prop.spp.OrgN.NAHON=sum(na.omit(NitrogenUptakeMetabolism=='NAHON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NAHON - species
                  prop.spp.OrgN.NALON=sum(na.omit(NitrogenUptakeMetabolism=='NALON'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NALON - species
                  prop.spp.OrgN.NHHONF=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF - species
                  prop.spp.OrgN.NHHONForNHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONF'|NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONF and NHHONO - species
                  prop.spp.OrgN.NHHONO=sum(na.omit(NitrogenUptakeMetabolism=='NHHONO'))/length(na.omit(NitrogenUptakeMetabolism)), #proportion of NHHONO - species
                  prop.spp.OxyReq.DO_100=sum(na.omit(OxygenRequirements=='DO_100'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 - species
                  prop.spp.OxyReq.DO_100orDO_75=sum(na.omit(OxygenRequirements=='DO_100'|OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_100 and DO_75 - species
                  prop.spp.OxyReq.DO_75=sum(na.omit(OxygenRequirements=='DO_75'))/length(na.omit(OxygenRequirements)), #proportion of DO_75 - species
                  prop.spp.OxyReq.DO_50=sum(na.omit(OxygenRequirements=='DO_50'))/length(na.omit(OxygenRequirements)), #proportion of DO_50 - species
                  prop.spp.OxyReq.DO_atleast50=sum(na.omit(OxygenRequirements %in% c('DO_50','DO_75','DO_100')))/length(na.omit(OxygenRequirements)), #proportion of at least DO_50 - species
                  prop.spp.OxyRed.DO_30=sum(na.omit(OxygenRequirements=='DO_30'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 - species
                  prop.spp.OxyReq.DO_30orDO_10=sum(na.omit(OxygenRequirements=='DO_30'|OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_30 and DO_10 - species
                  prop.spp.OxyReq.DO_10=sum(na.omit(OxygenRequirements=='DO_10'))/length(na.omit(OxygenRequirements)), #proportion of DO_10 - species
                  prop.spp.Salinity.B=sum(na.omit(Salinity=='B'))/length(na.omit(Salinity)), #proportion of B - species
                  prop.spp.Salinity.BF=sum(na.omit(Salinity=='BF'))/length(na.omit(Salinity)), #proportion of BF - species
                  prop.spp.Salinity.F=sum(na.omit(Salinity=='F'))/length(na.omit(Salinity)), #proportion of F - species
                  prop.spp.Salinity.FB=sum(na.omit(Salinity=='FB'))/length(na.omit(Salinity)), #proportion of FB - species
                  prop.spp.Saprobic.AM=sum(na.omit(Saprobity=='AM'))/length(na.omit(Saprobity)), #proportion of AM - species
                  prop.spp.Saprobic.AMorAMPS=sum(na.omit(Saprobity=='AM'|Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AM and AMPS - species
                  prop.spp.Saprobic.AMPS=sum(na.omit(Saprobity=='AMPS'))/length(na.omit(Saprobity)), #proportion of AMPS - species
                  prop.spp.Saprobic.BM=sum(na.omit(Saprobity=='BM'))/length(na.omit(Saprobity)), #proportion of BM - species
                  prop.spp.Saprobic.OS=sum(na.omit(Saprobity=='OS'))/length(na.omit(Saprobity)), #proportion of OS - species
                  prop.spp.Saprobic.OSorPS=sum(na.omit(Saprobity=='OS'|Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of OS and PS - species
                  prop.spp.Saprobic.PS=sum(na.omit(Saprobity=='PS'))/length(na.omit(Saprobity)), #proportion of PS - species
                  prop.spp.Trophic.E=sum(na.omit(TrophicState=='E'))/length(na.omit(TrophicState)), #proportion of E - species
                  prop.spp.Trophic.EorI=sum(na.omit(TrophicState=='E'|TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of E and I - species
                  prop.spp.Trophic.I=sum(na.omit(TrophicState=='I'))/length(na.omit(TrophicState)), #proportion of I - species
                  prop.spp.Trophic.M=sum(na.omit(TrophicState=='M'))/length(na.omit(TrophicState)), #proportion of M - species
                  prop.spp.Trophic.ME=sum(na.omit(TrophicState=='ME'))/length(na.omit(TrophicState)), #proportion of ME - species
                  prop.spp.Trophic.O=sum(na.omit(TrophicState=='O'))/length(na.omit(TrophicState)), #proportion of O - species
                  prop.spp.Trophic.OorOMorPH=sum(na.omit(TrophicState=='O'|TrophicState=='OM'|TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of O, OM, and PH - species
                  prop.spp.Trophic.OM=sum(na.omit(TrophicState=='OM'))/length(na.omit(TrophicState)), #proportion of OM
                  prop.spp.Trophic.PH=sum(na.omit(TrophicState=='PH'))/length(na.omit(TrophicState)), #proportion of PH
                  
                  prop.Achnanthes=sum(na.omit(Genus=='Achnanthes'))/length(na.omit(Genus)), #proportion of genus Achnanthes
                  prop.Amphora=sum(na.omit(Genus=='Amphora'))/length(na.omit(Genus)), #proportion of genus Amphora
                  prop.Cocconeis=sum(na.omit(Genus=='Cocconeis'))/length(na.omit(Genus)), #proportion of genus Cocconeis
                  prop.Cyclotella=sum(na.omit(Genus=='Cyclotella'))/length(na.omit(Genus)), #proportion of genus Cyclotella
                  prop.Cymbella=sum(na.omit(Genus=='Cymbella'))/length(na.omit(Genus)), #proportion of genus Cymbella
                  prop.Epithemia=sum(na.omit(Genus=='Epithemia'))/length(na.omit(Genus)), #proportion of genus Epithemia
                  prop.Eunotia=sum(na.omit(Genus=='Eunotia'))/length(na.omit(Genus)), #proportion of genus Eunotia
                  prop.Fragilaria=sum(na.omit(Genus=='Fragilaria'))/length(na.omit(Genus)), #proportion of genus Fragilaria
                  prop.Frustulia=sum(na.omit(Genus=='Frustulia'))/length(na.omit(Genus)), #proportion of genus Frustulia
                  prop.Gomphonema=sum(na.omit(Genus=='Gomphonema'))/length(na.omit(Genus)), #proportion of genus Gomphonema
                  prop.Navicula=sum(na.omit(Genus=='Navicula'))/length(na.omit(Genus)), #proportion of genus Navicula
                  prop.Nitzschia=sum(na.omit(Genus=='Nitzschia'))/length(na.omit(Genus)), #proportion of genus Nitzschia
                  prop.Rhoicosphenia=sum(na.omit(Genus=='Rhoicosphenia'))/length(na.omit(Genus)), #proportion of genus Rhoicosphenia
                  prop.Rhopalodia=sum(na.omit(Genus=='Rhopalodia'))/length(na.omit(Genus)), #proportion of genus Rhopalodia
                  prop.Surirella=sum(na.omit(Genus=='Surirella'))/length(na.omit(Genus)), #proportion of genus Surirella
                  prop.Synedra=sum(na.omit(Genus=='Synedra'))/length(na.omit(Genus)), #proportion of genus Synedra
                  prop.AchOverAchPlusNav=sum(na.omit(Genus=='Achnanthes'))/sum(na.omit(Genus=='Achnanthes'|Genus=='Navicula')), #proportion of genus Achnanthes divided by genus Achnanthes or Navicula
                  prop.CymOverCymPlusNav=sum(na.omit(Genus=='Cymbella'))/sum(na.omit(Genus=='Cymbella'|Genus=='Navicula')), #proportion of genus Cymbella divided by genus Cymbella or Navicula
                  prop.EpiOverEpiPlusRho=sum(na.omit(Genus=='Epithemia'))/sum(na.omit(Genus=='Epithemia'|Genus=='Rhopalodia')), #proportion of genus Epithemia divided by genus Epithemia or Rhopalodia
                  
                  prop.NNS=sum(na.omit(Genus=='Navicula'|Genus=='Nitzschia'|Genus=='Surirella'))/length(na.omit(Genus)), #proportion of NNS
                  prop.spp.HighMotility=sum(na.omit(Motility=='HM'))/length(na.omit(Motility)), #proportion of high motility - species
                  prop.spp.ModMotility=sum(na.omit(Motility=='MM'))/length(na.omit(Motility)), #proportion of moderate motility - species
                  prop.spp.MinMotility=sum(na.omit(Motility=='LM'))/length(na.omit(Motility)), #proportion of low motility - species
                  prop.spp.NonMotile=sum(na.omit(Motility=='NM'))/length(na.omit(Motility)), #proportion of no motility - species
                  
                  prop.spp.least.tol=sum(na.omit(designation=='Reference'))/length(na.omit(designation)), #proportion of species that are least tolerant
                  prop.spp.most.tol=sum(na.omit(designation=='Stressed'))/length(na.omit(designation)), #proportion of species that are most tolerant
                  
                  prop.spp.BCG1=sum(na.omit(BCG=='1'))/length(na.omit(BCG)), #proportion of BCG 1 - species
                  prop.spp.BCG2=sum(na.omit(BCG=='2'))/length(na.omit(BCG)), #proportion of BCG 2 - species
                  prop.spp.BCG3=sum(na.omit(BCG=='3'))/length(na.omit(BCG)), #proportion of BCG 3 - species
                  prop.spp.BCG4=sum(na.omit(BCG=='4'))/length(na.omit(BCG)), #proportion of BCG 4 - species
                  prop.spp.BCG5=sum(na.omit(BCG=='5'))/length(na.omit(BCG)), #proportion of BCG 5 - species
                  prop.spp.BCG6=sum(na.omit(BCG=='6'))/length(na.omit(BCG)), #proportion of BCG 6 - species
                  prop.spp.BCG12=sum(na.omit(BCG %in% c('1', "2")))/length(na.omit(BCG)), 
                  prop.spp.BCG45=sum(na.omit(BCG %in% c('4', "5")))/length(na.omit(BCG)), 
                  
                  prop.spp.IndicatorClass_TP_high=sum(na.omit(IndicatorClass_TP=='high'))/length(na.omit(IndicatorClass_TP)), 
                  prop.spp.IndicatorClass_TP_low=sum(na.omit(IndicatorClass_TP=='low'))/length(na.omit(IndicatorClass_TP)), 
                  prop.spp.IndicatorClass_TN_low=sum(na.omit(IndicatorClass_TN=='low'))/length(na.omit(IndicatorClass_TN)), 
                  prop.spp.IndicatorClass_Cu_high=sum(na.omit(IndicatorClass_Cu=='high'))/length(na.omit(IndicatorClass_Cu)), 
                  prop.spp.IndicatorClass_DOC_high=sum(na.omit(IndicatorClass_DOC=='high'))/length(na.omit(IndicatorClass_DOC)), 
                  prop.spp.IndicatorClass_Ref=sum(na.omit(IndicatorClass_Ref=='RF'))/length(na.omit(IndicatorClass_Ref)), 
                  prop.spp.IndicatorClass_NonRef=sum(na.omit(IndicatorClass_Ref=='NRF'))/length(na.omit(IndicatorClass_Ref)), 
                  prop.spp.IndicatorClass_TN_high=sum(na.omit(IndicatorClass_TN=='high'))/length(na.omit(IndicatorClass_TN)), 
                  prop.spp.Heterocy=sum(na.omit(Heterocystous=='yes'))/length(na.omit(Heterocystous)),
                  prop.spp.NHeterotroph=sum(na.omit(NitrogenUptakeMetabolism2=='Heterotroph'))/length(na.omit(NitrogenUptakeMetabolism2)), 
                  prop.spp.Halo=sum(na.omit(Salinity2=='Halo'))/length(na.omit(Salinity2)),
                  prop.spp.ZHR=sum(na.omit(ZHR=='yes'))/length(na.omit(ZHR)), 
                  prop.spp.CRUS=sum(na.omit(CRUS=='yes'))/length(na.omit(CRUS)), 
                  prop.spp.Green=sum(na.omit(Green=='yes'))/length(na.omit(Green)), 
                  prop.spp.Planktonic=sum(na.omit(Habitat=='P'))/length(na.omit(Habitat))
                
) }}


if ( diatoms == T ) { write.csv(metrics.d, "R/mmi/diatoms/diatoms.all.metrics.csv", row.names = F)}
if ( sba == T ) { write.csv(metrics.sba, "R/mmi/sba/sba.all.metrics.csv", row.names = F)}
if ( hybrid == T ) { write.csv(metrics.h, "R/mmi/hybrid/hybrid.all.metrics.csv", row.names = F)}

# plots

if(diatoms==T){
final_metrics<-metrics.d
z<-length(colnames(final_metrics))
stations.sub<-subset(stations, stations$SeqID %in% final_metrics$SeqID)
for (i in colnames(final_metrics[,2:z])) { # i<-"prop.Cyclotella"
ggplot(final_metrics, aes(x=stations.sub$SiteStatus, y=final_metrics[,i])) + geom_boxplot() +
  scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) + geom_point(position=position_jitter(height=0))
saveas<-paste0("R/mmi/diatoms/diatoms.plots/",i, ".pdf")
ggsave(saveas, width = 6, height = 6)
}}

if(sba==T){
  final_metrics<-metrics.sba
  z<-length(colnames(final_metrics))
  for (i in colnames(final_metrics[,2:z])) { # i<-"prop.Cyclotella"
    ggplot(final_metrics, aes(x=stations$SiteStatus, y=final_metrics[,i])) + geom_boxplot() +
      scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) + 
      geom_point(position=position_jitter(height=0))
    saveas<-paste0("R/mmi/sba/sba.plots/",i, ".pdf")
    ggsave(saveas, width = 6, height = 6)
  }}


if(hybrid==T){
  final_metrics<-metrics.h
  z<-length(colnames(final_metrics))
  for (i in colnames(final_metrics[,2:z])) { # i<-"prop.Cyclotella"
    ggplot(final_metrics, aes(x=stations$SiteStatus, y=final_metrics[,i])) + geom_boxplot() +
      scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) + 
      geom_point(position=position_jitter(height=0))
    saveas<-paste0("R/mmi/hybrid/hybrid.plots/",i, ".pdf")
    ggsave(saveas, width = 6, height = 6)
  }}

