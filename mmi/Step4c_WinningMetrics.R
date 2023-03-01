#find best custom combination of metrics for pMMI 
#step 1:  import your data from Step4b so it is in your global environment 


### Winning metrics --------------------------------------------------------------

if (assemblage == "diatoms") {
  if (ASV==F) {
  win.metrics<-c("OxyReq.DO_10.richness",
  "Saprobic.AM.richness",
  "prop.spp.OrgN.NHHONF_raw",
  "Salinity.F.richness_raw",
  "cnt.spp.IndicatorClass_TN_low_raw",
  "prop.spp.BCG3_raw",
               NULL ) 
  write.csv(win.metrics, "win.metrics.csv")

  win.metrics.raw<-c("OxyReq.DO_10.richness_raw",
                     "Saprobic.AM.richness_raw",
                     "prop.spp.OrgN.NHHONF_raw",
                     "Salinity.F.richness_raw",
                     "cnt.spp.IndicatorClass_TN_low_raw",
                     "prop.spp.BCG3_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }


if (ASV==T) {
  win.metrics<-c("trait_afdm.2.prop.spp",
                 "trait_tn.2.richness",
                 "trait_tp.2.prop.spp_raw",
                 "trait_doc.2.richness_raw",
                 NULL ) 
  write.csv(win.metrics, "win.metrics.csv")
  
  win.metrics.raw<-c("trait_afdm.2.prop.spp_raw",
                     "trait_tn.2.richness_raw",
                     "trait_tp.2.prop.spp_raw",
                     "trait_doc.2.richness_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }
  
  }





if (assemblage == "sba") { 
  
  win.metrics<-c("prop.spp.IndicatorClass_TP_high_raw",
                 "prop.spp.IndicatorClass_DOC_high_raw",
                 "prop.spp.IndicatorClass_NonRef_raw",
                 "prop.spp.ZHR_raw",
                 NULL ) 
  write.csv(win.metrics, "win.metrics.csv")
  
  win.metrics.raw<-c("prop.spp.IndicatorClass_TP_high_raw",
                     "prop.spp.IndicatorClass_DOC_high_raw",
                     "prop.spp.IndicatorClass_NonRef_raw",
                     "prop.spp.ZHR_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }

if (assemblage == "hybrid") { 
  
  win.metrics<-c("cnt.spp.most.tol",
                  "prop.spp.ZHR_raw",
                 "Salinity.BF.richness",
                 "prop.spp.Planktonic",
                 "cnt.spp.IndicatorClass_TP_high",
                 "prop.spp.Trophic.E",
                 "EpiRho.richness",
                 "OxyRed.DO_30.richness",
                 NULL ) 
  write.csv(win.metrics, "win.metrics.csv")
  
  win.metrics.raw<-c(
                    "cnt.spp.most.tol_raw",
                      "prop.spp.ZHR_raw",
                     "Salinity.BF.richness_raw",
                     "prop.spp.Planktonic_raw",
                     "cnt.spp.IndicatorClass_TP_high_raw",
                     "prop.spp.Trophic.E_raw",
                     "EpiRho.richness_raw",
                     "OxyRed.DO_30.richness_raw",
                     NULL ) 
  write.csv(win.metrics.raw, "win.metrics.raw.csv") }


