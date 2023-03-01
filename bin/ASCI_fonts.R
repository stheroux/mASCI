
# ASCI calc fonts 

# StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID

if(length(sta$StationCode) == 0) {sta$StationCode <- sta$stationcode }
if(length(sta$SampleDate) == 0) {sta$SampleDate <- sta$sampledate }
if(length(sta$Replicate) == 0) {sta$Replicate <- sta$replicate }
if(length(sta$PSA6C) == 0) {sta$PSA6C <- sta$psa6c }
if(length(sta$CondQR50) == 0) {sta$CondQR50 <- sta$condqr50 }

if(length(sta$AREA_SQKM) == 0) {sta$AREA_SQKM <- sta$area_sqkm }
if(length(sta$KFCT_AVE) == 0) {sta$KFCT_AVE <- sta$kfct_ave }
if(length(sta$PPT_00_09) == 0) {sta$PPT_00_09 <- sta$ppt_00_09 }
if(length(sta$TMAX_WS) == 0) {sta$TMAX_WS <- sta$tmax_ws }
if(length(sta$MAX_ELEV) == 0) {sta$MAX_ELEV <- sta$max_elev }
if(length(sta$SITE_ELEV) == 0) {sta$SITE_ELEV <- sta$site_elev }
if(length(sta$AtmCa) == 0) {sta$AtmCa <- sta$atmca }

if(length(sta$CaO_Mean) == 0) {sta$CaO_Mean <- sta$cao_mean }
if(length(sta$MgO_Mean) == 0) {sta$MgO_Mean <- sta$mgo_mean }
if(length(sta$S_Mean) == 0) {sta$S_Mean <- sta$s_mean }
if(length(sta$UCS_Mean) == 0) {sta$UCS_Mean <- sta$ucs_mean }
if(length(sta$LPREM_mean) == 0) {sta$LPREM_mean <- sta$lprem_mean }
if(length(sta$AtmMg) == 0) {sta$AtmMg <- sta$atmmg }
if(length(sta$AtmSO4) == 0) {sta$AtmSO4 <- sta$atmso4 }
if(length(sta$MINP_WS) == 0) {sta$MINP_WS <- sta$minp_ws }
if(length(sta$MEANP_WS) == 0) {sta$MEANP_WS <- sta$meanp_ws }
if(length(sta$SumAve_P) == 0) {sta$SumAve_P <- sta$sumave_p }
if(length(sta$TMAX_WS) == 0) {sta$TMAX_WS <- sta$tmax_ws }
if(length(sta$XWD_WS) == 0) {sta$XWD_WS <- sta$xwd_ws }
if(length(sta$MAXWD_WS) == 0) {sta$MAXWD_WS <- sta$maxwd_ws }
if(length(sta$BDH_AVE) == 0) {sta$BDH_AVE <- sta$bdh_ave }
if(length(sta$KFCT_AVE) == 0) {sta$KFCT_AVE <- sta$kfct_ave }
if(length(sta$PRMH_AVE) == 0) {sta$PRMH_AVE <- sta$prmh_ave }


if(length(sta$LST32AVE) == 0) {sta$LST32AVE <- sta$lst32ave }



# tax
if(length(tax$StationCode) == 0) {tax$StationCode <- tax$stationcode }
if(length(tax$SampleDate) == 0) {tax$SampleDate <- tax$sampledate }
if(length(tax$Replicate) == 0) {tax$Replicate <- tax$replicate }
if(length(tax$SampleTypeCode) == 0) {tax$SampleTypeCode <- tax$sampletypecode }
if(length(tax$BAResult) == 0) {tax$BAResult <- tax$baresult }
if(length(tax$Result) == 0) {tax$Result <- tax$result }
if(length(tax$FinalID) == 0) {tax$FinalID <- tax$finalid }
