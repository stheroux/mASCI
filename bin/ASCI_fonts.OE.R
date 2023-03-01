
# ASCI calc fonts 

# StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID

if(length(stations$StationCode) == 0) {stations$StationCode <- stations$stationcode }
if(length(stations$SampleDate) == 0) {stations$SampleDate <- stations$sampledate }
if(length(stations$Replicate) == 0) {stations$Replicate <- stations$replicate }
if(length(stations$PSA6C) == 0) {stations$PSA6C <- stations$psa6c }
if(length(stations$CondQR50) == 0) {stations$CondQR50 <- stations$condqr50 }

if(length(stations$AREA_SQKM) == 0) {stations$AREA_SQKM <- stations$area_sqkm }
if(length(stations$KFCT_AVE) == 0) {stations$KFCT_AVE <- stations$kfct_ave }
if(length(stations$PPT_00_09) == 0) {stations$PPT_00_09 <- stations$ppt_00_09 }
if(length(stations$TMAX_WS) == 0) {stations$TMAX_WS <- stations$tmax_ws }
if(length(stations$MAX_ELEV) == 0) {stations$MAX_ELEV <- stations$max_elev }
if(length(stations$SITE_ELEV) == 0) {stations$SITE_ELEV <- stations$site_elev }
if(length(stations$AtmCa) == 0) {stations$AtmCa <- stations$atmca }

if(length(stations$CaO_Mean) == 0) {stations$CaO_Mean <- stations$cao_mean }
if(length(stations$MgO_Mean) == 0) {stations$MgO_Mean <- stations$mgo_mean }
if(length(stations$S_Mean) == 0) {stations$S_Mean <- stations$s_mean }
if(length(stations$UCS_Mean) == 0) {stations$UCS_Mean <- stations$ucs_mean }
if(length(stations$LPREM_mean) == 0) {stations$LPREM_mean <- stations$lprem_mean }
if(length(stations$AtmMg) == 0) {stations$AtmMg <- stations$atmmg }
if(length(stations$AtmSO4) == 0) {stations$AtmSO4 <- stations$atmso4 }
if(length(stations$MINP_WS) == 0) {stations$MINP_WS <- stations$minp_ws }
if(length(stations$MEANP_WS) == 0) {stations$MEANP_WS <- stations$meanp_ws }
if(length(stations$SumAve_P) == 0) {stations$SumAve_P <- stations$sumave_p }
if(length(stations$TMAX_WS) == 0) {stations$TMAX_WS <- stations$tmax_ws }
if(length(stations$XWD_WS) == 0) {stations$XWD_WS <- stations$xwd_ws }
if(length(stations$MAXWD_WS) == 0) {stations$MAXWD_WS <- stations$maxwd_ws }
if(length(stations$BDH_AVE) == 0) {stations$BDH_AVE <- stations$bdh_ave }
if(length(stations$KFCT_AVE) == 0) {stations$KFCT_AVE <- stations$kfct_ave }
if(length(stations$PRMH_AVE) == 0) {stations$PRMH_AVE <- stations$prmh_ave }


if(length(stations$LST32AVE) == 0) {stations$LST32AVE <- stations$lst32ave }


# 
# # tax
# if(length(tax$StationCode) == 0) {tax$StationCode <- tax$stationcode }
# if(length(tax$SampleDate) == 0) {tax$SampleDate <- tax$sampledate }
# if(length(tax$Replicate) == 0) {tax$Replicate <- tax$replicate }
# if(length(tax$SampleTypeCode) == 0) {tax$SampleTypeCode <- tax$sampletypecode }
# if(length(tax$BAResult) == 0) {tax$BAResult <- tax$baresult }
# if(length(tax$Result) == 0) {tax$Result <- tax$result }
# if(length(tax$FinalID) == 0) {tax$FinalID <- tax$finalid }
