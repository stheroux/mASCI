# Step 4c Compile winning metrics 


#################################
# Import data ------------------
#################################


combined<-read.csv(paste0("combined.not.scaled.csv"), header=T, row.names=1)


############################ 
# winning pMMIs ###########
############################ 

win.metrics<-read.csv("win.metrics.csv", row.names=1)
win.metrics<-win.metrics$x
  
win.metrics.null<-read.csv("win.metrics.raw.csv", row.names=1)
win.metrics.null<-win.metrics.null$x


# select winners, get scores, and SCALE ---------------------------------
foo<-which(colnames(combined) %in% win.metrics)
combined.win<-combined[,foo] 
combined.win$MMI_raw<-rowMeans(combined.win, na.rm = T)
combined.win.rc<-subset(combined.win, row.names(combined.win) %in% rc$SeqID)
refcalmean<-mean(combined.win.rc$MMI_raw, na.rm = T)
combined.win$MMI_scaled<-combined.win$MMI_raw / refcalmean

foo<-which(colnames(combined) %in% win.metrics.null)
combined.win.null<-combined[,foo] 
combined.win.null$MMI_raw<-rowMeans(combined.win.null, na.rm = T)
combined.win.null.rc<-subset(combined.win.null, row.names(combined.win.null) %in% rc$SeqID)
refcalmean.null<-mean(combined.win.null.rc$MMI_raw, na.rm = T)
combined.win.null$MMI_scaled<-combined.win.null$MMI_raw / refcalmean.null

# output refcal means 
write.csv(refcalmean, paste0(assemblage, "_refcalmean.csv"))

saveas<-paste0(assemblage, ".combined.win.metrics.csv")
write.csv(combined.win, saveas)

saveas<-paste0(assemblage, ".combined.win.metrics.null.csv")
write.csv(combined.win.null, saveas)


foo<-subset(combined.win, row.names(combined.win) %in% rc$SeqID)
zoop<-mean(foo$MMI_scaled, na.rm = T)
 print(paste("RefCalMean equals",zoop))

