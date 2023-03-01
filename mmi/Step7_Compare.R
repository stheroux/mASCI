
# diatoms 

d.scores<-read.csv("diatoms.combined.win.metrics.csv")
names(d.scores)[names(d.scores) == 'X'] <- 'SeqID'

dat<-merge(d.scores, env, by="SeqID", all.x=T, all.y=F, sort=F)

ggplot(dat, aes(x=SiteStatus,y = MMI_scaled)) + geom_boxplot() + 
  scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) + 
  geom_point(position=position_jitter(height=0)) +
  ggsave("diatoms.boxplot.pdf") + ylim(0,1.5)

ggplot(dat, aes(x=log10(TN), y=MMI_scaled)) + geom_point() + stat_smooth(method="glm") + 
  ggsave("diatoms.TN.pdf")
