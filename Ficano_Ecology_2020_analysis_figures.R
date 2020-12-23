
#librarys#####
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("car")
library("lme4")

#data files####
amf <- read.csv(("LS_AMF_data_REAL.csv"),header=TRUE)
dat <- read.csv(("LS_AMF_SNFdata_REAL.csv"),header=TRUE)
#merging AMF and SNF data
data<-merge(amf,dat,by.x="sample", by.y="ï..P.num")  

#cleaning data files ####
data[data$percent_colonization=="N/A",]$percent_coloinzation <- NA
data$S.num = as.factor(data$especies)
data$nitro = as.factor(data$nitro)
data$views = as.integer(data$views)
data$tot.nod.mass = as.numeric(data$tot.nod.mass)
data = data %>%
  filter(views != 0)
data<-data%>%
  filter(!is.na(especies))

#Adding variables ####

data$tot.mass<-data$l.mass+data$s.mass+data$r.mass
data$shootroot<-(data$s.mass+data$l.mass)/data$r.mass
data$rootplant<-data$r.mass/data$tot.mass
data$abv = data$s.mass + data$l.mass
data$nod.al = data$tot.nod.mass/data$tot.mass
data$nod.alr = data$tot.nod.mass/data$r.mass
data$X.col.per.g<-data$X.col/data$tot.mass


#testing for normaility of AMF data####
hist(data$X.col)
shapiro.test(data$X.col)
#testing for normaility of biomass data
hist(data$tot.mass)
shapiro.test(data$tot.mass)
#testing for normaility of nodule biomass data
hist(data$tot.nod.mass)
shapiro.test(data$tot.nod.mass)


#Subsetting data #### 

#by species: data.p = Pentaclethra macroloba, data.z = Zygia longifolia, data.s = Stryphnodendrom microstachyum
data.p = data %>%
  filter(S.num == 1)
data.z = data %>%
  filter(S.num == 2)
data.s = data %>%
  filter(S.num == 3)
data.notp = data%>%
  filter(S.num != 1)

#by nitrogen treatment: data.n = added N, data.non = no added n
data.n = data %>%
  filter(nitro == 1)
data.non = data %>%
  filter(nitro == 0)


#AMF versus gap/understory########################################

#testing for significant differences in AMF root colonization ("X.col") between the gap and understory ("type')
wilcox.test(X.col~type, data=data)

#by species
wilcox.test(X.col~type, data=data.p)

wilcox.test(X.col~type, data=data.z)

wilcox.test(X.col~type, data=data.s)


#controlling for seedling biomass ("tot.mass")
wilcox.test((X.col/tot.mass)~type, data=data)

#by species
wilcox.test((X.col/tot.mass)~type, data=data.p)

wilcox.test((X.col/tot.mass)~type, data=data.z)

wilcox.test((X.col/tot.mass)~type, data=data.s)



#Figure 1: AMF in canopy gaps versus intact understory (controlling for seedling biomass, and not)
AMF.gap.understory<-ggplot(data,aes(x=type,y=X.col))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=18))+
  xlab("")

AMF.per.g.gap.understory<-ggplot(data,aes(x=type,y=X.col.per.g))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name=expression('AMF colonization g'^'-1'*' seedling'))+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=18))+
  xlab("")

ggarrange(AMF.gap.understory,AMF.per.g.gap.understory,
          ncol=2,nrow=1, labels=c("a","b"),
          font.label=list(size=20))

#AMF versus light and nitrogen treatments ########################

#generalized linear model: AMF root colonization ("X.col") ~ species ("S.num"), light ("light.per"), N treatment ("nitro"),
#seedling biomass ("tot.mass"), environment ("type")
m1<-glm(X.col~S.num+light.per*nitro+tot.mass+type,data=data, family=binomial(),weights=views)
summary(m1)

#by species
#light * N interaction removed if not significant

m1.p<-glm(X.col~light.per*nitro+tot.mass+type,data=data.p, family=binomial(), weights=views)
summary(m1.p)

m1.z<-glm(X.col~light.per+nitro+tot.mass+type,data=data.z, family=binomial(), weights=views)
summary(m1.z)

m1.s<-glm(X.col~light.per+tot.mass+type,data=data.s.non, family=binomial(), weights=views)
summary(m1.s)

#Figure 2: AMF versus light and nitrogen treatments
ggplot(data, aes(x=light.per, y=X.col, weight=views))+
  geom_point(aes(shape=nitro),size=3,alpha=0.7)+
  scale_shape_manual(values=c(1, 16))+
  xlab("")+
  theme_pubr()+
  scale_y_continuous(name="AMF root colonization (%)", labels=scales::percent)+
  theme(axis.title.y=element_text(size=18,angle=90),
        axis.text=element_text(size=14))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  theme(legend.position="none")+
  theme(plot.margin=unit(c(.1,.1,2.1,.1),"cm"))+
  coord_cartesian(ylim=c(.25,1),clip="off")+
  annotate("segment",x=.05,xend=.10, y=.14,yend=.14,color="grey30",size=3)+
  annotate("segment",x=.09,xend=.20, y=.13,yend=.13,color="grey45",size=3)+
  annotate("segment",x=.25,xend=.46, y=.13,yend=.13,color="grey70",size=3)+
  annotate("text",x=.25,y=0.04,label="Light \n(% total transmittance)",size=5.75)+
  geom_smooth(data=data.n,aes(x=light.per,y=X.col),color="black",size=1.75,method="glm",
              method.args=list(family="binomial"),se=FALSE)+
  geom_smooth(data=data.non,aes(x=light.per,y=X.col),color="black",size=1.75,linetype=2,method="glm",
              method.args=list(family="binomial"),se=FALSE)

#Evidence of Tripartite Symbiosis ################################

#generalized linear model: seedling biomass ("tot.mass") ~ species ("S.num"), %AMF ("X.col"), nodule biomass ("tot.nod.mass"),
# nitrogen treatment ("nitro"), light treatment ("light.per")
m2<-glm(tot.mass~S.num+nitro+light.per+X.col+tot.nod.mass+type,data=data, family= gaussian(link="log"))
summary(m2)

#by species
m2.p<-glm(tot.mass~light.per+nitro+X.col+tot.nod.mass+type,data=data.p, family= gaussian(link="log"))
summary(m2.p)

m2.z<-glm(tot.mass~light.per+nitro+X.col+tot.Nfixd+type,data=data.z, family= gaussian(link="log"))
summary(m2.z)

m2.s<-glm(tot.mass~light.per+nitro+X.col+tot.Nfixd+type,data=data.s, family= gaussian(link="log"))
summary(m2.s)

#Spearman correlation test between %AMF colonization and nodule biomass
cor.test(data$tot.nod.mass, data$X.col, method = "spearman")

#by species
cor.test(data.p$tot.nod.mass, data.p$X.col, method = "spearman")
cor.test(data.z$tot.nod.mass, data.z$X.col, method = "spearman")
cor.test(data.s$tot.nod.mass, data.s$X.col, method = "spearman")

#Figure 3a: AMF & nodule biomass
AMF.SNF<-ggplot(data, aes(x=tot.nod.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=.7)+
  theme_pubr()+
  stat_cor(data=data, aes(x=tot.nod.mass, y=X.col,group = 1),method = "spearman",cor.coef.name="rho",label.x.npc = "centre",
           label.y.npc = "bottom")+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Nodule biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data,aes(x=tot.nod.mass,y=X.col),method="lm",formula=y~log(x+1),size=1.75,
              se=FALSE,color="black")

#Spearman correlation test between Seedling biomass and Nodule biomass
cor.test(data$tot.nod.mass, data$tot.mass, method = "spearman")

#by species
cor.test(data.p$tot.nod.mass, data.p$X.col, method = "spearman")
cor.test(data.z$tot.nod.mass, data.z$X.col, method = "spearman")
cor.test(data.s$tot.nod.mass, data.s$X.col, method = "spearman")

#Figure 3b: Seedlng biomass versus SNF
SNF.Biomass<-ggplot(data,aes(x=tot.mass,y=tot.nod.mass,color=S.num))+
  geom_point(size=2.5,alpha=.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="Nodule biomass (g)\n")+
  xlab("")+
  coord_cartesian(ylim=c(0,.6),clip="on")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data,aes(x=tot.mass,y=tot.nod.mass),method="lm",size=1.75,formula=log(y+1)~log(x),
              se=FALSE,color="black")

#Spearman correlation test between Seedling biomass and %AMF
cor.test(data$X.col, data$tot.mass, method = "spearman")

#by species
cor.test(data.p$X.col, data.p$X.col, method = "spearman")
cor.test(data.z$X.col, data.z$X.col, method = "spearman")
cor.test(data.s$X.col, data.s$X.col, method = "spearman")

#Figure 3c: Seedling biomass versus AMF
AMF.Biomass<-ggplot(data,aes(x=tot.mass,y=X.col,color=S.num))+
  geom_point(size=2.5, alpha=.7)+
  theme_pubr()+
  stat_cor(data=data, aes(x=tot.mass, y=X.col,group = 1),method = "spearman",cor.coef.name="rho",label.x.npc = "centre",
           label.y.npc = "bottom")+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Seedling biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data,aes(x=tot.mass,y=X.col),method="lm",size=1.75,formula=y~log(x),
              se=FALSE,color="black")

#Figure 3: Tripartite symbiosis
ggarrange(AMF.SNF,
          ggarrange(SNF.Biomass,AMF.Biomass, nrow=2, labels = c("b","c"),font.label=list(size=20),hjust=1,vjust=1),
          ncol=2,labels = "a", vjust=1,
          font.label=list(size=20))

#Biomass Allocation tradeoffs ####################################

#Spearman correlation test for above to belowground biomass ratio ("shootroot") and %AMF ("X.col) 
cor.test(data$shootroot,data$X.col, method="spearman")
#by species
cor.test(data.p$shootroot,data.p$X.col, method="spearman")
cor.test(data.z$shootroot,data.z$X.col, method="spearman")
cor.test(data.s$shootroot,data.s$X.col, method="spearman")


#Spearman correaltion test for above to belowground biomass ratio ("shootroot") and Nodule biomass ("tot.nod.mass")
cor.test(data$shootroot,data$tot.nod.mass, method="spearman")
#by species
cor.test(data.p$shootroot,data.p$tot.nod.mass, method="spearman")
cor.test(data.z$shootroot,data.z$tot.nod.mass, method="spearman")
cor.test(data.s$shootroot,data.s$tot.nod.mass, method="spearman")


#Spearman correaltion test for above to belowground biomass ratio ("shootroot") and seedling biomass ("tot.mass")
cor.test(data$shootroot,data$tot.mass, method="spearman")
#by species
cor.test(data.p$shootroot,data.p$tot.mass, method="spearman")
cor.test(data.z$shootroot,data.z$tot.mass, method="spearman")
cor.test(data.s$shootroot,data.s$tot.mass, method="spearman")


#SUPPLEMENTAL FIGURES / ANALYSIS####

#separating data by species and treatment for creation of supplemental figures

#P. macroloba
#nitrgoen treatment
data.p.n = data.p %>%
  filter(nitro == 1)
data.p.non = data.p %>%
  filter(nitro == 0)

#Z. longifolia
#nitrogen treatment
data.z.n = data.z %>%
  filter(nitro == 1)
data.z.non = data.z %>%
  filter(nitro == 0)

#S. microstachyum
#nitrogen treatment
data.s.n = data.s %>%
  filter(nitro == 1)
data.s.non = data.s %>%
  filter(nitro == 0)


# Supplemental Figure 1 (Fig. 1 by species)
p.AMF.gap.understory<-ggplot(data.p,aes(x=type,y=X.col))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="",labels=scales::percent)+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank())+
  xlab("")

p.AMF.per.g.gap.understory<-ggplot(data.p,aes(x=type,y=X.col.per.g))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="")+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank())+
  xlab("")

z.AMF.gap.understory<-ggplot(data.z,aes(x=type,y=X.col))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="AMF root colonization (%)\n  ",labels=scales::percent)+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank())+
  xlab("")

z.AMF.per.g.gap.understory<-ggplot(data.z,aes(x=type,y=X.col.per.g))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name=" \n  ")+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank())+
  xlab("")

s.AMF.gap.understory<-ggplot(data.s,aes(x=type,y=X.col))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="",labels=scales::percent)+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=18))+
  xlab("")

s.AMF.per.g.gap.understory<-ggplot(data.s,aes(x=type,y=X.col.per.g))+
  geom_boxplot(fill=c("grey70","grey30"))+
  theme_pubr()+
  scale_y_continuous(name="")+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=18))+
  xlab("")


Fig.1.species1<-ggarrange(p.AMF.gap.understory,z.AMF.gap.understory,s.AMF.gap.understory,
                          labels=c("a","c","e"),
                          ncol=1,nrow=3,
                          hjust=c(-3,-3,-3),
                          vjust=c(1.5,1.5,1.5),
                          align="hv",
                          font.label=list(size=20))

Fig.1.species2<-ggarrange(p.AMF.per.g.gap.understory,z.AMF.per.g.gap.understory,s.AMF.per.g.gap.understory,
                          labels=c("b","d","f"),
                          ncol=1,nrow=3,
                          hjust=c(-3,-3,-8),
                          vjust=c(1.5,1.5,1.5),
                          align="hv",
                          font.label=list(size=20))

ggarrange(Fig.1.species1,NULL,Fig.1.species2,
          ncol=3,nrow=1,
          widths=c(1,0.15,1))

#Supplemental Figure 2 (Fig. 2 by species)
fig2.p<-ggplot(data.p, aes(x=light.per, y=X.col, weight=views))+
  geom_point(aes(shape=nitro),size=2.5,alpha=0.7)+
  scale_shape_manual(values=c(1, 16))+
  xlab("")+
  theme_pubr()+
  scale_y_continuous(name=" \n", labels=scales::percent)+
  theme(axis.title.y=element_text(size=18,angle=90),
        axis.text=element_text(size=14))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4))+
  theme(legend.position="none")+
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm"))+
  coord_cartesian(ylim=c(.25,1),clip="off")+
  geom_smooth(data=data.p.n,aes(x=light.per,y=X.col),color="black",size=1.5,method="glm",
              method.args=list(family="binomial"),se=FALSE)+
  geom_smooth(data=data.p.non,aes(x=light.per,y=X.col),color="black",size=1.75,linetype=2,method="glm",
              method.args=list(family="binomial"),se=FALSE)

fig2.z<-ggplot(data.z, aes(x=light.per, y=X.col, weight=views))+
  geom_point(aes(shape=nitro),size=2.5,alpha=0.7)+
  scale_shape_manual(values=c(1, 16))+
  xlab("")+
  theme_pubr()+
  scale_y_continuous(name="AMF root colonization (%)\n", labels=scales::percent)+
  theme(axis.title.y=element_text(size=18,angle=90),
        axis.text=element_text(size=14))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4))+
  theme(legend.position="none")+
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm"))+
  coord_cartesian(ylim=c(.25,1),clip="off")+
  geom_smooth(data=data.z,aes(x=light.per,y=X.col),color="black",size=1.75,method="glm",
              method.args=list(family="binomial"),se=FALSE)

fig2.s<-ggplot(data.s, aes(x=light.per, y=X.col, weight=views))+
  geom_point(aes(shape=nitro),size=2.5,alpha=0.7)+
  scale_shape_manual(values=c(1, 16))+
  xlab("")+
  theme_pubr()+
  scale_y_continuous(name=" \n", labels=scales::percent)+
  theme(axis.title.y=element_text(size=18,angle=90),
        axis.text=element_text(size=14))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  theme(legend.position="none")+
  theme(plot.margin=unit(c(.1,.1,.1,.1),"cm"))+
  coord_cartesian(ylim=c(.25,1),clip="off")+
  annotate("segment",x=.05,xend=.10, y=.12,yend=.12,color="grey30",size=3)+
  annotate("segment",x=.09,xend=.20, y=.11,yend=.11,color="grey45",size=3)+
  annotate("segment",x=.25,xend=.46, y=.11,yend=.11,color="grey70",size=3)+
  annotate("text",x=.25,y=0,label="Light \n(% total transmittance)",size=5.75)+
  geom_smooth(data=data.s.n,aes(x=light.per,y=X.col),color="black",size=1.75,method="glm",
              method.args=list(family="binomial"),se=FALSE)+
  geom_smooth(data=data.s.non,aes(x=light.per,y=X.col),color="black",size=1.5,linetype=2,method="glm",
              method.args=list(family="binomial"),se=FALSE)

ggarrange(fig2.p,fig2.z,fig2.s,NULL,
          ncol=1,nrow=4,
          labels=c("a","b","c",""),hjust=-2,vjust=.7,
          heights=c(1,1,1,.2),
          font.label = list(size=20))

#Supplemental Figure 3 (fig. 3 for P. mac)
AMF.SNF.p.mac<-ggplot(data.p, aes(x=tot.nod.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Nodule biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))

AMF.biomass.p.mac<-ggplot(data.p,aes(x=tot.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Seedling biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))

SNF.biomass.p.mac<-ggplot(data.p,aes(x=tot.mass,y=tot.nod.mass,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="Nodule biomass (g)\n  ")+
  xlab("")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.p,aes(x=tot.mass,y=tot.nod.mass),method="lm",size=1.75,formula=log(y+1)~log(x),
              se=FALSE,color="black")+
  coord_cartesian(ylim=c(0,.6),clip="on")

ggarrange(AMF.SNF.p.mac,
          ggarrange(SNF.biomass.p.mac,AMF.biomass.p.mac, nrow=2, labels = c("b","c"),hjust=1,vjust=1,font.label = list(size=20)),
          ncol=2,labels = "a", vjust=1,font.label = list(size=20))

#Supplemental Figure 4 (Fig. 3 for Z. longifolia)

AMF.SNF.z.long<-ggplot(data.z, aes(x=tot.nod.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Nodule biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))

AMF.biomass.z.long<-ggplot(data.z,aes(x=tot.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Seedling biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.z,aes(x=tot.mass,y=X.col),method="lm",size=1.75,formula=y~log(x),
              se=FALSE,color="black")

SNF.biomass.z.long<-ggplot(data.z,aes(x=tot.mass,y=tot.nod.mass,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="Nodule biomass (g)")+
  xlab("")+
  coord_cartesian(ylim=c(0,.1),xlim=c(0,6.5),clip="on")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.z,aes(x=tot.mass,y=tot.nod.mass),method="lm",size=1.75,formula=log(y+1)~log(x),
              se=FALSE,color="black")

ggarrange(AMF.SNF.z.long,
          ggarrange(SNF.biomass.z.long,AMF.biomass.z.long, nrow=2, labels = c("b","c"),hjust=1,vjust=1,font.label = list(size=20)),
          ncol=2,labels = "a", vjust=1,font.label = list(size=20))

#Supplemental Figure 5 (Fig. 3 for S. microstachyum)

AMF.SNF.s.micro<-ggplot(data.s, aes(x=tot.nod.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Nodule biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.s,aes(x=tot.nod.mass,y=X.col),method="lm",formula=y~log(x+1),size=1.75,
              se=FALSE,color="black")

AMF.biomass.s.micro<-ggplot(data.s,aes(x=tot.mass,y=X.col,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="AMF root colonization (%)",labels=scales::percent)+
  scale_x_continuous(name="Seedling biomass (g)")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.s,aes(x=tot.mass,y=X.col),method="lm",size=1.75,formula=y~log(x),
              se=FALSE,color="black")

SNF.biomass.s.micro<-ggplot(data.s,aes(x=tot.mass,y=tot.nod.mass,color=S.num))+
  geom_point(size=2.5,alpha=0.7)+
  theme_pubr()+
  scale_color_manual(values=c("1"="navyblue","2"="skyblue4","3"="skyblue"))+
  scale_y_continuous(name="Nodule biomass (g)")+
  xlab("")+
  coord_cartesian(ylim=c(0,.12),xlim=c(0,11),clip="on")+
  theme(legend.position="none")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text=element_text(size=14))+
  geom_smooth(data=data.s,aes(x=tot.mass,y=tot.nod.mass),method="lm",size=1.75,formula=log(y+1)~log(x),
              se=FALSE,color="black")

ggarrange(AMF.SNF.s.micro,
          ggarrange(SNF.biomass.s.micro,AMF.biomass.s.micro, nrow=2, labels = c("b","c"),hjust=1,vjust=1,font.label = list(size=20)),
          ncol=2,labels = "a", vjust=1,font.label = list(size=20))

#Calculating N fixed to compare glm results with nodule biomass glm results. Results found in Table S1.
#Calculating tot.Nfixd (Total N fixed per plant/hr using isotope data)####

#Changing total N in isotope samples from micrograms to miligrams 

data$Total.N.e <- data$total.N.15n*.001

data$Total.N.c <- data$total.N.c*.001

#Calculating %N for each sample from the total N and sample weight of enriched samples #

data$XN.e<-((data$Total.N.e)/data$mg.tin.e)*100

#Calculating %N for each sample from the total N and sample weight of control samples #

data$XN.c<-((data$total.N.c)/data$mg.tin.c)*100

#Syringe Volume replaced (proportion), Enrichment of normal atmosphere, enrichment of added gas #

svr<-.2

X15Nair<-.003663

X15Ngas<-.99

##svr.c 

svr.c <- 0

#Calculating the expected enrichment of syringe headspace during incubation.

#this will be very close to the volume replaced times the %15Ngas (how close depends on how close the gas enrichment is to 100%)

syr.enrich<-(X15Ngas*svr)+(X15Nair*(1-svr))

syr.enrich.c <-(X15Ngas*svr.c)+(X15Nair*(1-svr.c))

#Calculating the Percent of N atoms fixed in each sample #

#calculated by taking the difference between the enriched X15N and control X15N and (?dividing by expected enrichment)

data$XNfixdX<-with(data, (X15n.atX.15n-X15n.atX.c)/syr.enrich)

#Some samples only ended up having values for controls because enriched samples were below detection limit, 

#to avoid negatives I made them zeros

data$XNfixdX[data$XNfixdX < 0] = 0 

#Calculating the total amount of N (g) in the incubated sample amount. n15.mass is in grams

data$tot.N<-with(data, (n15.mass*XN.e))

#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period)

data$Nfixd.sample<-with(data, tot.N*XNfixdX)

#Calculating N fixation to units of (gN fixed per gram of nodule per hour) "SNF efficiency"##

#gN fixed per sample divided by the mass of nodules that were incubated (n15.mass)

data$Nfixd<-with(data, (Nfixd.sample/n15.mass)*2)

#calculating total N fixed per plant #

#Taking "SNF efficiency" and multiplying by the total nodule mass 

data$c.mass[is.na(data$c.mass)] <- 0

data$n15.mass[is.na(data$n15.mass)] <- 0

data$nod.ex.mass[is.na(data$nod.ex.mass)] <- 0

data$totnodmass = with(data, (nod.ex.mass+n15.mass+c.mass))

data$tot.Nfixd <- with(data, (Nfixd*totnodmass))

#Replacing Nodule biomass ("tot.nod.mass") with N fixed ("tot.Nfixd") in glms

m3<-glm(tot.mass~S.num+nitro+light.per+X.col+tot.nod.mass+type,data=data, family= gaussian(link="log"))
summary(m3)

#by species
m3.p<-glm(tot.mass~light.per+nitro+X.col+tot.nod.mass+type,data=data.p, family= gaussian(link="log"))
summary(m3.p)

m3.z<-glm(tot.mass~light.per+nitro+X.col+tot.Nfixd+type,data=data.z, family= gaussian(link="log"))
summary(m3.z)

m3.s<-glm(tot.mass~light.per+nitro+X.col+tot.Nfixd,data=data.s, family= gaussian(link="log"))
summary(m3.s)
#'type' is removed from the model for S. microstachyum because there was no Nfixed data for S.microstachyum understory seedlings