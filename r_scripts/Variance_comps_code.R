library(nlme)
library(sjstats)
library(lme4)
library(lme4)
library(ggplot2)

setwd("C:/Users/jdwr47/Documents/For_New_Comp/dense_UAV") # set to the root repository directory ("dense_UAV")


df <- read.csv("data/FinalC5ab.csv")

str(df)
df$Range <- as.factor(df$Range)
df$Row <- as.factor(df$Row)
df$Rep <- as.factor(df$Rep)
df$Pedigree <- as.factor(df$Pedigree)
df$Flight <- as.factor(df$Flight)

#### Raw data viz ######

p<-ggplot(df, aes(x=Flight, y=Elev.P99.Predictor1)) +
  geom_boxplot(size=0.5, outlier.shape = NA)+ 
  theme_bw()
p

p<-ggplot(df, aes(x=Flight, y=Elev.P99)) +
  geom_boxplot(size=0.5, outlier.shape = NA)+ 
  theme_bw()
p

hist(df$Asymptote)
hist(df$Asymptote1)

hist(df$Growth.Rate)
hist(df$Growth.Rate1)

hist(df$Inflection.Point)
hist(df$Inflection.Point1)

##############################################
####### ANOVA for  Asymptote #########
##############################################
df2 <- df[!duplicated(df$Plot_ID), ]
head(df2)

Asym <- lmer( Asymptote1  ~ (1|Pedigree)+
                #(1|Flight)+
                #(1|Management)+
                #(1|Pedigree:Flight)+
                #(1|Pedigree:Management)+
                #(1|Flight:Management)+
                #(1|Pedigree:Flight:Management)+
                (1|Range)+
                (1|Row)+
                (1|Rep), df2 )  


rmse(Asym)
summary(Asym)
cv(Asym)
MuMIn::r.squaredGLMM(Asym)
VC<-as.data.frame(print(VarCorr(Asym), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(Asym)[1,2]
VC$CV <- cv(Asym)
VC$Rmse <- rmse(Asym)
VC$Trait <- "Asymptote"
VT_Asym <- VC
VT_Asym


Blup_asym <- coef(Asym)$'Pedigree'
hist(Blup_asym$`(Intercept)`)

##############################################
###### ANOVA for  Growth.Rate ########
##############################################
Grw <- lmer( Growth.Rate1  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(Grw)
summary(Grw)
cv(Grw)
MuMIn::r.squaredGLMM(Grw)
VC<-as.data.frame(print(VarCorr(Grw), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(Grw)[1,2]
VC$CV <- cv(Grw)
VC$Rmse <- rmse(Grw)
VC$Trait <- "Growth rate"
VT_Grw <- VC
VT_Grw

Blup_grw <- coef(Grw)$'Pedigree'
hist(Blup_grw$`(Intercept)`)

##############################################
###### ANOVA for temporal Growth.Rate ########
##############################################
IP <- lmer( Inflection.Point1  ~ (1|Pedigree)+
              #(1|Flight)+
              #(1|Management)+
              #(1|Pedigree:Flight)+
              #(1|Pedigree:Management)+
              #(1|Flight:Management)+
              #(1|Pedigree:Flight:Management)+
              (1|Range)+
              (1|Row)+
              (1|Rep), df2 )  


rmse(IP)
summary(IP)
cv(IP)
MuMIn::r.squaredGLMM(IP)
VC<-as.data.frame(print(VarCorr(IP), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(IP)[1,2]
VC$CV <- cv(IP)
VC$Rmse <- rmse(IP)
VC$Trait <- "Inflection point"
VT_IP <- VC
VT_IP


Blup_IP <- coef(IP)$'Pedigree'
hist(Blup_IP$`(Intercept)`)


##############################################
########## ANOVA for  DTA ############
##############################################
DTA <- lmer( DTA  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(DTA)
summary(DTA)
cv(DTA)
MuMIn::r.squaredGLMM(DTA)
VC<-as.data.frame(print(VarCorr(DTA), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(DTA)[1,2]
VC$CV <- cv(DTA)
VC$Rmse <- rmse(DTA)
VC$Trait <- "DTA"
VT_DTA <- VC
VT_DTA


Blup_DTA <- coef(DTA)$'Pedigree'
hist(Blup_DTA$`(Intercept)`)



##############################################
########## ANOVA for  DTS ############
##############################################
DTS <- lmer( DTS  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(DTS)
summary(DTS)
cv(DTS)
MuMIn::r.squaredGLMM(DTS)
VC<-as.data.frame(print(VarCorr(DTS), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(DTS)[1,2]
VC$CV <- cv(DTS)
VC$Rmse <- rmse(DTS)
VC$Trait <- "DTS"
VT_DTS <- VC
VT_DTS


Blup_DTS <- coef(DTS)$'Pedigree'
hist(Blup_DTS$`(Intercept)`)


##############################################
########## ANOVA for  ASI ####################
##############################################

df2$ASI <- df2$DTS-df2$DTA

hist(df2$ASI)

ASI <- lmer( ASI  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(ASI)
summary(ASI)
cv(ASI)
MuMIn::r.squaredGLMM(ASI)
VC<-as.data.frame(print(VarCorr(ASI), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(ASI)[1,2]
VC$CV <- cv(ASI)
VC$Rmse <- rmse(ASI)
VC$Trait <- "ASI"
VT_ASI <- VC
VT_ASI


Blup_ASI <- coef(ASI)$'Pedigree'
hist(Blup_ASI$`(Intercept)`)
##############################################
########## ANOVA for  PHT ############
##############################################
hist(df2$Plant.Height..cm.)
PHT <- lmer( Plant.Height..cm.  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(PHT)
summary(PHT)
cv(PHT)
MuMIn::r.squaredGLMM(PHT)
VC<-as.data.frame(print(VarCorr(PHT), comp=c("Variance")))

VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(DTS)[1,2]
VC$CV <- cv(PHT)
VC$Rmse <- rmse(PHT)
VC$Trait <- "PHT"
VT_PHT <- VC
VT_PHT


Blup_PHT <- coef(PHT)$'Pedigree'
hist(Blup_PHT$`(Intercept)`)


##############################################
########## ANOVA for  YLD ############
##############################################
df2 <- read.csv("data/Masternote Washburn with yield.csv")

str(df2)
df2$Range <- as.factor(df2$Range)
df2$Row <- as.factor(df2$Row)
df2$Rep <- as.factor(df2$Rep)
df2$Pedigree <- as.factor(df2$Pedigree)

YLD <- lmer(Yield  ~ (1|Pedigree)+
               #(1|Flight)+
               #(1|Management)+
               #(1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df2 )  


rmse(YLD)
summary(YLD)
cv(YLD)
MuMIn::r.squaredGLMM(YLD)
VC<-as.data.frame(print(VarCorr(YLD), comp=c("Variance")))

VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(YLD)[1,2]
VC$CV <- cv(YLD)
VC$Rmse <- rmse(YLD)
VC$Trait <- "YLD"
VT_YLD <- VC
VT_YLD

str(df2)
Blup_YLD <- coef(YLD)$'Pedigree'
hist(Blup_YLD$`(Intercept)`)

##############################################
########## ANOVA for  EHT ############
##############################################
df2 <- read.csv("data/Masternote Washburn with yield.csv")

str(df2)
df2$Range <- as.factor(df2$Range)
df2$Row <- as.factor(df2$Row)
df2$Rep <- as.factor(df2$Rep)
df2$Pedigree <- as.factor(df2$Pedigree)

#df2$EHTdPHT <- df2$Ear.Height..cm./df2$Plant.Height..cm.

EHT <- lmer(Ear.Height..cm. ~ (1|Pedigree)+
              #(1|Flight)+
              #(1|Management)+
              #(1|Pedigree:Flight)+
              #(1|Pedigree:Management)+
              #(1|Flight:Management)+
              #(1|Pedigree:Flight:Management)+
              (1|Range)+
              (1|Row)+
              (1|Rep), df2 )  


rmse(EHT)
summary(EHT)
cv(EHT)

MuMIn::r.squaredGLMM(EHT)
VC<-as.data.frame(print(VarCorr(EHT), comp=c("Variance")))

VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[1,5] / (VC[1,5] + (VC[5,5]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(EHT)[1,2]
VC$CV <- cv(EHT)
VC$Rmse <- rmse(EHT)
VC$Trait <- "EHT"
VT_EHT <- VC
VT_EHT

str(df2)
Blup_EHT <- coef(EHT)$'Pedigree'
hist(Blup_EHT$`(Intercept)`)







##############################################
######### ANOVA for temporal PHTs ############
##############################################

Tpht<- lmer( Elev.P99.Predictor1  ~ (1|Pedigree)+
               (1|Flight)+
               #(1|Management)+
               (1|Pedigree:Flight)+
               #(1|Pedigree:Management)+
               #(1|Flight:Management)+
               #(1|Pedigree:Flight:Management)+
               (1|Range)+
               (1|Row)+
               (1|Rep), df )  


rmse(Tpht)
summary(Tpht)
cv(Tpht)
MuMIn::r.squaredGLMM(Tpht)
VC<-as.data.frame(print(VarCorr(Tpht ), comp=c("Variance")))
VC$Percent<-VC$vcov/sum(VC$vcov)
VC$Percent<-round(VC$Percent*100,1)
heritability <- VC[2,6] / (VC[2,6] + (VC[7,6]/2))
VC$Heritability <- heritability
VC$Rsquared <- MuMIn::r.squaredGLMM(Tpht)[1,2]
VC$CV <- cv(Tpht)
VC$Rmse <- rmse(Tpht)
VC$Trait <- "Temporal plant height"
VT_Tpht <- VC
VT_Tpht


Flight <- ranef(Tpht)$'Flight'
Pedigree <- ranef(Tpht)$'Pedigree'
fp <- coef(Tpht)$'Pedigree:Flight'
fp[,2] <- rownames(fp)


fp$Pedigree <- lapply(strsplit(as.character(fp$V2), "\\:"), "[", 1)
fp$DAP <- lapply(strsplit(as.character(fp$V2), "\\:"), "[", 2)
head(fp)
fp <- as.data.frame(lapply(fp, unlist))
head(fp)
names(fp)[1:2] <- c("Height","interaction" )

for (i in 1:nrow(fp)){
  fp$Height[i]<-fp$Height[i]+Flight$`(Intercept)`[row.names(Flight)==fp$DAP[i]]
}

for (i in 1:nrow(fp)){
  fp$Height[i]<-fp$Height[i]+Pedigree$`(Intercept)`[row.names(Pedigree)==fp$Pedigree[i]]
}

head(fp)


####################################
############ PCA of TPHTs ##########
####################################
library(tidyverse)
library(factoextra)
dt <- fp[,c(1,3,4)]
head(dt)
data <- reshape(dt, idvar = "Pedigree", timevar = "DAP", direction = "wide")
row.names(data) <- data$Pedigree
data <- data[,-1]
data <- data %>% mutate_if(is.numeric, function(x) ifelse(is.na(x), median(x, na.rm = T), x))
pc <- prcomp(data, , center = TRUE, scale. = TRUE)
df  <-pc$x
write.csv(df, "data/pca.csv")


df <- read.csv("data/pca.csv")
head(df)
df1 <- df[(!df$Group %in% c("Canadian", "Other")),]



p_pca <- ggplot(df1, aes(x=PC1, y=PC2, group=Group)) +
  geom_point(aes(shape=Group, color=Group))+theme_classic()

jpeg("data/pca_Tphts.jpeg", width = 6,height = 4,units = "in", res=300)
p_pca

dev.off()



#######################
#### Viz. of TPHTs ####
#######################



library("RColorBrewer")
fp$DAP <-  as.numeric(fp$DAP)
library(dplyr)
head(df)
head(fp)

fp$Height2 <- ifelse(fp$Height <0 , 0, fp$Height)

fp$Tester <- lapply(strsplit(as.character(fp$Pedigree), "\\X "), "[", 2)
fp <- as.data.frame(lapply(fp, unlist))

#fp_subset <- subset(fp, Tester=="PHZ51" | Tester=="PHK76" | Tester=="PHP02")

fp_subset <- fp[grepl('W10004', fp$Pedigree), ]
write.csv(fp_subset, "data/Tpht_subset.csv")


library(stringr)
library(dplyr)

#fp_subset <- fp %>% filter(str_detect(fp$Pedigree , "W10004"))

p1 <- ggplot( fp_subset, aes(x=DAP, y=Height, group=interaction(DAP, Tester),  fill=Tester)) + 
  # geom_smooth(aes(x=DAP, y=CHM), alpha=0.1,se = FALSE)+ 
  #scale_color_brewer(palette = "Dark2")+
  geom_boxplot(position=position_dodge(), width=1.2, size=0.4, outlier.size = 0.5)+ 
  #  geom_vline(xintercept = 59, linetype="dashed",  color = "#4DAF4A", size=1.5)+
  #  geom_vline(xintercept = 68, linetype="dashed", color = "#984EA3", size=1.5)+
  #  scale_color_manual(values = c("#4DAF4A", "#984EA3" ))+
  theme_bw()+ 
  #  facet_grid(~Management) +
  scale_y_continuous(name= "Temporal plant height")+
  theme(legend.position="right") + 
  #  labs(title = "A-) Pedigree:Flight:Management interaction in 2017")+
  scale_x_continuous(name= "Days after planting (DAP)",breaks=  as.numeric(names(table(fp$DAP))) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position =c(0.85, 0.25),
        #  legend.key.height= unit(0.8, 'cm'),
        #  legend.key.width= unit(0.5, 'cm'),
        #  legend.background = element_rect(fill="White",
        #                                   size=0.5,
        #                                   colour ="black",
        #                                   linetype="solid")
  )
#  guides(color = guide_legend(override.aes = list(size = 1,alpha = 1) ) )

p1 

getwd()

jpeg("Figures/Tphts.jpeg", width = 15,height = 6,units = "in", res=300)
p1
#grid.arrange(p, p1, p2, p3, ncol=2)

dev.off()


VC <- rbind(
  VT_DTA,
  VT_DTS,
  VT_Asym,
  VT_Grw,
  VT_IP,
  VT_PHT,
  VT_YLD,
  VT_EHT,
  VT_ASI)

VC$Type <- "Single"
VT_Tpht$Type <- "Temporal"

df <- rbind(VC,VT_Tpht)

ylim.prim <- c(0, 100)   # in this example, precipitation
ylim.sec <- c(0,1)    # in this example, temperature

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

df$grp <- factor(df$grp, levels = c( "Pedigree", "Flight",
                                     "Pedigree:Flight",
                                     "Range", 
                                     "Row", 
                                     "Rep",
                                     "Residual"))

df$Trait  <- factor(df$Trait, levels = c( "Asymptote", 
                                          "Growth rate",
                                          "Inflection point",
                                          "YLD",
                                          "DTA", 
                                          "DTS",
                                          "ASI",
                                          "EHT",
                                          "PHT",
                                          "Temporal plant height"))


head(df)
write.csv(df, "data/Var_comps_manual_phenos_Tpht.csv")

p2 <-  ggplot(data=df, aes(x= Trait, y=Percent)) +
  geom_col(mapping=aes(fill=grp), width= 0.9) +
  geom_point(mapping=aes(y = a +Heritability*b,shape="Repeatability"), fill="white",  color="black",size=4, alpha=1) +
  geom_point(mapping=aes(y = a +Rsquared*b,shape="Rsquared"), fill="white",  color="black",size=4, alpha=1) +
  geom_point(mapping=aes(y = a +CV*b,shape="CV"), fill="white",  color="black",size=4, alpha=1) +
  #  geom_point(mapping=aes(y = a + `Rsquared`*b,shape="Rsquared"),  alpha=1) +
  #  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous("Explained percent variation (%)", sec.axis = sec_axis(~ (. - a)/b, name = "Repeatability, Rsquared and CV")) + 
  #  scale_fill_manual(values = c("#B15928","#A6CEE3","#6A3D9A","#CAB2D6","#FF7F00","#FDBF6F","#E31A1C","#FB9A99","#33A02C" ,"#B2DF8A","#1F78B4" ))+
  #  scale_fill_brewer(palette = "Paired")+ 
  #  scale_colour_manual(values =c("black","black"))+
  scale_shape_manual(values=c(22,25,24))+
  #scale_fill_viridis_d()+
  theme_bw() + 
  scale_x_discrete("")+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 12.5)
        #      legend.key.width = unit(0.2, 'cm'),
        #      legend.key.height = unit(0.1, 'cm')
  )+
  facet_grid(~Type, scales = "free", space = "free")+
  #  labs(title = "A-) Multispectral HTP platform")+
  guides(fill=guide_legend(title="Variance \ncomponent")) 

p2



Blup_asym$Hybrid <-  rownames(Blup_asym)
Blup_asym$Trait <- "Asymptote"

Blup_PHT$Hybrid <-  rownames(Blup_PHT)
Blup_PHT$Trait <- "Plant height"

Blup_EHT$Hybrid <-  rownames(Blup_EHT)
Blup_EHT$Trait <- "Ear height"

Blup_DTA$Hybrid <-  rownames(Blup_DTA)
Blup_DTA$Trait <- "DTA"

Blup_DTS$Hybrid <-  rownames(Blup_DTS)
Blup_DTS$Trait <- "DTS"

Blup_ASI$Hybrid <-  rownames(Blup_ASI)
Blup_ASI$Trait <- "ASI"

Blup_grw$Hybrid <-  rownames(Blup_grw)
Blup_grw$Trait <- "Growth rate"

Blup_IP$Hybrid <-  rownames(Blup_IP)
Blup_IP$Trait <- "Inflection point"

Blup_YLD$Hybrid <-  rownames(Blup_YLD)
Blup_YLD$Trait <- "Yield"

dff <- rbind(Blup_YLD, Blup_ASI, Blup_EHT, Blup_asym,Blup_PHT, Blup_DTA, Blup_DTS, Blup_grw, Blup_IP)
head(dff)

dff$Tester <- lapply(strsplit(as.character(dff$Hybrid), "\\X "), "[", 2)
dff <- as.data.frame(lapply(dff, unlist))
head(dff)

dff <- dff[grepl('W10004', dff$Hybrid ), ]
names(dff)[1] <- "Data"
head(dff)

Blup <- ggplot(dff, aes(x=Data, color=Tester)) +
  #  geom_histogram(fill="white", position="dodge")+
  geom_density(alpha=0.6)+
  theme(legend.position="top") +
  geom_rug()+
  facet_wrap(~Trait, scales = "free")+
  theme_bw()+
  theme(legend.position = "top",
        legend.background =  element_rect(fill = alpha("gray80", 0.3) , size=0.1, linetype="solid"),
        legend.key.height = unit(0.1, 'cm'),
        legend.key.width =  unit(0.4, 'cm'))+
  scale_x_continuous("") 
Blup

write.csv(dff, "data/blups_manual_phenos.csv")

library(reshape2)

cordata <- dcast(dff, Tester+ Hybrid~Trait, value.var = "Data")
head(cordata)
library(GGally)

ggpairs(cordata, columns = 3:9, ggplot2::aes(colour=Tester)) 
library(gridExtra)
grid.arrange(p2, Blup, p1, layout_matrix = rbind(c(1,2),
                                                 c(3)))




jpeg("Figures/Fig3_r.jpeg", width =12,height = 8,units = "in", res=500)

grid.arrange(p2, Blup, p1, layout_matrix = rbind(c(1,1,1,2,2),
                                                 c(3)))

dev.off()
