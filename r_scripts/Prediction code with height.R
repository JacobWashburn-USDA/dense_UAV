
#setwd("C:/Users/jdwr47/Documents/For_New_Comp/dense_UAV") # set to the root repository directory ("dense_UAV")
setwd("~/Desktop/washburn_repos/dense_UAV/data/")

df <- read.csv("Phenomic_scaled_centered.csv", row.names = 1)

setwd("~/Desktop/washburn_repos/dense_UAV/data/Prediction ability/temp_files/")
df[1:5,1:5]
df_sub <- df
flight <- as.data.frame(colnames(df_sub))
flight$DAP <- lapply(strsplit(as.character(flight$`colnames(df_sub)`), "\\."), "[", 2)
flight <- as.data.frame(lapply(flight, unlist))
length(unique(flight$DAP))
unique(flight$DAP)

ncol(df_sub)
P43 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P43,"P43.csv")

df_sub <- dplyr::select(df, -contains( c("126") ))
ncol(df_sub)
P42 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P42,"P42.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122") ))
ncol(df_sub)
P41 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P41,"P41.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118") ))
P40 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P40,"P40.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114") ))
P39 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P39,"P39.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111") ))
P38 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P38,"P38.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107") ))
P37 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P37,"P37.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104") ))
P36 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P36,"P36.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101") ))
P35 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P35,"P35.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99") ))
P34 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P34,"P34.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96"  ) ))
P33 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P33,"P33.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93"  ) ))
P32 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P32,"P32.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90"  ) ))
P31 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P31,"P31.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ) ))
P30 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P30,"P30.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85" ) ))
P29 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P29,"P29.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83") ))
P28 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P28,"P28.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80"  ) ))
P27 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P27,"P27.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78"  ) ))
P26 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P26,"P26.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76"  ) ))
P25 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P25,"P25.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73"  ) ))
P24 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P24,"P24.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71"  ) ))
P23 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P23,"P23.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69" ) ))
P22 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P22,"P22.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69","66"  ) ))
P21 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P21,"P21.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64" ) ))
P20 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P20,"P20.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62"  ) ))
P19 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P19,"P19.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60"  ) ))
P18 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P18,"P18.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90","87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56"   ) ))
P17 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P17,"P17.csv")

df_sub <- dplyr::select(df, -contains(c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                        "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53"   ) ))
P16 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P16,"P16.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51"   )  ))
P15 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P15,"P15.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48"   )  ))
P14 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P14,"P14.csv")

df_sub <- dplyr::select(df,-contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                        "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45"   )  ))
P13 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P13,"P13.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41"   )  ))
P12 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P12,"P12.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39"  )  ))
P11 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P11,"P11.csv")

df_sub <- dplyr::select(df,  -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                          "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37"  )  ))
P10 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P10,"P10.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34")  ))
P9 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P9,"P9.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32")  ))
P8 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P8,"P8.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30")  ))
P7 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P7,"P7.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30","27")  ))
P6 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P6,"P6.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30","27","24")  ))
P5 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P5,"P5.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30","27","24","22")  ))
P4 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P4,"P4.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30","27","24","22","17")  ))
P3 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P3,"P3.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48","45","41" ,"39","37",
                                         "34","32","30","27","24","22","17","15")  ))
P2 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P2,"P2.csv")

df_sub <- dplyr::select(df, -contains( c("126", "122", "118" ,"114", "111", "107", "104", "101","99","96","93","90",
                                         "87" ,"85","83","80","78","76" ,"73" ,"71","69","66","64","62","60","56","53","51","48",
                                         "45","41" ,"39","37","34","32","30","27","24","22","17","15","10")  ))
P1 <- tcrossprod(as.matrix(df_sub))/ncol(as.matrix(df_sub))
write.csv(P1,"P1.csv")


library(dplyr)
library(BGLR)

setwd("~/Desktop/washburn_repos/dense_UAV/data/")
AM <-as.matrix(read.csv("AM.csv",row.names = 1))
AM[1:10,1:10]

DM <-as.matrix(read.csv("DM.csv",row.names = 1))
DM[1:10,1:10]

setwd("~/Desktop/washburn_repos/dense_UAV/data/Prediction ability/temp_files/")
P1 <-as.matrix(read.csv("P1.csv",row.names = 1))
P2 <-as.matrix(read.csv("P2.csv",row.names = 1))
P3 <-as.matrix(read.csv("P3.csv",row.names = 1))
P4 <-as.matrix(read.csv("P4.csv",row.names = 1))
P5 <-as.matrix(read.csv("P5.csv",row.names = 1))
P6 <-as.matrix(read.csv("P6.csv",row.names = 1))
P7 <-as.matrix(read.csv("P7.csv",row.names = 1))
P8 <-as.matrix(read.csv("P8.csv",row.names = 1))
P9 <-as.matrix(read.csv("P9.csv",row.names = 1))
P10 <-as.matrix(read.csv("P10.csv",row.names = 1))
P11 <-as.matrix(read.csv("P11.csv",row.names = 1))
P12 <-as.matrix(read.csv("P12.csv",row.names = 1))
P13 <-as.matrix(read.csv("P13.csv",row.names = 1))
P14 <-as.matrix(read.csv("P14.csv",row.names = 1))
P15 <-as.matrix(read.csv("P15.csv",row.names = 1))
P16 <-as.matrix(read.csv("P16.csv",row.names = 1))
P17 <-as.matrix(read.csv("P17.csv",row.names = 1))
P18 <-as.matrix(read.csv("P18.csv",row.names = 1))
P19 <-as.matrix(read.csv("P19.csv",row.names = 1))
P20 <-as.matrix(read.csv("P20.csv",row.names = 1))
P21 <-as.matrix(read.csv("P21.csv",row.names = 1))
P22 <-as.matrix(read.csv("P22.csv",row.names = 1))
P23 <-as.matrix(read.csv("P23.csv",row.names = 1))
P24 <-as.matrix(read.csv("P24.csv",row.names = 1))
P25 <-as.matrix(read.csv("P25.csv",row.names = 1))
P26 <-as.matrix(read.csv("P26.csv",row.names = 1))
P27 <-as.matrix(read.csv("P27.csv",row.names = 1))
P28 <-as.matrix(read.csv("P28.csv",row.names = 1))
P29 <-as.matrix(read.csv("P29.csv",row.names = 1))
P30 <-as.matrix(read.csv("P31.csv",row.names = 1))
P31 <-as.matrix(read.csv("P31.csv",row.names = 1))
P32 <-as.matrix(read.csv("P32.csv",row.names = 1))
P33 <-as.matrix(read.csv("P33.csv",row.names = 1))
P34 <-as.matrix(read.csv("P34.csv",row.names = 1))
P35 <-as.matrix(read.csv("P35.csv",row.names = 1))
P36 <-as.matrix(read.csv("P36.csv",row.names = 1))
P37 <-as.matrix(read.csv("P37.csv",row.names = 1))
P38 <-as.matrix(read.csv("P38.csv",row.names = 1))
P39 <-as.matrix(read.csv("P39.csv",row.names = 1))
P40 <-as.matrix(read.csv("P40.csv",row.names = 1))
P41 <-as.matrix(read.csv("P41.csv",row.names = 1))
P42 <-as.matrix(read.csv("P42.csv",row.names = 1))
P43 <-as.matrix(read.csv("P43.csv",row.names = 1))

#####

setwd("~/Desktop/washburn_repos/dense_UAV/data/")
Pheno <- read.csv("Pheno_blues.csv")
head(Pheno)


# Kernels #
Eta1<-list(A=list(K=AM,model="RKHS"), D=list(K=DM,model="RKHS"))

Eta2<-list(P=list(K=P1,model="RKHS"))
Eta3<-list(P=list(K=P2,model="RKHS"))
Eta4<-list(P=list(K=P3,model="RKHS"))
Eta5<-list(P=list(K=P4,model="RKHS"))
Eta6<-list(P=list(K=P5,model="RKHS"))
Eta7<-list(P=list(K=P6,model="RKHS"))
Eta8<-list(P=list(K=P7,model="RKHS"))
Eta9<-list(P=list(K=P8,model="RKHS"))
Eta10<-list(P=list(K=P9,model="RKHS"))
Eta11<-list(P=list(K=P10,model="RKHS"))
Eta12<-list(P=list(K=P11,model="RKHS"))
Eta13<-list(P=list(K=P12,model="RKHS"))
Eta14<-list(P=list(K=P13,model="RKHS"))
Eta15<-list(P=list(K=P14,model="RKHS"))
Eta16<-list(P=list(K=P15,model="RKHS"))
Eta17<-list(P=list(K=P16,model="RKHS"))
Eta18<-list(P=list(K=P17,model="RKHS"))
Eta19<-list(P=list(K=P18,model="RKHS"))
Eta20<-list(P=list(K=P19,model="RKHS"))
Eta21<-list(P=list(K=P20,model="RKHS"))
Eta22<-list(P=list(K=P21,model="RKHS"))
Eta23<-list(P=list(K=P22,model="RKHS"))
Eta24<-list(P=list(K=P23,model="RKHS"))
Eta25<-list(P=list(K=P24,model="RKHS"))
Eta26<-list(P=list(K=P25,model="RKHS"))
Eta27<-list(P=list(K=P26,model="RKHS"))
Eta28<-list(P=list(K=P27,model="RKHS"))
Eta29<-list(P=list(K=P28,model="RKHS"))
Eta30<-list(P=list(K=P29,model="RKHS"))
Eta31<-list(P=list(K=P30,model="RKHS"))
Eta32<-list(P=list(K=P31,model="RKHS"))
Eta33<-list(P=list(K=P32,model="RKHS"))
Eta34<-list(P=list(K=P33,model="RKHS"))
Eta35<-list(P=list(K=P34,model="RKHS"))
Eta36<-list(P=list(K=P35,model="RKHS"))
Eta37<-list(P=list(K=P36,model="RKHS"))
Eta38<-list(P=list(K=P37,model="RKHS"))
Eta39<-list(P=list(K=P38,model="RKHS"))
Eta40<-list(P=list(K=P39,model="RKHS"))
Eta41<-list(P=list(K=P40,model="RKHS"))
Eta42<-list(P=list(K=P41,model="RKHS"))
Eta43<-list(P=list(K=P42,model="RKHS"))
Eta44<-list(P=list(K=P43,model="RKHS"))


#### Cross-Validation ####
hybrid<-as.character(unique(Pheno$Hybrid))
head(Pheno)
#### Within-Environment Prediction ####
set.seed(789)
CV1 <- list()
#CV2 <- list()
library(foreach)
library(BGLR)

Models <- list(Eta1,Eta2,Eta3,Eta4,Eta5,Eta6,Eta7,Eta8,Eta9,
               Eta10,Eta11,Eta12,Eta13,Eta14,Eta15,Eta16,Eta17,
               Eta18,Eta19,Eta20,Eta21,Eta22,Eta23,Eta24,Eta25,
               Eta26,Eta27,Eta28,Eta29,Eta30,Eta31,Eta32,Eta33,
               Eta34,Eta35,Eta36,Eta37,Eta38,Eta39,Eta40,Eta41,
               Eta42,Eta43,Eta44
)

library(parallel)
library(doParallel)
library(MASS)

numCores <- detectCores()-1
numCores
registerDoParallel(numCores)

setwd("~/Desktop/washburn_repos/dense_UAV/data/Prediction ability/")

foreach(MODEL = 1:length(Models), .packages = c("BGLR", "dplyr")) %dopar% {
  for (rep_num in 1:25)
  {
    
    x <- (rep_num - 1) * 5 # Allows saving in lists as rep_num increases
    
    k <- 5 # Set the number of folds
    obs_per_fold <- ceiling(nrow(Pheno)/k) # Calculate the number of observations per fold
    folds <-sample(rep(1:k, each = obs_per_fold, length.out = nrow(Pheno)))
    # indices into k approximately equal-sized subsets
    df <- as.data.frame(cbind(hybrid, folds))
    
    for(fold_num in 1:max(folds))
    {
      
      test_geno <- subset(df, grepl(fold_num, folds))$hybrid
      train_geno <- subset(df, !grepl(fold_num, folds))$hybrid
      
      yield <- Pheno[,c(1,2)] # Yield
      yield$Y2<-NA
      yield$Y2[yield$Hybrid%in%train_geno] <- yield[,2][yield$Hybrid%in%train_geno]
      y_t<-as.numeric(yield$Y2)
      
      
      fit<-BGLR(y=y_t,ETA=Models[[MODEL]],nIter=5000,burnIn=1000, thin=10)
      yield$yhat <- fit$yHat
      # CV1_M1 #
      df_test <- subset(yield, yield$Hybrid %in% test_geno)
      CV1[[(fold_num+x)]] <- as.data.frame(df_test %>% dplyr::summarize(cor=cor(df_test[,2], yhat,use = "complete.obs")))
      # CV2_M1 #
      #df_train <- subset(yield, yield$Hybrid %in% train_geno)
      #CV2[[(fold_num+x)]] <- as.data.frame(df_train %>% dplyr::summarize(cor=cor(df_train[,2], yhat,use = "complete.obs")))
      
    }
    
    if (rep_num == 25) {
      CV1out <- plyr::ldply(CV1, data.frame)
      #CV2out <- plyr::ldply(CV2, data.frame)
      write.csv(CV1out,file = paste("CV1_", "Eta", MODEL, ".csv", sep = ""), row.names = F)
      #write.csv(CV2out,file = paste("CV2_", "Eta", MODEL, ".csv", sep = ""), row.names = F)
      
    }
  }
}

#setwd("D:\\Washburn\\Prediction new\\Pred ability")
#setwd("~/Desktop/washburn_repos/dense_UAV/data/Prediction ability/")
setwd("C:/Users/jdwr47/Documents/For_New_Comp/dense_UAV/data/Prediction ability")

library(readr)
list_csv_files <- list.files(path = ".")
list_csv_files
df <-as.data.frame(readr::read_csv(list_csv_files, id = "file_name"))
head(df)

df$file_name <- gsub(".csv", "",df$file_name)
df$file_name <- gsub("Eta", "M",df$file_name)


df$CV <- lapply(strsplit(as.character(df$file_name ), "\\_"), "[", 1)
df$Model <- lapply(strsplit(as.character(df$file_name), "\\_"), "[", 2)
df        <- as.data.frame(lapply(df, unlist))

head(df)


write.csv(df,"Prediction ability with hieght rerun.csv")
library(ggplot2)


df <- read.csv("Prediction ability.csv")
head(df)
df <- as.data.frame(df %>%  dplyr::group_by(Model,CV) %>% 
                      dplyr::summarise(M = mean(cor, na.rm=TRUE),
                                       SD = sd(cor, na.rm=TRUE)))
df$Type <- ifelse(df$Model=="M1","Genomic", "Phenomic")
df$Model <- factor(df$Model, levels =  c("M1",
                                         "M2",
                                         "M3",
                                         "M4",
                                         "M5",
                                         "M6",
                                         "M7",
                                         "M8",
                                         "M9",
                                         "M10",
                                         "M11",
                                         "M12",
                                         "M13",
                                         "M14",
                                         "M15",
                                         "M16",
                                         "M17",
                                         "M18",
                                         "M19",
                                         "M20",
                                         "M21",
                                         "M22",
                                         "M23",
                                         "M24",
                                         "M25",
                                         "M26",
                                         "M27",
                                         "M28",
                                         "M29",
                                         "M30",
                                         "M31",
                                         "M32",
                                         "M33",
                                         "M34",
                                         "M35",
                                         "M36",
                                         "M37",
                                         "M38",
                                         "M39",
                                         "M40",
                                         "M41",
                                         "M42",
                                         "M43",
                                         "M44"))

p <- ggplot(data=subset(df, df$CV=="CV1")  , aes(x=Model, y=M, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=format(round(M,3),nsmall = 3)  ), hjust=3, color="white",
            position = position_dodge(0.9), angle = 90,size=3.5)+
  geom_errorbar(aes(ymin=M-SD, ymax=M+SD), width=.2,
                position=position_dodge(.9))+
  #  scale_fill_brewer(palette="Paired")+
  theme_bw()+
  scale_y_continuous("Prediction ability")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

jpeg("Pred_Yield.jpeg", 
     width = 10,height = 4,units = "in", res=300)
p
dev.off()


