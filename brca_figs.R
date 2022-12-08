setwd("/Users/alirezarasoulzadeh/Desktop/Roth Lab")

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)


brca_cosmic<-
  read.csv("/Users/alirezarasoulzadeh/Desktop/Roth Lab/brca1data_cosmic.csv")


brca_cosmic2<-brca_cosmic[,c(1,8,12,22,33,35)] 
brca_cosmic2 <- rename(brca_cosmic2,hgvs_pro=HGVSP)
brca_comsic_joined <- left_join(brca_cosmic2,CHEK2_Imputed_Refined,by="hgvs_pro")

#only keep those with known mutations
chek2_cosmic_cleanjoin<-chek2_comsic_joined[chek2_comsic_joined$Mutation_description!="Unknown",]
#keep those with scores
chek2_cosmic_cleanjoin<-chek2_comsic_joined[!is.na(chek2_comsic_joined$refined_score),]
