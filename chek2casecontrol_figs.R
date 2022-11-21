setwd("/Users/alirezarasoulzadeh/Desktop/Roth Lab")

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
CARRIERS_populationbased_CHEK2_variant_carriers <- read_excel("CARRIERS populationbased CHEK2 variant carriers.xlsx", 
                                                              sheet = "CHEK2 carriers")
View(CARRIERS_populationbased_CHEK2_variant_carriers)

CHEK2_Imputed_Refined <- read.csv("~/Desktop/Roth Lab/CHEK2_Imputed_Refined.csv")

#extract rows with missense, and nonsense variants from 'CHEK2 carriers'
#there is an issue with synonymous mutations where there is amino acid sequence
#only nucleotide sequence is provided

chek2_carrier_filtered<-
  CARRIERS_populationbased_CHEK2_variant_carriers[
    CARRIERS_populationbased_CHEK2_variant_carriers$`variant type`=="stop_gained"
      |CARRIERS_populationbased_CHEK2_variant_carriers$`variant type`=="missense_variant",]


#add hgvs_pro column
chek2_carrier_filtered$hgvs_pro = str_extract(chek2_carrier_filtered$`CHEK2 variant`,"p.\\w*")

#match scores to case control dataset depending on the amino acid change observed
chek2_matched<-left_join(chek2_carrier_filtered,CHEK2_Imputed_Refined,by = "hgvs_pro")

write.csv(chek2_matched,file = "chek2_matched_carrierdataset.csv")

chek2_matched_clean<-na.omit(chek2_matched)
#222 different amino acid change are present in this data
length(levels(as.factor(chek2_matched_clean$hgvs_pro)))



chek2_matched_clean[chek2_matched_clean$`Age (case:AgeDx; control:Enroll)`==888,7]<-NA

#na omit the chek2_matched_clean again to remove rows that have unknown age
temp<-na.omit(chek2_matched_clean)
temp%>%group_by(hgvs_pro)%>%
  summarise(num_in_hgvspro = n(),average_age = mean(`Age (case:AgeDx; control:Enroll)`))


ggplot(temp, aes(temp$`Age (case:AgeDx; control:Enroll)`,temp$refined_score)) +
  geom_point() + theme_bw()+
  xlab("Patient Age") + ylab("Refined Score") + 
  labs(title = "Refined Score vs Patient Age - Missense Variants",color="Estrogen Receptor Status")



#filter out the 777 and 888 from the estrogen receptor status

temp<-chek2_matched_clean[chek2_matched_clean$`ER_status (1=positive,0=negative,777=NA,888=unknown)`==888|
     chek2_matched_clean$`ER_status (1=positive,0=negative,777=NA,888=unknown)`==777,]    

temp<-anti_join(chek2_matched_clean,temp,by="ER_status (1=positive,0=negative,777=NA,888=unknown)")
temp<-rename(temp, er_status = `ER_status (1=positive,0=negative,777=NA,888=unknown)`)
temp$er_status<-as.factor(temp$er_status)

ggplot(temp, aes(temp$`Age (case:AgeDx; control:Enroll)`,temp$refined_score,color=er_status)) +
  geom_point() + theme_bw()+
  xlab("Patient Age") + ylab("Refined Score") + 
  labs(title = "Refined Score vs Patient Age - Missense Variants",color="Estrogen Receptor Status")
  





