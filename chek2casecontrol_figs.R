setwd("/Users/alirezarasoulzadeh/Desktop/Roth Lab")

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
CARRIERS_populationbased_CHEK2_variant_carriers <- read_excel("CARRIERS populationbased CHEK2 variant carriers.xlsx", 
                                                              sheet = "CHEK2 carriers")
CARRIERS_populationbased_CHEK2_variant_carriers<-read.csv("/Users/alirezarasoulzadeh/Desktop/Roth Lab/chek2_variant_cleaned.csv")

View(CARRIERS_populationbased_CHEK2_variant_carriers)

CHEK2_Imputed_Refined <- read.csv("~/Desktop/Roth Lab/CHEK2_Imputed_Refined.csv")

#extract rows with missense, and nonsense variants from 'CHEK2 carriers'
#there is an issue with synonymous mutations where there is amino acid sequence
#only nucleotide sequence is provided
#write.csv(CARRIERS_populationbased_CHEK2_variant_carriers,file="chek2variant_carriers.csv")
#write.csv(aa_seq_chek2,file="aa_seq_chek2.csv",row.names = FALSE)

chek2_carrier_filtered<-
  CARRIERS_populationbased_CHEK2_variant_carriers[
    CARRIERS_populationbased_CHEK2_variant_carriers$`variant.type`=="stop_gained"
      |CARRIERS_populationbased_CHEK2_variant_carriers$`variant.type`=="missense_variant"
    |CARRIERS_populationbased_CHEK2_variant_carriers$`variant.type`=="synonymous_variant",]


#add hgvs_pro column
chek2_carrier_filtered$hgvs_pro = str_extract(chek2_carrier_filtered$`CHEK2.variant`,"p.\\w*")

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

ggplot(temp, aes(x=temp$`Age (case:AgeDx; control:Enroll)`,y=temp$refined_score))+
  geom_histogram(stat="identity")

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
  
chek2_cosmic <- read.csv("/Users/alirezarasoulzadeh/Desktop/Roth Lab/chek2data_cosmic.csv")

chek2_cosmic2<-chek2_cosmic[,c(1,8,12,22,33,35)] 
levels(as.factor(chek2_cosmic$Primary_site)) #36 different tissue types in data
chek2_cosmic2 <- rename(chek2_cosmic2,hgvs_pro=HGVSP)

chek2_comsic_joined <- left_join(chek2_cosmic2,CHEK2_Imputed_Refined,by="hgvs_pro")
#chek2 data with tissue distribution from cosmic joined to the imputed scores for chek2 
#variant effect map

#only keep those with known mutations
chek2_cosmic_cleanjoin<-chek2_comsic_joined[chek2_comsic_joined$Mutation_description!="Unknown",]
#keep those with scores
chek2_cosmic_cleanjoin<-chek2_comsic_joined[!is.na(chek2_comsic_joined$refined_score),]

temp<-chek2_cosmic_cleanjoin%>%group_by(Primary_site)%>%summarise(count_tissue_occurance = n())
#use these as the organs to look at 
filter(temp,count_tissue_occurance>50)[,1]
myvector<-c("breast","central_nervous_system","endometrium",
            "haematopoietic_and_lymphoid_tissue","kidney","large_intestine",
            "liver","lung","meninges","oesophagus","ovary","prostate",
            "skin","soft_tissue","stomach","thyroid","upper_aerodigestive_tract",
            "urinary_tract"
            )
temp2<-chek2_cosmic_cleanjoin %>% filter(Primary_site %in% myvector)

ggplot(temp2, aes(y=refined_score))+
  geom_boxplot() + facet_wrap(~Primary_site) + xlab(NULL)
##################################################################


missense_chek2 <- read.csv("~/Desktop/Roth Lab/missense_chek2.csv")
synonymous_chek2 <- read.csv("~/Desktop/Roth Lab/synonymous_chek2.csv")
nonsense_chek2 <- read.csv("~/Desktop/Roth Lab/nonsense_chek2.csv")

synonymous_chek2$aa_position<-NULL

cleaned_chek2_cases <-rbind(missense_chek2,synonymous_chek2)
cleaned_chek2_cases<-rbind(cleaned_chek2_cases,nonsense_chek2)

cleaned_chek2_matched<-left_join(cleaned_chek2_cases,CHEK2_Imputed_Refined,by = "hgvs_pro")

cleaned_chek2_matched<-cleaned_chek2_matched[!is.na(cleaned_chek2_matched$refined_score),]

ggplot(temp2, aes(x=refined_score))+
  geom_histogram() + facet_wrap(~Primary_site) + xlab(NULL) +
  theme_bw()

chek2_cosmic_cleanjoin$lower_score<-chek2_cosmic_cleanjoin$refined_score-
  chek2_cosmic_cleanjoin$refined_score_se

chek2_cosmic_cleanjoin$upper_score<-chek2_cosmic_cleanjoin$refined_score+
  chek2_cosmic_cleanjoin$refined_score_se

range_checker<-function(row){
  if(as.numeric(row[17])<0.5 && as.numeric(row[18])>0.5){
    return(1) #if the upper range and lower range cross 0.5 assign a value for subsetting
  }else{0} #score we are more sure about 
}

#for some reason the apply function did not work for what I was doing 
range_mask<-vector("logical",4283)
for(i in (1:4283)){
  row<-chek2_cosmic_cleanjoin[i,]
  range_mask[i]<-range_checker(row)
  
}

chek2_cosmic_cleanjoin<-cbind(chek2_cosmic_cleanjoin,range_mask)
chek2_cosmic_tissue<-chek2_cosmic_cleanjoin %>% filter(range_mask == FALSE)


ggplot(chek2_cosmic_cleanjoin, aes(x=refined_score))+
  geom_histogram() + facet_wrap(~Primary_site) + xlab(NULL) +
  theme_bw() +geom_text(data=tissue_summary,mapping = aes(x=-Inf,y=Inf, label=labels,
                                      hjust = -0.5,vjust   = 1.5)) +geom_vline(
                      xintercept = 0.5, linetype = "dotted"
                                      ) + labs(title="Chek2 Variant Tissue Distribution",
                                               caption = "scores that crossed 0.5 were included")

tissue_summary<-chek2_cosmic_cleanjoin %>% group_by(Primary_site) %>% summarise(tissue_occurence = n())
tissue_summary<-as.data.frame(tissue_summary)
tissue_summary$labels <-apply(tissue_summary,1,label_maker)

label_maker<-function(row){
  return(sprintf("n = %s",row[2]))
}

ggplot(chek2_cosmic_tissue, aes(x=type ,y=refined_score))+
  geom_boxplot(outlier.colour="white") + facet_wrap(~Primary_site) + xlab(NULL) +
  theme_bw() +geom_text(data=tissue_summary,mapping = aes(x=-Inf,y=Inf, label=labels,
                                         hjust = -2,vjust   = 2.5))+
  geom_hline(yintercept = 0.5, linetype= "dotted")+
                            labs(title="Chek2 Variant Tissue Distribution",
                         caption = "scores that crossed 0.5 were omitted,
                         lower scores - worse",color="Mutation Type") +
 geom_jitter(width = 0.1,aes(color=chek2_cosmic_tissue$type))





#################################################################
#filter out scores we are not too certain of in the case control data as well

cleaned_chek2_matched$lower_score<-cleaned_chek2_matched$refined_score-
  cleaned_chek2_matched$refined_score_se

cleaned_chek2_matched$upper_score<-cleaned_chek2_matched$refined_score+
  cleaned_chek2_matched$refined_score_se


range_checker2<-function(row){
  if(as.numeric(row[22])<0.5 && as.numeric(row[23])>0.5){
    return(1) #if the upper range and lower range cross 0.5 assign a value for subsetting
  }else{0} #score we are more sure about 
}

#for some reason the apply function did not work for what I was doing 
range_mask2<-vector("numeric",1818)
for(i in (1:1818)){
  row<-cleaned_chek2_matched[i,]
  range_mask2[i]<-range_checker2(row)
  
}

cleaned_chek2_matched<-cbind(cleaned_chek2_matched,range_mask2)

temp <- cleaned_chek2_matched %>% filter(range_mask2<1)

#sort the case control data based on refined score to assign them to a bin
sorted_chek2cc<-cleaned_chek2_matched[order(cleaned_chek2_matched$refined_score),]

#create a vector that represents the different bins in the sorted chek2cc
group_vector <-vector("numeric",1818)
group_counter<-0
for(i in 1:1818){
  group_vector[i]<-group_counter
  if(i%%92 == 0){ #modulo will determine the group size
    group_counter <- group_counter + 1
  }

}
#sorted_chek2cc$group_vector <- NULL
sorted_chek2cc<-cbind(sorted_chek2cc,group_vector)

write.csv(sorted_chek2cc,file="grouped_cc_chek2")

#look at if scores in some tissues are systematically higher than others 
#brca1 map to compare look at an example of sitatuon where tissue distribution is
#skewed
#look into more data about cancer cases with lower scores in tissue distribution
#(tcga,pcag....cancer genome atlas,cbioportal)

#coefficient of variation
