require(dplyr)
require(tidyverse)
require(ggplot2)
require(ggpubr)
require(gridExtra)


`%!in%` = Negate(`%in%`)

set.seed(333)

#######################################

Species_Data <- read.csv("../Data/Species_Data_20222.csv")
Species <- read.csv("../Data/Corrected_SpeciesHistory_June222022.csv", header=T, stringsAsFactors = F)
Historic_assess <- read.csv("../Data/Corrected_SpeciesHistory_June222022.csv", header = T, stringsAsFactors = F)

# Grab the species with enough data to model post-cleaning
#Final_Species_List <- Historic_assess$scientific_name

split_bm <- function(bm_data, split_no = 3){
  # initialise vecotr
  splits <- c()
  for(i in 1:split_no){
    # create a vector of the splits
    splits <- append(splits, i/split_no)
  }
  #remove any that are 0 or one to avoid errors
  splits <- splits[! splits %in% c(0,1)]
  
  #threshold for the split
  bm_threshold <- quantile(bm_data$body_mass, splits)
  
  return(bm_threshold)
}

########### MAMMALS ######

# mammal_thirds <- split_bm(mammal_bm, 3)
# # split into 3
# mammal_bottom <- Mammal_bm_assessments[which(Mammal_bm_assessments$body_mass < mammal_thirds[1]),]
# mammal_middle <- Mammal_bm_assessments[which(Mammal_bm_assessments$body_mass > mammal_thirds[1] & Mammal_bm_assessments$body_mass < mammal_thirds[2]),]
# mammal_top <- Mammal_bm_assessments[which(Mammal_bm_assessments$body_mass > mammal_thirds[2]),]
# 
# ############ BIRDS ##########
# 
# bird_third <- split_bm(birds_bm, 3)
# # split into three
# birds_bottom <- bird_bm_assessments[which(bird_bm_assessments$body_mass < bird_third[1]),]
# birds_middle <- bird_bm_assessments[which(bird_bm_assessments$body_mass > bird_third[1] & bird_bm_assessments$body_mass < bird_third[2]),]
# birds_top <- bird_bm_assessments[which(bird_bm_assessments$body_mass > bird_third[2]),]
# 
# ###############################
# 
# ## Run model as per the workflow ##
# bird_heavy_boot <- read.csv("../Data/Bird_heavy_boot.csv")
# bird_light_boot <- read.csv("../Data/Bird_light_boot.csv")
# mammal_heavy_boot <- read.csv("../Data/Mammal_heavy_boot.csv")
# mammal_light_boot <- read.csv("../Data/Mammal_light_boot.csv")

####

Boot_100 <- function(Boot_Probs){
  cats <- c("LC","NT","VU", "EN","CR", "EX")
  
  Boot_means <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, mean)
  Boot_top <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.975)))
  Boot_bottom <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.025)))
  
  hundred_year <- rbind(Boot_means[100,2:6], Boot_bottom[100,2:6], Boot_top[100,2:6])
  hundred_year["Source"] <- c("Mean", "Bottom", "Top")
  hundred_year <- as.data.frame(hundred_year)
  hundred_year <- gather(hundred_year, key = "Threat_level", value = "Probability", 1:5)
  hundred_year$Threat_level <- recode(hundred_year$Threat_level, "LC" = 1, "NT" = 2, "VU" = 3, "EN" = 4, "CR" = 5)
  hundred_year$Threat_level <- as.factor(hundred_year$Threat_level)
  hundred_year$Probability <- as.numeric(hundred_year$Probability)
  return(hundred_year)
}


### graphing ###

Plot_100 <- function(hundred_year, ylabel = "Probability of extinction at t=100", xlabel = "IUCN Species Threat Assessment", leg_pos = c(0.2,0.8), y_limits = ylim(0,0.25), plottag = ""){
  p <- ggplot(data = subset(hundred_year, Source %in% c("Mean")), aes(x = Threat_level, y = Probability, fill = Mass)) + scale_fill_manual(values = c("darkcyan", "darkorange", "darkred"), name = "Body Mass")
  p <- p + geom_bar(stat = "identity", position = "dodge") + scale_x_discrete(breaks = 1:5, labels=c("LC","NT","VU", "EN","CR")) + y_limits
  p <- p + labs(y = ylabel, x = xlabel, tag = plottag)
  p <- p + geom_errorbar(aes(ymin= hundred_year$Probability[hundred_year$Source == "Bottom"], ymax=hundred_year$Probability[hundred_year$Source == "Top"]), width=.2, position=position_dodge(.9)) 
  p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
                 axis.text.y = element_text(size=16), axis.title = element_text(size=20), axis.text.x = element_text(size=16), legend.position = leg_pos, legend.text = element_text(size=12), 
                 legend.title = element_text(size=14), strip.text = element_text(size=14), plot.tag.position = "topleft")
  return(p)
}

#Plot_100(birds_bm_100)

#t.test(Bird_Heavy_boot[which(Bird_Heavy_boot$Time == 100), "EN"], Bird_Light_boot[which(Bird_Light_boot$Time == 100), "EN"])

significance_test_2way <- function(heavy, light, cats = c("LC","NT","VU", "EN","CR", "EX")){
  heavy <- heavy[which(heavy$Time==100),]
  heavy[,ncol(heavy)+1] <- "Heavy"
  light <- light[which(light$Time==100),]
  light[,ncol(light)+1] <- "Light"
  bodymass <- rbind(heavy, light)
  lapply(bodymass, function(x) t.test(x~bodymass$Mass)[c("LC","NT","VU", "EN","CR")])
  test <- bodymass %>% map_df(~ t.test(. ~ bodymass$Mass), .id = 'var')
  t.test(light$LC, heavy$LC)
  t.test(light$NT, heavy$NT)
  t.test(light$VU, heavy$VU)
  t.test(light$EN, heavy$EN)
  t.test(light$CR, heavy$CR)
}



### 3 way split 

# bird_bottom_100 <- Boot_100(bird_bottom_boot)
# bird_middle_100 <- Boot_100(bird_middle_boot)
# bird_top_100 <- Boot_100(bird_top_boot)
# 
# bird_bottom_100[,4] <- "Bottom"
# bird_middle_100[,4] <- "Middle"
# bird_top_100[,4] <- "Top"
# 
# birds_bm_100 <- rbind(bird_bottom_100, bird_middle_100, bird_top_100)
# Plot_100(birds_bm_100)

# 
# 
# Mammal_bm_100$
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# Final_Species_List <- Historic_assess$scientific_name
# 
# ## read in a list of viable species to save processing time
# 
# iucnmammals <- Species_Data[which(Species_Data$class_name == "MAMMALIA"),"scientific_name"]
# ## reduce this list to only viable species - i.e. ones that pass the cleaning steps.
# 
# mammal_bm <- read.csv("../Data/mammal_bodymass_2022.csv")
# 
# matched_mammal <- mammal_bm[which(mammal_bm$IUCN_name %in% iucnmammals),]
# 
# matched_mammal <- matched_mammal[,c(1,3,5)]


### checked and corrected issues with out of date matches
# unmatched_mammal <- mammal_bm[which(mammal_bm$IUCN_name %!in% iucnmammals),]
# mammal_bm$IUCN_name[2394] <- "Micronomus norfolkensis"
# mammal_bm$IUCN_name[3545] <- "Desmalopex leucopterus"
# mammal_bm$IUCN_name[4012] <- "Sorex monticola"
# mammal_bm$IUCN_name[4200] <- "Tamiops mcclellandii"
# mammal_bm$IUCN_name[4222] <- "Cephalopachus bancanus"
# 
# write.csv(mammal_bm, "../Data/mammal_bodymass_2022.csv", row.names = FALSE)
