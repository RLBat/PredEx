###########################
#### Read in files ######

# Read in the table 7 data from 07-19
#Table7 <- read.csv("../Data/Table7.csv", header=TRUE, stringsAsFactors = FALSE)

set.seed(333)
setwd("Documents/PredEx/Code/")

###########################

###### DATA COLLECTION #####

# source("Data_Collection.R")
# 
# # Set API Key
# API_key = "" #use your personal key here
# 
# # Collect all species in the Red List and their base information
# Species_Data <- Species_Info_Collect(API_key = API_key)###### DATA COLLECTION #####

# 
# # Collect historical assessment data for each species
# Species_History<-Species_History_Collect()

### OR SKIP AND USE THE FOLLOWING ##

#Species_History <- read.csv("../Data/Species_History_IDs_20223.csv", stringsAsFactors = T)
#Species_Data <- read.csv("../Data/Species_Data_20222.csv", stringsAsFactors = F)

#########################

##### DATA PROCESSING #####

source("Data_Processing.R")

# Does initial data cleaning, removing invalid or unusable assessments
Species_History <- Assess_Clean(Species_History)

# Cleans up table 7 data for use
Cat_Changes <- Add_table7_tags(Table7, Species_History)

# Assigns any known tags
##### Takes a long time!!! #######
Species_History <- Assign_known_tags(Cat_Changes, Species_History)

# Works out the probabilities of an assessment being gen or non-gen
Cat_probs <- Define_probabilities(Species_History)

###############

### If running once, takes a good while
Species_History_Tags <- Generate_tags(Species_History, Cat_probs)

# Final cleaning steps
Corrected_cats <- Final_clean(Species_History_Tags)

#######################

######### MODELLING #########

source("Bootstrapping.R")

# Make Transition matrix
Q <- Transition_intensity_matrix(Categories <- c("LC", "NT", "VU", "EN", "CR", "EX"))

### Run the bootstrapped model
Boot_Probs <- Run_bootmarkov(Historic_assess = Corrected_cats, Q)

######### PLOTTING ##########

cats <- c("LC","NT","VU", "EN","CR", "EX")
Boot_median <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, median)
Boot_top <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.975)))
Boot_bottom <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.025)))

# Bind them together into one df for graphing
Boot_output <- bind_rows(Boot_median, Boot_bottom, Boot_top, .id = "Type")
Boot_output$Type[Boot_output$Type == 1] <- "Median"; Boot_output$Type[Boot_output$Type == 2] <- "Bottom"; Boot_output$Type[Boot_output$Type == 3] <- "Top"
Boot_output <- Boot_output[,1:7]
# Convert to long format
Boot_output <- gather(Boot_output, key = "Threat_level", value = "Probability", LC:CR)
Boot_output <- spread(Boot_output, key = "Type", value = "Probability")

Boot_output$Threat_level <- factor(Boot_output$Threat_level, levels = c("CR", "EN", "VU", "NT", "LC"))

######### Graphing ##############

p <- ggplot(data = Boot_output, aes(x = Time, y = Median, colour = Threat_level, xmax = 100)) + scale_color_manual(values = c("darkred", "darkorange", "gold", "darkcyan", "lightblue"))
p <- p + geom_line(size=1.2) + scale_y_continuous(breaks = seq(0,1,0.1), minor_breaks = seq(0,1,0.01))
p <- p + geom_ribbon(aes(ymin=Bottom, ymax=Top, alpha=0.5),fill="lightgrey", linetype = 2, show.legend = FALSE)
p <- p + labs(y = "Probability of extinction", x= "Time (years)", colour = "Red List Cateogry") 
p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
               axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), axis.title = element_text(size=20), 
               legend.position = c(0.2,0.8), legend.text = element_text(size=14), legend.title = element_text(size=16), strip.text = element_text(size=14))
p <- p + ylim(0,0.2)
p


###########################

########### SPECIES ATTRIBUTES  ############

#### TAXA ####


source("Taxa_processing.R")
# Create a vector of all species remaining after processing
Final_Species <- unique(Corrected_cats$taxonid)

Taxa_index<-Process_taxa(Species_Data, Final_Species)
# Removes misc if fewer than 500 species as it wouldn't be able to run. If more than 500 species,
# investigate to see what is in there (probably fungi)
if (length(Taxa_index$Misc) < 500){
  Taxa_index <- within(Taxa_index, rm(Misc))
} else {
  print("Large Misc category, investigation is advised.")
}

Species_wTaxa <- Assign_taxa(Corrected_cats, Taxa_index)

### Run a markov model for each species group, and non-species group.
### Some memory issues with the following
Taxa <- unique(na.omit(Species_wTaxa$Taxon))
taxonomic_models <- c()
not_taxonomic_models <- c()
for (i in 1:length(Taxa)){
  taxon <- Species_wTaxa[which(Species_wTaxa$Taxon==Taxa[i]),]
  not_taxon <- Species_wTaxa[which(Species_wTaxa$Taxon!=Taxa[i]),]
  taxonomic_models_i <- Run_bootmarkov(taxon, Q)
  taxonomic_models[[i]] <- taxonomic_models_i[which(taxonomic_models_i$Time == 100),]
  rm(taxonomic_models_i); gc()
  not_taxonomic_models_i <- Run_bootmarkov(not_taxon, Q)
  not_taxonomic_models[[i]] <- not_taxonomic_models_i[which(not_taxonomic_models_i$Time == 100),]
  rm(not_taxonomic_models_i); gc()
}

### get uncertainties

taxa_summary <- function(taxa, nottaxa){
  if (ncol(taxa) == 5){
    taxa$column_label <- 1:nrow(taxa)
    taxa <- taxa[,c(6,1,2,3,4,5)]
  }
  if (ncol(nottaxa) == 5){
    nottaxa$column_label <- 1:nrow(nottaxa)
    nottaxa <- nottaxa[,c(6,1,2,3,4,5)]
  }
  taxadb <- full_join(pivot_longer(taxa, cols = c(2:6), names_to = "Threat_level", values_to = "Probability"), pivot_longer(nottaxa, cols = c(2:6), names_to = "Threat_level", values_to = "Probability"), by = "column_label")
  taxadb <- taxadb[which(taxadb$Threat_level.x == taxadb$Threat_level.y),]
  taxadb$ratio <- log(taxadb$Probability.x/taxadb$Probability.y)
  Median <- taxadb %>% group_by(Threat_level.x) %>% summarise_at("ratio", median)
  Top <-  taxadb %>% group_by(Threat_level.x) %>% summarise_at("ratio", ~quantile(.x, c(.975)))
  Bottom <-  taxadb %>% group_by(Threat_level.x) %>% summarise_at("ratio", ~quantile(.x, c(.025)))
  taxa_summ <- gdata::combine(Median, Top, Bottom)
  names(taxa_summ) <- c("Threat_level", "Probability", "Source")
  return(taxa_summ)
}

# get significance values
p_values_taxa <- c()
Taxa <- unique(na.omit(Species_wTaxa$Taxon))[c(1,3:7)]
for (i in 1:(length(cats)-1)){
  p_values_taxa[i]<-t.test(mammal[, cats[i]], notmammal[, cats[i]])[[3]]
}
for (i in 1:(length(cats)-1)){
  p_values_taxa[i+5]<-t.test(fish[, cats[i]], notfish[, cats[i]])[[3]]
}
for (i in 1:(length(cats)-1)){
  p_values_taxa[i+10]<-t.test(invert[, cats[i]], notinvert[, cats[i]])[[3]]
}
for (i in 1:(length(cats)-1)){
  p_values_taxa[i+15]<-t.test(amphib[, cats[i]], notamphib[, cats[i]])[[3]]
}
for (i in 1:(length(cats)-1)){
  p_values_taxa[i+20]<-t.test(plant[, cats[i]], notplant[, cats[i]])[[3]]
}
for (i in 1:(length(cats)-1)){
  p_values_taxa[i+25]<-t.test(bird[, cats[i]], notbird[, cats[i]])[[3]]
}
p_values_taxa <- p.adjust(p_values_taxa, method = "bonferroni")

all_taxa <- matrix(p_values_taxa, nrow = 5)
all_taxa <- as.data.frame(all_taxa)
names(all_taxa) <- Taxa

## Get for each group
Mammals <- taxa_summary(mammal, notmammal)
Reptiles <- taxa_summary(reptile, notreptile)
Fish <- taxa_summary(fish, notfish)
Invertebrates <- taxa_summary(invert, notinvert)
Amphibians <- taxa_summary(amphib, notamphib)
Plants <- taxa_summary(plant, notplant)
Birds <- taxa_summary(bird, notbird)

#combine into one df for graphing
Taxa_comp <- gdata::combine(Mammals, Birds, Amphibians, Fish, Invertebrates, Plants)

# write.csv(Taxa_comp, file = "../Data/taxa_full.csv", row.names = FALSE)
Taxa_comp <- read.csv("../Data/taxa_full.csv", header = T)

catcols <- c("lightblue", "darkcyan", "gold", "darkorange2", "darkred")
Taxa_comp$Threat_level <- factor(Taxa_comp$Threat_level, levels = c("LC","NT", "VU", "EN", "CR"))
leg <- ggplot(data = Taxa_comp, aes(Probability, fill = Threat_level)) + geom_bar() +
  scale_fill_manual(values = catcols, name = "Red List Category", 
                    labels = c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered")) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 18))
leg <- cowplot::get_legend(leg)
p1 <- plot_taxa(threat_level = 1, xaxis = "", xlabels = element_blank()) + 
  annotate("text", x=3,y=2,label = "†")
p2 <- plot_taxa(threat_level = 2, xaxis = "", xlabels = element_blank(), yaxis = "") +
  annotate("text", x=3,y=2,label = "†")
p3 <- plot_taxa(threat_level = 3)
p4 <- plot_taxa(threat_level = 4, yaxis = "")
p5 <- plot_taxa(threat_level = 5, yaxis = "")

cairo_pdf(file = "../Output/Taxa_figure_140323.pdf", width = 12, height=8)
grid.arrange(grobs = list(p1, p2, p3, p4, p5, leg), layout_matrix = rbind(c(1,2, 6), c(3, 4,5)), heights = 4:5)
dev.off()

######################

#### BODY MASS ######

#####################

source("Body_Mass_Processing.R")
# Grab the species with enough data to model post-cleaning
Final_Species_List <- unique(Corrected_cats$scientific_name)


## MAMMALS ##


# List all mammals
iucnmammals <- Species_Data[which(Species_Data$class_name == "MAMMALIA"),"scientific_name"]
# Generate a list of all viable mammals
viable_mammals <- iucnmammals[which(iucnmammals %in% Final_Species_List)]
# import body mass data (which was manually taxon matched by Yuheng)
mammal_bm <- read.csv("../Data/mammal_bodymass_2022.csv")
# reduce body mass data down to viable species
matched_mammal <- mammal_bm[which(mammal_bm$IUCN_name %in% viable_mammals),c(1,3,5)]
# make a df with only the required species and add the body mass to each
Mammal_bm_assessments <- Historic_assess[which(Historic_assess$scientific_name %in% matched_mammal$IUCN_name),]
# rename cols for merge
names(matched_mammal) <- c("taxonid", "scientific_name", "body_mass")
# merge dfs to add body mass to the historic assessment data 
Mammal_bm_assessments <- merge(Mammal_bm_assessments, matched_mammal)


# split into sections for 2way comparison
mammal_median <- split_bm(mammal_bm, 2)
mammal_heavy <- Mammal_bm_assessments[which(Mammal_bm_assessments$body_mass > mammal_median),]
mammal_light <- Mammal_bm_assessments[which(Mammal_bm_assessments$body_mass < mammal_median),]

## Run the bootstrapped models for mammals
mammal_light_boot <- Run_bootmarkov(Historic_assess = mammal_light, Q)
mammal_heavy_boot <- Run_bootmarkov(Historic_assess = mammal_heavy, Q)

mammal_light_boot <- read.csv("../Data/Mammal_light_boot.csv", header = T)
mammal_heavy_boot <- read.csv("../Data/Mammal_heavy_boot.csv", header = T)

# extract values at 100 years
mammal_light_100 <- Boot_100(mammal_light_boot)
mammal_heavy_100 <- Boot_100(mammal_heavy_boot)

#rename columns for merge
mammal_heavy_100[,4] <- "Heavy"
mammal_light_100[,4] <- "Light"
# merge the dataframes
mammal_bm_100 <- rbind(mammal_heavy_100, mammal_light_100)
names(mammal_bm_100) <- c("Source", "Threat_level", "Probability", "Mass")
#birds_bm_100 <- birds_bm_100[which(birds_bm_100$Source == "Mean"),]
mammal_bm_100 <- mammal_bm_100 %>% mutate(Mass = fct_relevel(Mass, "Light", "Heavy"))

p_values_mammal <- c()
for (i in 1:(length(cats)-1)){
  p_values_mammal[i]<-t.test(mammal_heavy_boot[which(mammal_heavy_boot$Time == 100), cats[i]], 
    mammal_light_boot[which(mammal_light_boot$Time == 100), cats[i]])[[3]]
}

p1<-Plot_100(mammal_bm_100,ylabel = "P(EX) at t=100", xlabel = "", leg_pos = "none", plottag = "A") +
  annotate("text", x=5, y=0.24, label = "†")
leg <- Plot_100(mammal_bm_100) %>% get_legend() %>% as_ggplot()

## BIRDS ##


# import bird body mass data
birds_bm <- read.csv("../Data/Yuheng/viablebirds_mass_complete.csv")
birds_bm <- birds_bm[,c(2:5)]
names(birds_bm)[4] <- "body_mass"
# List all birds
iucnbirds <- Species_Data[which(Species_Data$class_name == "AVES"),"scientific_name"]
# Generate a list of all viable birds
viable_birds <- iucnbirds[which(iucnbirds %in% Final_Species_List)]
# reduce body mass data down to viable species
matched_birds <- birds_bm[which(birds_bm$IUCN_name %in% viable_birds),c(1,3,4)]
unmatched_birds <- birds_bm[which(birds_bm$IUCN_name %!in% viable_birds),c(1,3,4)]
# make a df with only the required species and add the body mass to each
bird_bm_assessments <- Historic_assess[which(Historic_assess$scientific_name %in% matched_birds$IUCN_name),]
# rename cols for merge
names(matched_birds) <- c("taxonid", "scientific_name", "body_mass")
# merge dfs to add body mass to the historic assessment data 
bird_bm_assessments <- merge(bird_bm_assessments, matched_birds)

### split into sections for comparisons ## 2 way
bird_median <- split_bm(birds_bm, 2)
birds_light <- bird_bm_assessments[which(bird_bm_assessments$body_mass < bird_median),]
birds_heavy <- bird_bm_assessments[which(bird_bm_assessments$body_mass > mammal_median),]

## Run the bootstrapped models for birds
bird_light_boot <- Run_bootmarkov(Historic_assess = birds_light, Q)
bird_heavy_boot <- Run_bootmarkov(Historic_assess = birds_heavy, Q)

bird_light_boot <- read.csv("../Data/Bird_light_boot.csv", header = T)
bird_heavy_boot <- read.csv("../Data/Bird_heavy_boot.csv", header = T)

## Extract values at 100 years
birds_light_100 <- Boot_100(bird_light_boot)
birds_heavy_100 <- Boot_100(bird_heavy_boot)

#rename columns for merge
birds_heavy_100[,4] <- "Heavy"
birds_light_100[,4] <- "Light"
# merge the dataframes
birds_bm_100 <- rbind(birds_heavy_100, birds_light_100)
names(birds_bm_100) <- c("Source", "Threat_level", "Probability", "Mass")
#birds_bm_100 <- birds_bm_100[which(birds_bm_100$Source == "Mean"),]
birds_bm_100 <- birds_bm_100 %>% mutate(Mass = fct_relevel(Mass, "Light", "Heavy"))

## run significance tests 

p_values_bird <- c()
for (i in 1:(length(cats)-1)){
  p_values_bird[i]<-t.test(bird_heavy_boot[which(bird_heavy_boot$Time == 100), cats[i]], 
    bird_light_boot[which(bird_light_boot$Time == 100), cats[i]])[[3]]
}

p2 <- Plot_100(birds_bm_100, xlabel = "", ylab = "", leg_pos = c(0.2,0.8), plottag = "B")

#### CR + EX ####

source("Model_to_CR.R")

## Mammals ##

## Run the bootstrapped models for mammals
mammal_light_boot2 <- Bootstrapped_probs_CRandEX(Historic_assess = mammal_light, Q)
mammal_heavy_boot2 <- Bootstrapped_probs_CRandEX(Historic_assess = mammal_heavy, Q)

mammal_light_boot2 <- read.csv("../Data/mammal_light_boot2_140323.csv", header = T)
mammal_heavy_boot2 <- read.csv("../Data/mammal_heavy_boot2_140323.csv", header = T)

# extract values at 100 years
mammal_light_100 <- Boot_100(mammal_light_boot2)
mammal_heavy_100 <- Boot_100(mammal_heavy_boot2)

#rename columns for merge
mammal_heavy_100[,4] <- "Heavy"
mammal_light_100[,4] <- "Light"
# merge the dataframes
mammal_bm_100 <- rbind(mammal_heavy_100, mammal_light_100)
names(mammal_bm_100) <- c("Source", "Threat_level", "Probability", "Mass")
mammal_bm_100 <- mammal_bm_100 %>% mutate(Mass = fct_relevel(Mass, "Light", "Heavy"))

p_values_mammalCREX <- c()
for (i in 1:(length(cats)-1)){
  p_values_mammalCREX[i]<-t.test(mammal_heavy_boot2[which(mammal_heavy_boot2$Time == 100), cats[i]], 
                      mammal_light_boot2[which(mammal_light_boot2$Time == 100), cats[i]])[[3]]
}

p3<-Plot_100(mammal_bm_100, ylabel = "P(CR or EX) at t=100", 
             xlabel = "Red List Category\nMammals", leg_pos = "none", y_limits = ylim(0,0.5), plottag = "C")

## Birds ##

## Run the bootstrapped models for birds
bird_light_boot2 <- Bootstrapped_probs_CRandEX(Historic_assess = birds_light, Q)
bird_heavy_boot2 <- Bootstrapped_probs_CRandEX(Historic_assess = birds_heavy, Q)

bird_light_boot2 <-read.csv("../Data/bird_light_boot2.csv", header = T)
bird_heavy_boot2 <- read.csv("../Data/bird_heavy_boot2.csv", header = T)

## Extract values at 100 years
birds_light_100 <- Boot_100(bird_light_boot2)
birds_heavy_100 <- Boot_100(bird_heavy_boot2)

#rename columns for merge
birds_heavy_100[,4] <- "Heavy"
birds_light_100[,4] <- "Light"
# merge the dataframes
birds_bm_100 <- rbind(birds_heavy_100, birds_light_100)
names(birds_bm_100) <- c("Source", "Threat_level", "Probability", "Mass")
#birds_bm_100 <- birds_bm_100[which(birds_bm_100$Source == "Mean"),]
birds_bm_100 <- birds_bm_100 %>% mutate(Mass = fct_relevel(Mass, "Light", "Heavy"))

## run significance tests 

p_values_birdCREX <- c()
for (i in 1:(length(cats)-1)){
  p_values_birdCREX[i]<-t.test(bird_heavy_boot2[which(bird_heavy_boot2$Time == 100), cats[i]], 
                      bird_light_boot2[which(bird_light_boot2$Time == 100), cats[i]])[[3]]
}

p_values_mass <- c(p_values_bird, p_values_birdCREX, p_values_mammal, p_values_mammalCREX)
p_values_mass <- p.adjust(p_values_mass, method = "bonferroni")

all_mass <- matrix(p_values_mass, nrow = 5)
all_mass <- as.data.frame(all_mass)
names(all_mass) <- c("Birds", "Birds CR/EX", "Mammals", "Mammals CR/EX")

p4<-Plot_100(birds_bm_100, ylabel = "", xlabel = "Red List Category\nBirds", leg_pos = "none", y_limits = ylim(0,0.5), plottag = "D")


cairo_pdf(file = "../Output/Body_Mass_figure_240124.pdf", width = 12, height=8)
grid.arrange(grobs = list(p1,p2,p3,p4), widths = c(2,2), layout_matrix = rbind(c(1,2), c(3,4)))

dev.off()

####################

### HABITAT ###

##################

source("Model_to_CR.R")

###########

#source("Habitat_data_collection.R")
All_Habitat <- read.csv("../Data/Habitat_Data.csv")
source("Habitat_data.R")

Habitat_data <- Process_habitat_data(All_Habitat)

Habitat_assess <- merge(Corrected_cats, Habitat_data, by.x = "scientific_name", by.y = "Species")

Habitat_general <- Habitat_assess[which(Habitat_assess$type == "Generalist"),]
Habitat_special <- Habitat_assess[which(Habitat_assess$type == "Specialist"),]

## Run the bootstrapped models for mammals
General_boot <- Run_bootmarkov(Historic_assess = Habitat_general, Q)
Special_boot <- Run_bootmarkov(Historic_assess = Habitat_special, Q)

#significance values
p_values_habitat <- c()
for (i in 1:(length(cats)-1)){
  p_values_habitat[i]<-t.test(General_boot[which(General_boot$Time == 100), cats[i]], 
                      Special_boot[which(Special_boot$Time == 100), cats[i]])[[3]]
}
p_values_habitat <- p.adjust(p_values_habitat, method = "bonferroni")

# extract values at 100 years
General_100 <- Boot_100(General_boot)
Special_100 <- Boot_100(Special_boot)

#rename columns for merge
General_100[,4] <- "Generalist"
Special_100[,4] <- "Specialist"
# merge the dataframes
Habitat_100 <- rbind(General_100, Special_100)
names(Habitat_100) <- c("Source", "Threat_level", "Probability", "Type")

Plot_habitat_100(Habitat_100)

### not in final paper

# Run for CR and EX
General_boot2 <- Bootstrapped_probs_CRandEX(Habitat_general, Q)
Special_boot2 <- Bootstrapped_probs_CRandEX(Habitat_special, Q)

# Extract values
General_1002 <- Boot_100(General_boot2)
Special_1002 <- Boot_100(Special_boot2)
# add column
General_1002[,4] <- "Generalist"
Special_1002[,4] <- "Specialist"

Habitat_1002 <- rbind(General_1002, Special_1002)
names(Habitat_1002) <- c("Source", "Threat_level", "Probability", "Type")

Plot_habitat_100(Habitat_1002, y_limits = ylim(0,0.5), ylabel = "P(EX or CR) at t=100")

################################

### Bird RLI analysis and comparisons ###
source("RLI2.R")

Bird_RLI <- read.csv("../Data/Bird_RLI_original.csv", header = T, stringsAsFactors = T)
Bird_RLI <- Format_BirdRLI(Bird_RLI)
Bird_RLI <- Bird_RLI[,c(5,2,3,4,6,7)]


Boot_BirdRLI <- Run_bootmarkov(Bird_RLI, Q)

Bird_RLI <- read.csv("../Data/Bird_RLI_original.csv", header = T, stringsAsFactors = T)

Bird_RLI_EX <- Format_BirdRLI_PEX(Bird_RLI)
Bird_RLI_EX <- Bird_RLI_EX[,c(5,2,3,4,6,7)]

Boot_BirdRLI_EX <- Run_bootmarkov(Bird_RLI_EX, Q)

#### Just birds to compare

# viable species
Final_Species_List <- unique(Corrected_cats$scientific_name)
# List all birds
iucnbirds <- Species_Data[which(Species_Data$class_name == "AVES"),"scientific_name"]
# Generate a list of all viable birds
viable_birds <- iucnbirds[which(iucnbirds %in% Final_Species_List)]
Bird_assessments <- Corrected_cats[which(Corrected_cats$scientific_name %in% viable_birds),]

Boot_birds <- Run_bootmarkov(Bird_assessments, Q)


source("Data_processing_noT7.R")
# grab all birds
iucnbirds <- Species_Data[which(Species_Data$class_name == "AVES"),"scientific_name"]
Final_Species_List <- unique(Species_History$scientific_name)
# Generate a list of all viable birds
viable_birds <- iucnbirds[which(iucnbirds %in% Final_Species_List)]
Bird_assessments_noT7 <- Species_History_noT7[which(Species_History_noT7$scientific_name %in% viable_birds),]

Boot_birds_notcorrected <- Run_bootmarkov(Bird_assessments_noT7, Q)

wT7 <- Boot_birds %>% Boot_100() %>% add_column(Data = "With T7")
woT7 <- Boot_birds_notcorrected %>% Boot_100() %>% add_column(Data = "Without T7")
RLICR <- Boot_BirdRLI %>% Boot_100() %>% add_column(Data = "RLI CR(PEX) as CR")
RLIEX <- Boot_BirdRLI_EX %>% Boot_100() %>% add_column(Data = "RLI CR(PEX) as EX")

Birds100 <- rbind(wT7, woT7, RLICR, RLIEX)
Birds100$Data <- factor(Birds100$Data, levels = c("With T7", "Without T7", "RLI CR(PEX) as CR", "RLI CR(PEX) as EX"))

Birds100$Threat_level <- recode(Birds100$Threat_level, "1" = "LC", "2" = "NT", "3" = "VU", "4" = "EN", "5" = "CR")

Plot_100_birds <- function(hundred_year){
  p <- ggplot(data = subset(hundred_year, Source %in% c("Median")), aes(x = Threat_level, y = Probability, fill = Data)) + scale_fill_manual(values = c("cyan3", "tomato3", "purple", "goldenrod"), name = "Body Mass")
  p <- p + geom_bar(stat = "identity", position = "dodge") #+ scale_x_discrete(breaks = 1:5, labels=c("LC","NT","VU", "EN","CR"))
  p <- p + labs(y = "Probability of Ex or CR at t=100", x = "IUCN Species Threat Assessment")
  p <- p + geom_errorbar(aes(ymin= hundred_year$Probability[hundred_year$Source == "Bottom"], ymax=hundred_year$Probability[hundred_year$Source == "Top"]), width=.2, position=position_dodge(.9)) 
  p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
                 axis.text.y = element_text(size=16), axis.title = element_text(size=20), axis.text.x = element_text(size=16), legend.position = c(0.2,0.8), legend.text = element_text(size=12), legend.title = element_text(size=14), strip.text = element_text(size=14))
  return(p)
}

Plot_100_birds(Birds100)

#####

# ANOVA testing

wT7 <- Boot_birds %>% add_column(Data = "With T7")
woT7 <- Boot_birds_notcorrected %>% add_column(Data = "Without T7")
RLICR <- Boot_BirdRLI %>% add_column(Data = "RLI CR(PEX) as CR")
RLIEX <- Boot_BirdRLI_EX %>% add_column(Data = "RLI CR(PEX) as EX")
Birds_all100 <- rbind(wT7, woT7, RLICR, RLIEX)

Birds_all100 <- Birds_all100[Birds_all100$Time==100,c(3:7,9)]

kruskal.test(LC ~ Data, data=Birds_all100)
pairwise.wilcox.test()


t.test()

read.csv()

### Checking for EW

# read in the EW as CR booted values
EWasCR <- read.csv("../Data/BootProbs_EWasCR.csv", header = T, stringsAsFactors = T)
EWasCR <- EWasCR[which(EWasCR$Time==100),]
EWasEX <- read.csv("../Data/Overall_Probabilities.csv", header = T, stringsAsFactors = T)
EWasEX <- EWasEX[which(EWasEX$Time==100),]
EWasEX$Treatment <- "EW_as_EX"
EWasCR$Treatment <- "EW_as_CR"

EW_testing <- rbind(EWasCR[,c(2:9)], EWasEX[,c(2:9)])

wilcox.test(EWasCR$LC, EWasEX$LC)
wilcox.test(EWasCR$NT, EWasEX$NT)
wilcox.test(EWasCR$VU, EWasEX$VU)
wilcox.test(EWasCR$EN, EWasEX$EN)
wilcox.test(EWasCR$CR, EWasEX$CR)

#####################################################################

### sandbox TO CLEAN ####

# Checking the Bird body mass data for CR + EX to work out the confusing outcomes with few EX examples

light_probs <- Bootstrapped_probs_CRandEX(birds_light, Q)
heavy_probs <- Bootstrapped_probs_CRandEX(birds_heavy, Q)

placeholder <- function(Boot_Probs){
  Boot_median <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, median)
  Boot_top <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.975)))
  Boot_bottom <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.025)))
  
  # Bind them together into one df for graphing
  Boot_output <- bind_rows(Boot_median, Boot_bottom, Boot_top, .id = "Type")
  Boot_output$Type[Boot_output$Type == 1] <- "Median"; Boot_output$Type[Boot_output$Type == 2] <- "Bottom"; Boot_output$Type[Boot_output$Type == 3] <- "Top"
  Boot_output <- Boot_output[,1:7]
  # Convert to long format
  Boot_output <- gather(Boot_output, key = "Threat_level", value = "Probability", LC:CR)
  Boot_output$Threat_level <- factor(Boot_output$Threat_level, levels = c("CR", "EN", "VU", "NT", "LC"))
  # only at 100 years
  Boot_output <- Boot_output[which(Boot_output$Time==100),]
  Boot_output <- within(Boot_output, rm(Time))
  return(Boot_output)
}

light_output <- placeholder(light_probs)
heavy_output <- placeholder(heavy_probs)

##### from body_mass-processing 
birds_heavy_100[,4] <- "Heavy"
birds_light_100[,4] <- "Light"

birds_bm_100 <- rbind(birds_heavy_100, birds_light_100)
names(birds_bm_100) <- c("Source", "Threat_level", "Probability", "V4")

birds_bm_100 <- birds_bm_100 %>% mutate(V4 = fct_relevel(V4, "Light", "Heavy"))
birds_bm_100$Threat_level <-factor(birds_bm_100$Threat_level, levels = cats) 

Plot_100 <- function(hundred_year){
  p <- ggplot(data = subset(hundred_year, Source %in% c("Median")), aes(x = Threat_level, y = Probability, fill = V4)) + scale_fill_manual(values = c("cyan3", "tomato3", "purple"), name = "Body Mass")
  p <- p + geom_bar(stat = "identity", position = "dodge")
  p <- p + labs(y = "Probability of Ex or CR at t=100", x = "IUCN Species Threat Assessment")
  p <- p + geom_errorbar(aes(ymin= hundred_year$Probability[hundred_year$Source == "Bottom"], ymax=hundred_year$Probability[hundred_year$Source == "Top"]), width=.2, position=position_dodge(.9)) 
  p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
                 axis.text.y = element_text(size=16), axis.title = element_text(size=20), axis.text.x = element_text(size=16), legend.position = c(0.2,0.8), legend.text = element_text(size=12), legend.title = element_text(size=14), strip.text = element_text(size=14))
  return(p)
}

Plot_100(birds_bm_100)



## 3 way
bird_bottom_100[,4] <- "Bottom"
bird_middle_100[,4] <- "Middle"
bird_top_100[,4] <- "Top"


#### Supp mats: redo for 3 way

light_probs <- Bootstrapped_probs_CRandEX(birds_light, Q)
bottom_probs <- Bootstrapped_probs_CRandEX(birds_bottom, Q)
middle_probs <- Bootstrapped_probs_CRandEX(birds_middle, Q)

### save as birds_x_100

birds_bm_100 <- rbind(bird_bottom_100, bird_middle_100, bird_top_100)
names(birds_bm_100) <- c("Source", "Threat_level", "Probability", "V4")
birds_bm_100 <- birds_bm_100 %>% mutate(V4 = fct_relevel(V4, "Bottom", "Middle", "Top"))
birds_bm_100$Threat_level <-factor(birds_bm_100$Threat_level, levels = cats) 
Plot_100(birds_bm_100)

#####################

# Extinction risk comparison plot

CritE <- c(0.0001, 0.01, 0.1, 0.667, 0.999)
RedList <- c(0,0.2,0.4,0.6,0.8)

Overall_probs <- read.csv("../Data/Overallbootoutput.csv", stringsAsFactors = T)
Overall_probs <- Overall_probs[Overall_probs$Time==100,c("Threat_level", "Median")]
Overall_probs <- Overall_probs %>% mutate(Threat_level = factor(Threat_level, levels = cats)) %>% arrange((Threat_level))
Bates <- Overall_probs
RedList_Probs <- data.frame(Bates, CritE, RedList)
names(RedList_Probs) <- c("Category", "Bates", "Criterion E", "Red List Index")

RedList_Probs <- pivot_longer(RedList_Probs, cols = c("Bates", "Criterion E", "Red List Index"), 
                              names_to = "Method", values_to = "Probability")

ggplot(RedList_Probs, aes(x=Category, y=Probability, group = Method)) +
  geom_line(aes(colour = Method), linewidth = 2) + scale_color_manual(values = c("red", "orange", "darkblue")) +
  theme(panel.grid.major = element_blank(), panel.background = element_blank(), 
        panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), 
        axis.line.x = element_line(colour = "black"), axis.text.y = element_text(size=16), 
        axis.title = element_text(size=20), axis.text.x = element_text(size=16), 
        legend.position = c(0.3,0.8), legend.text = element_text(size=12), 
        legend.title = element_text(size=14), strip.text = element_text(size=14))



#######################






