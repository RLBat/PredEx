require(ggplot2)
require(gdata)
require(gridExtra)
require(cowplot)

Corrected_cats <- read.csv("../Data/Corrected_SpeciesHistory_June222022.csv", header = T, stringsAsFactors = F)
#Species_Data <- read.csv("../Data/Species_Data.csv", stringsAsFactors = F)

########################

## Uncertainty function_

# Species_list <- unique(Species_History$taxonid)
# 
# Uncertainty <- function(Species_History, Species_list){
#   # Take a random half of the species
#   sub <- sample(Species_list, size = length(Species_list)/2)
#   # subset the overall df
#   sub_species <- subset(Species_History, Species_History$taxonid %in% sub)
#   # Randomly generate tags for subset
#   sub_species <- Generate_tags(sub_species, Cat_probs)
#   # Clean!
#   sub_species <- Reassign_Cats(sub_species)
#   sub_species <- Final_clean(sub_species)
#   return(sub_species)
# }


####################################

## Threat and Taxon processing


## Use Species_Data to subset by taxon.
Filter_taxa <-function(Species){
  classification <- NA
  if (Species["kingdom_name"] == "PLANTAE"){
    classification <- "Plant"
  } else if (Species["kingdom_name"] != "PLANTAE" && Species["kingdom_name"] != "ANIMALIA"){
    classification <- "Misc"
  } else if (Species["phylum_name"] != "CHORDATA" && Species["kingdom_name"] == "ANIMALIA"){
    classification <- "Invertebrate"    
  } else if (Species["class_name"] == "AMPHIBIA"){
    classification <- "Amphibian"
  } else if (Species["class_name"] == "AVES"){
    classification <- "Bird"
  } else if (Species["class_name"] == "MAMMALIA"){
    classification <- "Mammal"
  } else if (Species["class_name"] == "REPTILIA"){
    classification <- "Reptile"
  } else {
    classification <- "Fish"
  }
  return(classification)
}

# Filter to highest sensible taxa and get ids for each
Process_taxa <- function(Species_Data, Final_Species){
  # Filter the data down to species that are included in the modelling
  Species_Data <- dplyr::filter(Species_Data, taxonid %in% Final_Species)
  # Label each species with its highest taxon
  Species_Data$Taxon <- apply(Species_Data, 1, Filter_taxa)
  # List of used taxa
  Taxa <- c("Plant", "Misc", "Invertebrate", "Amphibian", "Bird", "Mammal", "Reptile", "Fish")
  # Get ids for each used taxa
  Taxa_index <- lapply(Taxa, function(i) unique(Species_Data$taxonid[which(Species_Data$Taxon==i)])) 
  # rename lists
  names(Taxa_index) <- Taxa
  return(Taxa_index)
}

Assign_taxa <- function(Species_History, Taxa_index){
  Species_History$Taxon <- NA
  for (i in 1:length(Taxa_index)){
    #Taxon_species <- which(Species_History$taxonid %in% Taxa_index[[i]])
    Species_History$Taxon[Species_History$taxonid %in% Taxa_index[[i]]] <- names(Taxa_index)[i]
  }
  return(Species_History)
}

plot_taxa <- function (threat_level, xaxis = "Taxonomic Group", xlabels = element_text(size=16, angle = 90,  hjust = 0.8, vjust = 0.5), yaxis = "Comparative Extinction Risk", ylabels = element_text(size=16)){
  i <- threat_level
  cats <- c("LC","NT","VU", "EN","CR")
  catcols <- c("lightblue", "darkcyan", "goldenrod1", "darkorange", "darkred")
  fillcolour <-catcols[i]
  categories <- c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered")
  p <- ggplot(data = subset(Taxa_comp, Source %in% c("Median") & Threat_level %in% cats[i]), aes(x = source, y = Probability))
  p <- p + geom_bar(stat = "identity", position = "dodge", fill = fillcolour) + #scale_x_discrete(labels=c("Mammal", "Reptile", "Fish", "Invertebrate", "Amphibian", "Plant", "Bird")) +
    labs(x = xaxis, tag = "", y = yaxis) + ylim(-4,3)
  #+  scale_fill_manual(values = c("lightblue", "darkcyan", "orange", "darkorange", "orangered3", "darkred"))
  p <- p + geom_errorbar(aes(ymin= Taxa_comp$Probability[Taxa_comp$Source == "Bottom" & Taxa_comp$Threat_level == cats[i]], ymax=Taxa_comp$Probability[Taxa_comp$Source == "Top" & Taxa_comp$Threat_level == cats[i]]), width=.2, position=position_dodge(.9)) 
  p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), 
                 axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
                 axis.text.y = ylabels, axis.title = element_text(size=16), axis.text.x = xlabels, axis.ticks.x = element_blank(),
                 legend.title = element_text(size=24), legend.text = element_text(size=20), strip.text = element_text(size=14), plot.tag.position = "top", plot.tag = element_text(size = 14))+
    geom_hline(yintercept = 0)
  p
  return(p)
}





