### Formatting the Bird RLI data so it can be fed into the model

require(dplyr)
require(tidyverse)
require(xtable)

### Call in the Bird RLI data
#Bird_RLI <- read.csv("../Data/Bird_RLI_original.csv", stringsAsFactors = T)

## Modified version of the assessment cleaner, renames categories so it'll assign correctly and removes sub-species
Assess_Clean_RLI <- function(Species_History){
  # Assigns an index for debugging
  Species_History <- Species_History %>% mutate(row_ID = row_number())
  # Remove all subspecies
  Species_Data <- Species_Data[which(is.na(Species_Data$infra_rank)),]
  Species_History <- Species_History[which(Species_History$taxonid %in% Species_Data$taxonid),]
  # Now rename codes where they have several names
  Species_History$category <- recode(Species_History$category, "Ex" = "EX", "Ex?" = "EX", "EW" = "EX", "LR/lc" = "LC", "LR/nt" = "NT", "LR/cd" = "NT", "CR(PE)" = "CR", "CR (PE)" = "CR", "CR(PEW)"= "CR")
  # Remove DD assessments
  Species_History <- Species_History[which(Species_History$category != "DD"),] 
  # Remove species with only one assessment remaining
  Species_History <- Species_History %>% group_by(taxonid) %>% filter(n()>1) %>% ungroup
  # Add species name to table
  Species_History <- inner_join(Species_History, Species_Data[,c(1,8)], by = "taxonid")
  Species_History$category<-droplevels(Species_History$category)
  return(Species_History)
}

Assess_Clean_RLI_PEXasEX <- function(Species_History){
  # Assigns an index for debugging
  Species_History <- Species_History %>% mutate(row_ID = row_number())
  # Remove all subspecies
  Species_Data <- Species_Data[which(is.na(Species_Data$infra_rank)),]
  Species_History <- Species_History[which(Species_History$taxonid %in% Species_Data$taxonid),]
  # Now rename codes where they have several names
  Species_History$category <- recode(Species_History$category, "Ex" = "EX", "Ex?" = "EX", "EW" = "EX", "LR/lc" = "LC", "LR/nt" = "NT", "LR/cd" = "NT", "CR(PE)" = "EX", "CR (PE)" = "CR", "CR(PEW)"= "CR")
  # Remove DD assessments
  Species_History <- Species_History[which(Species_History$category != "DD"),] 
  # Remove species with only one assessment remaining
  Species_History <- Species_History %>% group_by(taxonid) %>% filter(n()>1) %>% ungroup
  # Add species name to table
  Species_History <- inner_join(Species_History, Species_Data[,c(1,8)], by = "taxonid")
  Species_History$category<-droplevels(Species_History$category)
  return(Species_History)
}

Correct_False_Extinctions <- function(Species_History_Tags){
  # Any time where a species has a true extant category post extinction, 
  # that extinction should be labelled as false
  for (i in unique(Species_History_Tags$taxonid)){
    species <- Species_History_Tags[Species_History_Tags$taxonid == i,]
    if ("EX" %in% species$category) { #If there are any EX assessments
      EX_assess <- which(species$category == "EX")
      non_EX_assess <- which(species$category != "EX")
      if (length(non_EX_assess)>0){  
        # Determines EX assessments happening before extant ones
        False_assess <- species[EX_assess[min(non_EX_assess) > EX_assess],]
        if (nrow(False_assess)>0){
          # Assign tags if there are false de-extinctions
          Species_History_Tags[which(Species_History_Tags$row_ID %in% False_assess$row_ID),]$Verified <- "False"
        }
      } else if (length(non_EX_assess)==0){
        Species_History_Tags[which(Species_History_Tags$taxonid == i),]$Verified <- "False"
      }
    }
  }
  return(Species_History_Tags)
}

# Format the Bird RLI data for modelling
Format_BirdRLI <- function(Bird_RLI){
  Bird_RLI <- Bird_RLI[,2:10]
  # rename columns
  names(Bird_RLI) <- c("binomial", "taxonid", "1988", "1994", "2000","2004","2008","2012","2016")
  # pivot to long formatting
  Bird_RLI <- pivot_longer(Bird_RLI, cols = c("1988", "1994", "2000","2004","2008","2012","2016")
                           , names_to = "year", values_to = "category")
  Bird_RLI$year <- as.numeric(Bird_RLI$year)
  # Run the cleaner on it to remove subspecies and properly sort categories
  Bird_RLI <- Assess_Clean_RLI(Bird_RLI)
  Bird_RLI$Verified <- NA
  # Check for artificial de-extinctions and recode those assessments as CR. Also removes species with only EX assessments
  Bird_RLI <- Correct_False_Extinctions(Bird_RLI)
  Bird_RLI <- Bird_RLI[which(is.na(Bird_RLI$Verified)),]
  return(Bird_RLI)
}

Format_BirdRLI_PEX <- function(Bird_RLI){
  Bird_RLI <- Bird_RLI[,2:10]
  # rename columns
  names(Bird_RLI) <- c("binomial", "taxonid", "1988", "1994", "2000","2004","2008","2012","2016")
  # pivot to long formatting
  Bird_RLI <- pivot_longer(Bird_RLI, cols = c("1988", "1994", "2000","2004","2008","2012","2016")
                           , names_to = "year", values_to = "category")
  Bird_RLI$year <- as.numeric(Bird_RLI$year)
  # Run the cleaner on it to remove subspecies and properly sort categories
  Bird_RLI <- Assess_Clean_RLI_PEXasEX(Bird_RLI)
  Bird_RLI$Verified <- NA
  # Check for artificial de-extinctions and recode those assessments as CR. Also removes species with only EX assessments
  Bird_RLI <- Correct_False_Extinctions(Bird_RLI)
  Bird_RLI <- Bird_RLI[which(is.na(Bird_RLI$Verified)),]
  return(Bird_RLI)
}




