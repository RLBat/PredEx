require(dplyr)
require(tidyverse)

`%!in%` = Negate(`%in%`)

#### Read in your files here #####

Species_History <- read.csv("../Data/Species_History_IDs_20223.csv", stringsAsFactors = T)
Species_Data <- read.csv("../Data/Species_Data_20222.csv", stringsAsFactors = F)
Table7 <- read.csv("../Data/Table7.csv", header=TRUE, stringsAsFactors = FALSE)

###############################

Assess_Clean <- function(Species_History){
  # Assigns an index for debugging
  Species_History <- Species_History %>% mutate(row_ID = row_number())
  # Remove all subspecies
  Species_Data <- Species_Data[which(is.na(Species_Data$infra_rank)),]
  Species_History <- Species_History[which(Species_History$taxonid %in% Species_Data$taxonid),]
  # Remove all pre 1994 listings as that is when the current system was implemented
  Species_History <- subset(Species_History, Species_History$year >= 1994)
  # Select which codes to remove    
  lose_codes <- c("I","NR","K", "R", "CT", "E", "V", "Ex/E", "Ex?", "nt")
  # Remove those codes
  Species_History <- dplyr::filter(Species_History, !Species_History$category %in% lose_codes)
  # Now rename codes where they have several names
  Species_History$category <- recode(Species_History$category, "Ex" = "EX", "EW" = "EX", "LR/lc" = "LC", "LR/nt" = "NT", "LR/cd" = "NT")
  # Generate a df of only the years with two assessments
  Duplicates <- Species_History %>% group_by(taxonid) %>% filter(duplicated(year)|duplicated(year, fromLast=TRUE))
  # Remove these assessments as there is no way to know which order they were in
  Species_History <- Species_History %>% filter(!Species_History$row_ID %in% Duplicates$row_ID)
  # Remove species with only one assessment remaining
  Species_History <- Species_History %>% group_by(taxonid) %>% filter(n()>1) %>% ungroup
  # Add species name to table
  Species_History <- inner_join(Species_History, Species_Data[,c(1,8)], by = "taxonid")
  Species_History$category<-droplevels(Species_History$category)
  return(Species_History)
}

# read in Table7 csv
#Table7 <- read.csv("../Data/Table7.csv", header=TRUE, stringsAsFactors = FALSE)

#Species_History <- Assess_Clean(Species_History)

Add_table7_tags <- function(Table7, Species_History){
  Table7$year_sub <- Table7$year
  Table7$year <- floor(as.integer(readr::parse_number(as.character(Table7$year))))
  Species_History$scientific_name <- as.character(Species_History$scientific_name)
  Table7$scientific_name <- as.character(Table7$scientific_name)
  # Make a df of all assessment changes contained in Table 7
  ## This removes all species that have changed name
  Cat_Changes <- inner_join(Species_History, Table7[,c(1,3:6)], by = c("scientific_name", "year"), relationship = "many-to-many")
  # Cooerce to character
  Cat_Changes$category <- as.character(Cat_Changes$category)
  Cat_Changes$new_category <- as.character(Cat_Changes$new_category)
  Cat_Changes[Cat_Changes=="CR(PE)"|Cat_Changes=="CR(PEW)"|Cat_Changes=="CR (PE)"] <- "CR"
  # Make a df of non-matching assessments as a check
  check <- Cat_Changes %>% filter(category != Cat_Changes$new_category)
  test2 <- Cat_Changes %>% filter(row_ID %in% check$row_ID)# 10 that don't line up for no obvious reason
  # Remove those, assume clerical error
  Cat_Changes <- Cat_Changes %>% filter(!row_ID %in% check$row_ID)
  # Remove any non True/False reasons (E seems to stick around for some)
  Cat_Changes <- Cat_Changes[Cat_Changes$reason_for_change=="N"|Cat_Changes$reason_for_change=="G",]
  return(Cat_Changes)
}

#Cat_Changes <- Add_table7_tags(Table7, Species_History)

Same_cat_tag <- function(group_df){
  # Function to identify where two assignments in a row are the same and mark the older one as true
  for (i in 2:nrow(group_df)){
    if (group_df$category[i] == group_df$category[i-1]){
      group_df$Verified[i] <- "True"
    } else {
      next
    }
  }
  return(group_df)
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
        False_assess <- species[EX_assess[min(non_EX_assess) < EX_assess],]
        if (nrow(False_assess)>0){
          # Assign tags if there are false de-extinctions
          Species_History_Tags[which(Species_History_Tags$row_ID %in% False_assess$row_ID),]$Verified <- "False"
        }
      }
    }
  }
  return(Species_History_Tags)
}


Assign_known_tags <- function(Cat_Changes, Species_History){
  # Use Cat_changes to assign tags where we can, need to accces both assessments where a change has happened
  Species_History$Verified <- NA
  Species_History$Verified <- factor(Species_History$Verified, levels = c("Unknown", "True", "False"))
  # Do any changes that need to be done groupwise
  groupwise_df <- Species_History[NULL,]
  for (i in unique(Species_History$taxonid)){
    # subset to species
    species <- Species_History[(Species_History$taxonid == i),]
    # tag the newest entry as unknown
    species[(species$year==max(species$year)),]$Verified <- "Unknown"
    # Assign True tags where the assessment has been the same twice in a row
    species <- Same_cat_tag(species)
    # work out which years had DD classifications
    DD_years <- which(species$category=="DD")
    if (length(DD_years)>=1){
      # Assign from the DDyears backwards as false
      species[min(DD_years):nrow(species),]$Verified <- "False"
    }
    # Save to df
    groupwise_df <- rbind(groupwise_df, species)
  }
  Species_History <- groupwise_df
  # Create a reference df
  Reference <- Cat_Changes[,c("row_ID", "reason_for_change")]
  Reference$reason_for_change <- recode(Reference$reason_for_change, N = "False", G = "True")
  Reference <- Reference %>% mutate(row_ID = row_ID + 1) # Correct to be previous assessment
  # Add the True/False tags to the main df
  Species_History <- merge(Species_History, Reference, by = "row_ID", all.x = TRUE)
  x<-which(!is.na(Species_History$reason_for_change)) # select all rows with a true/false tag
  # Copy the tags to the main column
  Species_History$Verified[x] <- as.character(Species_History$reason_for_change[x])
  # Remove the column used to transfer the tags
  Species_History$reason_for_change <- NULL
  #### Tag all false extinctions (where it's later extant) as False
  Species_History <- Correct_False_Extinctions(Species_History)
  # Assign False tags to all DD assessments
  Species_History$Verified[which(Species_History$category=="DD")] <- "False" 
  # Assign unknown values to all other assessments
  Species_History$Verified[which(is.na(Species_History$Verified))] <- "Unknown"
  return(Species_History)
}

#Species_History <- Assign_known_tags(Cat_Changes, Species_History)

# checkpoint
# write.csv(Species_History, "../Data/SpeciesHistory_June21.csv", row.names = FALSE)
# Species_History <- read.csv("../Data/SpeciesHistory_June21.csv", header = T, stringsAsFactors = F)


Define_probabilities <- function(Species_History){
  # Probability of a change being true by cat # Very different! LC and EX are much more likely to be true
  Categories <- c("LC", "NT", "VU", "EN", "CR", "EX")
  Cat_probs <- c(rep(NA, length(Categories)))
  names(Cat_probs) <- Categories
  # Works out the prop(true) for each cat and assigns to named list
  for (i in Categories){
    Cat_subset <- Species_History[which(Species_History$category==i & Species_History$Verified!= "Unknown"),]
    Cat_probs[[i]] <- length(which(Cat_subset$Verified=="True"))/nrow(Cat_subset)
  }
  return(Cat_probs)
}

# Cat_probs <- Define_probabilities(Species_History)
# Last bit of deterministic stuff


#################################
####### Repeat this!!! #########
###############################

#### RANDOM BIT #####
Generate_tags <- function(Species_History, Cat_probs){
  # Randomly assigns T/F to unknown assessments based on category probabilities
  Categories <- c("LC", "NT", "VU", "EN", "CR", "EX")  
  for (i in Categories){
    # Finds the index values of all unknown assesments fora category
    Cat_unknown <- which(Species_History$category==i & Species_History$Verified=="Unknown")
    # Samples T/F based on the probability for that category
    Tags <- sample(c("True","False"), size=length(Cat_unknown), replace=TRUE, prob=c(Cat_probs[[i]], 1-Cat_probs[[i]]))
    # Places the T/F values where the unknown ones are
    Species_History$Verified[Cat_unknown] <- Tags
  }
  return(Species_History)
}

# Use different name to preserve the original df pre-random assignment
#Species_History_Tags <- Generate_tags(Species_History, Cat_probs)

Final_clean <- function(Corrected_cats){
  # Remove species with only EX assessments
  for (i in unique(Corrected_cats$taxonid)){
    species_cats <- unique(Corrected_cats[Corrected_cats$taxonid==i,]$category)
    if (length(species_cats)==1 && species_cats=="EX"){
      Corrected_cats[Corrected_cats$taxonid==i,]$Verified <- "False"
      }
  }
  # Any remaining as false must not meet the criteria so should be dropped
  Corrected_cats <- Corrected_cats[which(Corrected_cats$Verified!="False"),]
  # Remove species with only one assessment remaining
  Corrected_cats <- Corrected_cats %>% group_by(taxonid) %>% filter(n()>1) %>% ungroup
  return(Corrected_cats)
}

#Corrected_cats <- Final_clean(Species_History)

####################################

#Corrected_cats <- read.csv("../Data/Corrected_SpeciesHistory_June21.csv", header = T, stringsAsFactors = F)
#write.csv(Corrected_cats, "../Data/Corrected_SpeciesHistory_June222022.csv", row.names = FALSE)







