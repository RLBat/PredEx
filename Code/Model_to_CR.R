## Import packages
require(msm)
require(dplyr)
require(doParallel)
require(tidyverse)

`%!in%` = Negate(`%in%`)

############################

#### realised this doesn't actually work

# Transition_intensity_matrix_CR <- function(Categories = c("LC", "NT", "VU", "EN", "CR")){
#   # Build transition intensity matrix
#   # This defines that they have to pass through all states to reach death
#   # Slightly increases extinction chance but not much
#   Q <- rbind(c(0.5, 0.5, 0, 0, 0),
#              c(0.25, 0.5, 0.25, 0, 0),
#              c(0, 0.25, 0.5, 0.25, 0),
#              c(0, 0, 0.25, 0.5, 0.25),
#              c(0, 0, 0, 0.5, 0.5))
#   row.names(Q) <- Categories
#   colnames(Q) <- Categories
#   #Q[6,6] <- 0 #removes the transition from extinct to extinct
#   return(Q)
# }
# 
# Q <- Transition_intensity_matrix_CR()

## Remove any species with EX assessments
# Species_noEX <- Species_History %>% group_by(taxonid) %>% filter(!any(category == "EX")) %>% ungroup
# 
# Run_Markov <- function(Historic_assess, Q){
#   # Reverse year order for the state table
#   Historic_assess$category <- as.character(Historic_assess$category)
#   Historic_assess<-arrange(Historic_assess, taxonid, year)
#   
#   # Change cats to numbers for modelling
#   Historic_assess$category <- recode(Historic_assess$category, "LC" = "1", "NT" = "2", "VU" = "3", "EN"= "4", "CR" = "5", "EX" = "6")
#   Historic_assess$category <- as.integer(Historic_assess$category)
#   
#   # Create state table
#   state_table <- statetable.msm(state= Historic_assess$category, subject = Historic_assess$taxonid)
#   
#   # Add a Time since first assessment column (TSFA)
#   Historic_assess <- Historic_assess %>% group_by(taxonid) %>% mutate(TSFA = year-min(year))
#   
#   # Actual model!
#   msm_model <- msm(category ~ TSFA, subject = taxonid, data = Historic_assess, qmatrix = Q, gen.inits = TRUE, control=list(fnscale=60000,maxit=500))# covariates = ~ gen_count)
#   
#   return(msm_model)
# }
# 
# msm_noEX <- Run_Markov(Species_noEX, Q)

#################################

## trying new approach, ignore above

#Historic_assess <- read.csv("../Data/Corrected_SpeciesHistory_June222022.csv", header = TRUE)

#Q <- Transition_intensity_matrix(Categories <- c("LC", "NT", "VU", "EN", "CR", "EX"))

#probs_noEX <- Run_bootmarkov(Species_History, Q, years = 1:100, state = "CR")

Bootstrapped_probs_CRandEX <- function(Historic_assess, Q, years = 1:100){
  ### Overall model
  msm_model <- Run_Markov(Historic_assess, Q)
  
  ## Bootstrap the model
  Boot_models <- Bootstrap_msm(msm_model, repeats = 100)
  
  ### If not all converge
  Boot_models <- Boot_models[!sapply(Boot_models, function(x) class(x) == "try-error")]
  
  if (length(Boot_models)<100){
    print("Not all models converged")
  }
  
  cats <- c("LC","NT","VU", "EN","CR", "EX")
  #state = "CR"
  # Extract the probabilities to save or graph
  Boot_Probs_CR <- Boot_probs(Boot_models = Boot_models, years = years, state = "CR")
  Boot_Probs_EX <- Boot_probs(Boot_models = Boot_models, years = years, state = "EX")
  Boot_Probs_Sum <- Boot_Probs_CR
  Boot_Probs_Sum[,cats] <- Boot_Probs_CR[,cats] + Boot_Probs_EX[,cats]
  
  ## code to very quickly check the outcome of models
  test <- Boot_Probs_EX[which(Boot_Probs_EX$Time == max(years)),]
  summ <- test %>% summarise_all(range)
  if (0 %in% summ[1,]){
    print("Some models produced a 0% chance of extinction")
  }
  return(Boot_Probs_Sum)
}

#write.csv(Boot_Probs_Sum, "../Data/Birds_heavy_CRandEX_100", row.names = FALSE)






######### PLOTTING ##########




















