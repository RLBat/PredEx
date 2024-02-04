## Import packages
require(msm)
require(dplyr)
require(doParallel)
require(tidyverse)

`%!in%` = Negate(`%in%`)

set.seed(333)

############################

Transition_intensity_matrix <- function(Categories = c("LC", "NT", "VU", "EN", "CR", "EX")){
  # Build transition intensity matrix
  # This defines that they have to pass through all states to reach death
  # Slightly increases extinction chance but not much
  Q <- rbind(c(0.5, 0.5, 0, 0, 0, 0),
             c(0.25, 0.5, 0.25, 0, 0, 0),
             c(0, 0.25, 0.5, 0.25, 0, 0),
             c(0, 0, 0.25, 0.5, 0.25, 0),
             c(0, 0, 0, 0.25, 0.5, 0.25),
             c(0, 0, 0, 0, 0, 1))
  row.names(Q) <- Categories
  colnames(Q) <- Categories
  Q[6,6] <- 0 #removes the transition from extinct to extinct
  return(Q)
}

Run_Markov <- function(Historic_assess, Q){
  # Reverse year order for the state table
  Historic_assess$category <- as.character(Historic_assess$category)
  Historic_assess<-arrange(Historic_assess, taxonid, year)
  
  # Change cats to numbers for modelling
  Historic_assess$category <- recode(Historic_assess$category, "LC" = "1", "NT" = "2", "VU" = "3", "EN"= "4", "CR" = "5", "EX" = "6")
  Historic_assess$category <- as.integer(Historic_assess$category)
  
  # Create state table
  state_table <- statetable.msm(state= Historic_assess$category, subject = Historic_assess$taxonid)
  
  # Add a Time since first assessment column (TSFA)
  Historic_assess <- Historic_assess %>% group_by(taxonid) %>% mutate(TSFA = year-min(year))
  
  # Actual model!
  msm_model <- msm(category ~ TSFA, subject = taxonid, data = Historic_assess, qmatrix = Q, gen.inits = TRUE, control=list(fnscale=60000,maxit=500))# covariates = ~ gen_count)
  
  return(msm_model)
}

Bootstrap_msm <- function(msm_model, repeats = 100){
  # The bootstrapping, this takes a while (~15 mins)
  Boot_models <- boot.msm(msm_model, stat = NULL, B=repeats, cores = (detectCores()-1))
  ### Add an if function to check for non-converged models
  # Boot_models <- Boot_models[!sapply(Boot_models, function(x) class(x) == "try-error")]
  return(Boot_models)
}

Extract_probs <- function(model = Boot_models, years = 1:100, state = "EX"){
  Probabilities <- data.frame(Time = years, LC = NA, NT = NA, VU = NA, EN = NA, CR = NA, EX = NA)
  for (i in years){
    Probabilities[i,2:7] <- pmatrix.msm(model, t=i)[,state]
  }
  return(Probabilities)
}

Boot_probs <- function(Boot_models, years = c(1:100), state = "EX"){
  Boot_probabilities<-lapply(Boot_models, Extract_probs, years = years, state = state)
  Boot_probabilities <- bind_rows(Boot_probabilities, .id = "column_label")
  return(Boot_probabilities)
}


Run_bootmarkov <- function(Historic_assess, Q, years = 1:100, state = "EX"){
  ### Overall model
  msm_model <- Run_Markov(Historic_assess, Q)
  
  ## Bootstrap the model
  Boot_models <- Bootstrap_msm(msm_model, repeats = 100)
  
  ### If not all converge
  Boot_models <- Boot_models[!sapply(Boot_models, function(x) class(x) == "try-error")]
  
  if (length(Boot_models)<100){
    print("Not all models converged")
  }
  
  # Extract the probabilities to save or graph
  Boot_Probs <- Boot_probs(Boot_models = Boot_models, years = years, state = state)
  
  ## code to very quickly check the outcome of models
  test <- Boot_Probs[which(Boot_Probs$Time == max(years)),]
  summ <- test %>% summarise_all(range)
  if (0 %in% summ[1,]){
    print("Some models produced a 0% chance of extinction")
  }
  return(Boot_Probs)
}




Boot_100 <- function(Boot_Probs){
  cats <- c("LC","NT","VU", "EN","CR", "EX")
  
  Boot_median <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, median)
  Boot_top <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.975)))
  Boot_bottom <- Boot_Probs %>% group_by(Time) %>% summarise_at(cats, ~quantile(.x, c(.025)))
  
  hundred_year <- rbind(Boot_median[100,2:6], Boot_bottom[100,2:6], Boot_top[100,2:6])
  hundred_year["Source"] <- c("Median", "Bottom", "Top")
  hundred_year <- as.data.frame(hundred_year)
  hundred_year <- gather(hundred_year, key = "Threat_level", value = "Probability", 1:5)
  hundred_year$Threat_level <- recode(hundred_year$Threat_level, "LC" = 1, "NT" = 2, "VU" = 3, "EN" = 4, "CR" = 5)
  hundred_year$Threat_level <- as.factor(hundred_year$Threat_level)
  hundred_year$Probability <- as.numeric(hundred_year$Probability)
  return(hundred_year)
}
