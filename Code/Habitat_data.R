require(dplyr)
require(stringr)
require(ggplot2)
require(ggpubr)
require(rstatix)


## read in habitat data
## assign habitat specialist vs generalist to all species


Process_habitat_data <- function(All_Habitat){
  Habitat_data <- All_Habitat[,c("code", "Species")]
  Habitat_data <- Habitat_data[which(!is.na(Habitat_data$code)),]
  Habitat_data <- Habitat_data %>% mutate(code = str_extract(code, "([0-9]+)(?=\\.?)"))
  Habitat_data <- unique(Habitat_data)
  Habitat_data$type <- NA
  Habitat_data <- Habitat_data %>% group_by(Species) %>% mutate(type = if(n()>1) "Generalist" else "Specialist") %>% ungroup
  Habitat_data <- unique(Habitat_data[,c("Species", "type")])
  return(Habitat_data)
}


Plot_habitat_100 <- function(hundred_year, ylabel = "Probability of extinction at 100 years\n", xlabel = "\nIUCN Red List Category", leg_pos = c(0.2,0.8), y_limits = ylim(0,0.2), plottag = ""){
  p <- ggplot(data = subset(hundred_year, Source %in% c("Median")), aes(x = Threat_level, y = Probability, fill = Type)) + scale_fill_manual(values = c("darkcyan", "darkorange", "darkred"), name = "")
  p <- p + geom_bar(stat = "identity", position = "dodge") + scale_x_discrete(breaks = 1:5, labels=c("LC","NT","VU", "EN","CR")) + y_limits
  p <- p + labs(y = ylabel, x = xlabel, tag = plottag)
  p <- p + geom_errorbar(aes(ymin= hundred_year$Probability[hundred_year$Source == "Bottom"], ymax=hundred_year$Probability[hundred_year$Source == "Top"]), width=.2, position=position_dodge(.9)) 
  p <- p + theme(panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
                 axis.text.y = element_text(size=16), axis.title = element_text(size=20), axis.text.x = element_text(size=16), legend.position = leg_pos, legend.text = element_text(size=20), 
                 legend.title = element_text(size=24), strip.text = element_text(size=14), plot.tag.position = "topright")
  #p <- p + stat_pvalue_manual(data = p_values, x = "Threat_level")
  return(p)
}
