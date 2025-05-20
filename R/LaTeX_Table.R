rm(list = ls())
# setwd("./Thèse/Revision ICML/SYNTHETIC DATA/")
library(tidyverse)
library(MASS)
library(mgcv)
library(mvnfast)
library(dplyr)
library(tidyr)
library(parallel)
library(RColorBrewer)
num_cores <- detectCores() - 1 
print(num_cores)
source("./R/function.R")

# Collect and aggregate the results for the different settings
results_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_length_individuals <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))

i <- 1
for(n in c(12, 12**2, 12**3)){
  for (l in (3:4)){
    print(paste0("l : ",l," --- n : ",n))
    df <- readRDS(paste0("./Config_",i,"/Final_Output.rds"))
    m <- length(unique(df$variable))
    summarised_df <- df %>% group_by(simulation, method, métrique) %>% summarise(Averaged_values = mean(value), .groups = "drop")
    Coverage <- summarised_df %>% filter(métrique == "Coverage")
    Coverage <- Coverage %>% group_by(method) %>%
      summarise(
        Mean_Coverage = 100*mean(Averaged_values, na.rm = TRUE),
        SE_Coverage = 100*1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    Length <- summarised_df %>% filter(métrique == "Length") %>% mutate(Averaged_values = m*Averaged_values)# Je vais *m car j'ai pris la moyenne et non la somme avant
    Length <-  Length %>% group_by(method) %>%
      summarise(
        sqrt_Mean_Length = sqrt(mean(Averaged_values, na.rm = TRUE)),
        sqrt_SE_Length = sqrt(1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())),   
        Mean_Length = mean(Averaged_values, na.rm = TRUE),
        SE_Length = 1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    
    results_length <- results_length %>% cbind(Simu = Length$sqrt_Mean_Length)
    results_coverage <- results_coverage %>% cbind(Simu = Coverage$Mean_Coverage)
    
    uncertainty_length <- uncertainty_length %>% cbind(Simu = Length$sqrt_SE_Length)
    uncertainty_coverage <- uncertainty_coverage %>% cbind(Simu = Coverage$SE_Coverage)
    
    summarised_df_nodes <- df %>% group_by(method, métrique,variable) %>% summarise(Averaged_values = mean(value), .groups = "drop")
    Length_nodes <- summarised_df_nodes %>% filter(métrique == "Length")  %>% group_by(variable) %>%      
      summarise(method_min = method[which.min(Averaged_values)], .groups = 'drop')
    dataframe_methods <- summarised_df_nodes %>% distinct(method) 
    exponent <- log(n,base = 12)
    if(l == 3){
      Length_nodes <- dataframe_methods %>% left_join(Length_nodes, by = character() )%>%
        mutate(nb_occurences = ifelse(method == method_min, 1, 0),
               Groupe = case_when(
                 variable %in% c("V1") ~ "Layer 1", 
                 variable %in% paste0("V",seq(from = 2, to = 1+3**exponent)) ~ "Layer 2", 
                 variable %in% paste0("V",seq(from = 2+3**exponent , to = 2+3**exponent + 1 + 12**exponent))~ "Layer 3",
                 TRUE ~ "Autre")
        )
    }else{
      Length_nodes <- dataframe_methods %>% left_join(Length_nodes, by = character() )%>%
        mutate(nb_occurences = ifelse(method == method_min, 1, 0),
               Groupe = case_when(
                 variable %in% c("V1") ~ "Layer 1", 
                 variable %in% paste0("V",seq(from = 2, to = 1+2**exponent)) ~ "Layer 2", 
                 variable %in% paste0("V",seq(from = 2+2**exponent , to = 2+2**exponent + 1 + 4**exponent))~ "Layer 3",
                 variable %in% paste0("V",seq(from = 2**exponent + 1 + 4**exponent +1 , to =2+2**exponent + 1 + 4**exponent + 1 + 12**exponent))~ "Layer 4",
                 TRUE ~ "Autre")
        )
    }
    Length_nodes <- Length_nodes %>% group_by(method, Groupe)%>%  
      summarise(Occurence = sum(nb_occurences), .groups = 'drop') %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method) 
    Length_nodes <- Length_nodes %>%
      pivot_wider(names_from = method, values_from = Occurence, values_fill = list(Occurence = 0))
    results_length_individuals <- results_length_individuals %>% rbind(Length_nodes)
    i <- i + 1
  }
}
str(results_length_individuals)
config1 <- 100*colSums(results_length_individuals[1:3,-1])/16
config2 <- 100*colSums(results_length_individuals[4:7,-1])/19
config3 <- 100*colSums(results_length_individuals[8:10,-1])/154
config4 <- 100*colSums(results_length_individuals[11:14,-1])/165
config5 <- 100*colSums(results_length_individuals[15:17,-1])/1756
config6 <- 100*colSums(results_length_individuals[18:21,-1])/1801

frequences <- as.data.frame(rbind(config1,config2,config3,config4,config5,config6))
df <- frequences %>%
 tibble::rownames_to_column("configuration") %>%
  pivot_longer(
     cols = -configuration,
    names_to = "method",
    values_to = "freq"
  )
saveRDS(df,"./frequences.rds")
# Table of Coverages
sink("./fichier_LaTeX.txt")
cat(get_table(t(results_length),t(uncertainty_length)))
# Table of Lengths (sqrt)
cat(get_table(t(results_coverage),t(uncertainty_coverage),bold=FALSE))
# Table of Individual Length-minimizer Occurences
cat(get_table_occurence(results_length_individuals[,-1]))
sink()

#####################################################
################## Joint Coverage ###################
#####################################################

# Collect and aggregate the results for the different settings
results_Volume <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_Volume <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
# results_length_individuals <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))

i <- 21
for(n in c(12)){
  for (l in (3:4)){
    print(paste0("l : ",l," --- n : ",n))
    df <- readRDS(paste0("./Config_",i,"/Final_Output_Sigma_joint.rds"))
    summarised_df <- df %>% group_by(simulation, method, métrique) %>% summarise(Averaged_values = mean(V), .groups = "drop")
    Coverage <- summarised_df %>% filter(métrique == "Coverage")
    Coverage <- Coverage %>% group_by(method) %>%
      summarise(
        Mean_Coverage = 100*mean(Averaged_values, na.rm = TRUE),
        SE_Coverage = 100*1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    Volume <- summarised_df %>% filter(métrique == "Volume") 
    Volume <-  Volume %>% group_by(method) %>%
      summarise(
        Mean_Volume = mean(Averaged_values, na.rm = TRUE),
        SE_Volume = 1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    
    results_Volume <- results_Volume %>% cbind(Simu = Volume$Mean_Volume)
    results_coverage <- results_coverage %>% cbind(Simu = Coverage$Mean_Coverage)
    
    uncertainty_Volume <- uncertainty_Volume %>% cbind(Simu = Volume$SE_Volume)
    uncertainty_coverage <- uncertainty_coverage %>% cbind(Simu = Coverage$SE_Coverage)
    
    # summarised_df_nodes <- df %>% group_by(method, métrique,variable) %>% summarise(Averaged_values = mean(value), .groups = "drop")
    # Length_nodes <- summarised_df_nodes %>% filter(métrique == "Length")  %>% group_by(variable) %>%      
    #   summarise(method_min = method[which.min(Averaged_values)], .groups = 'drop')
    # dataframe_methods <- summarised_df_nodes %>% distinct(method) 
    # exponent <- log(n,base = 12)
    # if(l == 3){
    #   Length_nodes <- dataframe_methods %>% left_join(Length_nodes, by = character() )%>%
    #     mutate(nb_occurences = ifelse(method == method_min, 1, 0),
    #            Groupe = case_when(
    #              variable %in% c("V1") ~ "Layer 1", 
    #              variable %in% paste0("V",seq(from = 2, to = 2+3**exponent)) ~ "Layer 2", 
    #              variable %in% paste0("V",seq(from = 2+3**exponent + 1, to = 2+3**exponent + 1 + 12**exponent))~ "Layer 3",
    #              TRUE ~ "Autre")
    #     )
    # }else{
    #   Length_nodes <- dataframe_methods %>% left_join(Length_nodes, by = character() )%>%
    #     mutate(nb_occurences = ifelse(method == method_min, 1, 0),
    #            Groupe = case_when(
    #              variable %in% c("V1") ~ "Layer 1", 
    #              variable %in% paste0("V",seq(from = 2, to = 2+2**exponent)) ~ "Layer 2", 
    #              variable %in% paste0("V",seq(from = 2+2**exponent + 1, to = 2+2**exponent + 1 + 4**exponent))~ "Layer 3",
    #              variable %in% paste0("V",seq(from = 2+2**exponent + 1 + 4**exponent +1 , to =2+2**exponent + 1 + 4**exponent + 1 + 12**exponent))~ "Layer 4",
    #              TRUE ~ "Autre")
    #     )
    # }
    # Length_nodes <- Length_nodes %>% group_by(method, Groupe)%>%  
    #   summarise(Occurence = sum(nb_occurences), .groups = 'drop') %>%
    #   mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
    #   arrange(method) 
    # Length_nodes <- Length_nodes %>%
    #   pivot_wider(names_from = method, values_from = Occurence, values_fill = list(Occurence = 0))
    # results_length_individuals <- results_length_individuals %>% rbind(Length_nodes)
    i <- i + 1
  }
}

# Table of Coverages
sink("./fichier_LaTeX_joint.txt")
cat(get_table(t(results_Volume),t(uncertainty_Volume),digits = 3))
# Table of Lengths (sqrt)
cat(get_table(t(results_coverage),t(uncertainty_coverage),bold=FALSE))
# Table of Individual Length-minimizer Occurences
sink()
