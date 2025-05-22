rm(list = ls())
library(MASS)
library(mgcv)
library(mvnfast)
library(dplyr)
library(tidyr)
library(stats)
library(float)
library(RColorBrewer)

source("./R/function.R")

# Collect and aggregate the results for the different settings
results_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))

i <- 1
for(n in c(12, 12**2, 12**3)){
  for (l in (3:4)){
    print(paste0("l : ",l," --- n : ",n))
    df <- readRDS(paste0("./Config_",i,"/Final_Output.rds"))
    m <- length(unique(df$variable))
    summarised_df <- df %>% group_by(simulation, method, metric) %>% summarise(Averaged_values = mean(value), .groups = "drop")
    Coverage <- summarised_df %>% filter(metric == "Coverage")
    Coverage <- Coverage %>% group_by(method) %>%
      summarise(
        Mean_Coverage = 100*mean(Averaged_values, na.rm = TRUE),
        SE_Coverage = 100*1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    Length <- summarised_df %>% filter(metric == "Length") %>% mutate(Averaged_values = m*Averaged_values) # the "times m" is to get the sum from the average
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
  
    i <- i + 1
  }
}

# Table of Coverages
sink("./LaTeX_file.txt")
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

i <- 1
for(n in c(12)){
  for (l in (3:4)){
    print(paste0("l : ",l," --- n : ",n))
    df <- readRDS(paste0("./Config_",i,"/Final_Output_joint.rds"))
    summarised_df <- df %>% group_by(simulation, method, metric) %>% summarise(Averaged_values = mean(V), .groups = "drop")
    Coverage <- summarised_df %>% filter(metric == "Coverage")
    Coverage <- Coverage %>% group_by(method) %>%
      summarise(
        Mean_Coverage = 100*mean(Averaged_values, na.rm = TRUE),
        SE_Coverage = 100*1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
      arrange(method)
    Volume <- summarised_df %>% filter(metric == "Volume") 
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
    
    i <- i + 1
  }
}

# Table of Coverages
sink("./LaTeX_joint_file.txt")
cat(get_table(t(results_Volume),t(uncertainty_Volume),digits = 3))
# Table of Lengths (sqrt)
cat(get_table(t(results_coverage),t(uncertainty_coverage),bold=FALSE))
# Table of Individual Length-minimizer Occurences
sink()
