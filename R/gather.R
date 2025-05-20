rm(list = ls())
library(tidyverse)
library(MASS)
library(mgcv)
library(mvnfast)
library(dplyr)
library(tidyr)
library(stats)
library(float)
library(RColorBrewer)
library(parallel)

source("./R/function.R")

num_cores <- detectCores() - 1 
print(num_cores)

args <- commandArgs(trailingOnly = TRUE)
for (i in seq_along(args)) {
  if (args[i] == "--simulation_indice" && i < length(args)) {
    simulation_indice <- as.numeric(args[i + 1])
  }
  if (args[i] == "--config" && i < length(args)) {
    config_vals <- strsplit(args[i + 1], " ")[[1]]
    l <- as.numeric(config_vals[1])
    n <- as.numeric(config_vals[2])
    n_obs <- as.numeric(config_vals[3])
    num_simulations <- as.numeric(config_vals[4])
    num_config <- as.numeric(config_vals[5])
  }
}
print(paste0("Gather for config ",num_config))
df_list <- list()
df_list_joint <- list()
for (i in 1:num_simulations) {
  simulation_results <- readRDS(paste0("./Config_",num_config,"/simulation_",i-1,".rds"))
  df_list[[i]] <- simulation_results
  simulation_results_joint <- readRDS(paste0("./Config_",num_config,"/simulation_joint_",i-1,".rds"))
  df_list_joint[[i]] <- simulation_results_joint
  
}
results_proba <- bind_rows(df_list)
saveRDS(results_proba,paste0("./Config_",num_config,"/Final_Output.rds"))
results_proba_joint <- bind_rows(df_list_joint)
saveRDS(results_proba_joint,paste0("./Config_",num_config,"/Final_Output_joint.rds"))
