rm(list = ls())
library(tidyverse)
library(MASS)
library(mgcv)
library(mvnfast)
library(dplyr)
library(tidyr)
library(stats)
library(parallel)
library(float)
library(RColorBrewer)
# num_cores <- detectCores() - 1 
# print(num_cores)
# setwd("./Thèse/Revision ICML/SYNTHETIC DATA/")
source("./R/function.R")
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
# simulation_indice <- 1
# l <- 3
# n <- 12
# n_obs <- 1000000
# num_simulations <- 1000
# num_config <- 21
config <- "easy"
print(paste0("Running for Config ",num_config, ", l = ",l,", n = ",n,", n_obs = ", n_obs,", num_simulations = ", num_simulations, ", simulation indice = ",simulation_indice))
# Get the Structural Matrix for the given configuration
S <- generate_S(n = n, l = l, balanced=TRUE)
# S <- rbind(c(1,1),diag(c(1,1)))
m <- nrow(S)
alpha <- 0.1
model <- "gam"

print(paste0("Starting Simulation ",simulation_indice, " for l = ",l, " and n = ",n))
list_A <- c("Sigma","Diag","Id")
for (A in list_A){
  results_simu <- simulation_multi_plot(S, n_obs, model, A)
  sim_data_proba <- results_simu[[1]]
  sim_data_proba_joint <- do.call(rbind, results_simu[[2]]) 
  sim_data_proba_joint$method <- sub("\\..*$", "", rownames(sim_data_proba_joint))
  rownames(sim_data_proba_joint) <- NULL
  
  sim_data_proba$simulation <- factor(rep(simulation_indice, m*3*5)) 
  sim_data_proba_joint$simulation <- factor(rep(simulation_indice, 3*5)) 
  
  print(paste0("Ending Simulation ",simulation_indice, " for l = ",l, " and n = ",n))
  saveRDS(sim_data_proba, file = paste0("./Config_",num_config,"/A_",A,"simulation_",simulation_indice,".rds"))
  saveRDS(sim_data_proba_joint, file = paste0("./Config_",num_config,"/A_",A,"simulation_joint_",simulation_indice,".rds"))
}

