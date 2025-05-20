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
noise_type <- "A"
alpha <- 0.1
model <- "gam"
print(paste0("Running for Config ",num_config, ", l = ",l,", n = ",n,", n_obs = ", n_obs,", num_simulations = ", num_simulations, ", simulation indice = ",simulation_indice))

# Get the Structural Matrix for the given configuration
S <- generate_S(n = n, l = l, balanced=TRUE)
m <- nrow(S)

print(paste0("Starting Simulation ",simulation_indice, " for l = ",l, " and n = ",n))

# Run an instance of the experiment
results_simu <- simulation_multi(S, n_obs, model)

# Extract the component-wise results
sim_data_proba <- results_simu[[1]]

# Extract the joint coverage results
sim_data_proba_joint <- do.call(rbind, results_simu[[2]]) 
sim_data_proba_joint$method <- sub("\\..*$", "", rownames(sim_data_proba_joint))
rownames(sim_data_proba_joint) <- NULL

print(paste0("Ending Simulation ",simulation_indice, " for l = ",l, " and n = ",n))

# Save the results
sim_data_proba$simulation <- factor(rep(simulation_indice, m*3*5)) 
sim_data_proba_joint$simulation <- factor(rep(simulation_indice, 3*5)) 
saveRDS(sim_data_proba, file = paste0("./Config_",num_config,"/simulation_",simulation_indice,".rds"))
saveRDS(sim_data_proba_joint, file = paste0("./Config_",num_config,"/simulation_joint_",simulation_indice,".rds"))
