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
library(ggplot2)

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
m <- switch(num_config, 16, 19, 154, 165, 1756, 1801)
# Collect and aggregate the results for the different settings
results_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_length <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
uncertainty_coverage <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))
results_length_individuals <- data.frame(row.names = c("Direct", "OLS","WLS","Combi","MinT"))

data <- readRDS(paste0("./Config_",num_config,"/Final_Output.rds"))
for (j in 1:m){
  df <- data %>% filter(variable == paste0("V",j))
  summarised_df <- df %>% group_by(simulation, method, métrique) %>% summarise(Averaged_values = mean(value), .groups = "drop")
  Coverage <- summarised_df %>% filter(métrique == "Coverage")
  Coverage <- Coverage %>% group_by(method) %>%
    summarise(
      Mean_Coverage = 100*mean(Averaged_values, na.rm = TRUE),
      SE_Coverage = 100*1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
    ) %>%
    mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
    arrange(method)
  Length <- summarised_df %>% filter(métrique == "Length") %>% mutate(Averaged_values = Averaged_values)# Je vais *m car j'ai pris la moyenne et non la somme avant
  Length <-  Length %>% group_by(method) %>%
    summarise(
      sqrt_Mean_Length = sqrt(mean(Averaged_values, na.rm = TRUE)),
      sqrt_SE_Length = sqrt(1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())),   
      Mean_Length = mean(Averaged_values, na.rm = TRUE),
      SE_Length = 1.96*sd(Averaged_values, na.rm = TRUE) / sqrt(n())
    ) %>%
    mutate(method = factor(method, levels = c("Direct", "OLS","WLS","Combi","MinT"))) %>%
    arrange(method)
  mean_result <- data.frame(Mean_Coverage = Coverage$Mean_Coverage,
                            Method = Coverage$method,
                            Mean_Length = Length$Mean_Length,
                            SE_Length = Length$SE_Length,
                            SE_Coverage = Coverage$SE_Coverage)
  alpha = 0.1
  p1 <- ggplot(mean_result, aes(x = Mean_Coverage, y = sqrt(Mean_Length), color = Method, shape = Method)) +
    geom_point(size = 5) + # Points pour les moyennes
    geom_errorbar(aes(ymin = sqrt(Mean_Length - SE_Length), ymax = sqrt(Mean_Length + SE_Length)), width = 0.002, alpha = 1 ) +
    geom_errorbarh(aes(xmin = (Mean_Coverage - SE_Coverage), xmax = (Mean_Coverage + SE_Coverage)), height = (max(sqrt(mean_result$Mean_Length))-min(sqrt(mean_result$Mean_Length)))/15, alpha = 0.8 ) +
    geom_vline(xintercept =  100*(1 - alpha), col = "gray", linetype = "dashed") +
    scale_color_manual(name = "Methods",
                       values = c("Direct" = "gold",
                                  "OLS" = "#FFA500",
                                  "WLS" = "#FF7F7F",
                                  "MinT" = "#9C1111",
                                  "Combi" = "red")) +
    scale_shape_manual(
      name = "Methods",
      values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "Combi" = 7)) +
    scale_x_continuous(limits = c(89.985, 90.015),
                       breaks = seq(89.99, 90.01, by = 0.01),
                       labels = function(x) sprintf("%.2f", x)) +
    labs(title = "", x = "Coverage", y = "Length") +
    theme_minimal() +
    theme(axis.text = element_text(size = 26),
          axis.text.y = element_text(size = 26, margin = margin(r =  -20)),
          axis.title.y = element_text(size = 30, margin = margin(r =  20)),
          axis.title.x =  element_text(size = 29),
          plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0.5)) +
    theme(legend.position = "none")
  
  ggsave(filename = paste0("./Plot/Config_",num_config,"_CovLength_",j,".pdf"), plot = p1, width = 7, height = 6, units = "in")
}