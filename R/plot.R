rm(list = ls())
setwd("./Thèse/Revision ICML/SYNTHETIC DATA/")
library(tidyverse)
library(MASS)
library(mgcv)
library(mvnfast)
library(dplyr)
library(tidyr)
library(parallel)
library(RColorBrewer)
library(ggplot2)
library(cowplot) 

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
n <- 12
l <- 3
data <- readRDS(paste0("./Config_",i,"/Final_Output.rds"))
for (j in 1:19){
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
  p1
  ggsave(filename = paste0("./Plot/Config_",i,"_CovLength_",j,".pdf"), width = 7, height = 6, units = "in")
}

data_example <- data.frame(
  Mean_Coverage = c(90.003,90.002),
  Mean_Length = c(87^2,116^2),  # car dans ton graphique y = sqrt(Mean_Length)
  SE_Coverage = c(0.01,0.008),
  SE_Length = c(40^2,35^2),     # erreur verticale au niveau de sqrt(Mean_Length)
  Method = c("dot","cross")
)

# Création du graphique
p_example <- ggplot(data_example, aes(x = Mean_Coverage, y = sqrt(Mean_Length),shape = Method, color = Method)) +
  geom_errorbar(aes(ymin = sqrt(Mean_Length - SE_Length), ymax = sqrt(Mean_Length + SE_Length)), 
                width = 0.002, , color = "black", alpha = 1) +
  geom_errorbarh(aes(xmin = Mean_Coverage - SE_Coverage, xmax = Mean_Coverage + SE_Coverage), 
                 height = 6, color = "black", alpha = 1) +
  geom_point(size = 15) +
  geom_vline(xintercept =  90, col = "gray", linetype = "dashed") +
  scale_x_continuous(limits = c(89.985, 90.015),
                     breaks = seq(89.99, 90.01, by = 0.01),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(limits = c(70, 130)) +
  scale_shape_manual(name = "Method",
                     values = c("dot"=16,"cross"=15))+
  scale_color_manual(name = "Method",
                     values = c("dot"="blue","cross"="black"))+
  labs(x = "Coverage", y = "Length") +
  theme_minimal() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
        legend.position = "none")

print(p_example)
ggsave(filename = paste0("./Plot/Example.pdf"), p_example, width = 7, height = 6, units = "in")

p_legend <- ggplot(mean_result, aes(x = Mean_Coverage, y = sqrt(Mean_Length), color = Method, shape = Method)) +
  geom_point(size = 5) + # Points pour les moyennes
  geom_errorbar(aes(ymin = sqrt(Mean_Length - SE_Length), ymax = sqrt(Mean_Length + SE_Length)), width = 0.005, alpha = 0.8 ) +
  geom_errorbarh(aes(xmin = (Mean_Coverage - SE_Coverage), xmax = (Mean_Coverage + SE_Coverage)), height = (max(sqrt(mean_result$Mean_Length))-min(sqrt(mean_result$Mean_Length)))/5, alpha = 0.8 ) +
  geom_vline(xintercept =  100*(1 - alpha), col = "gray", linetype = "dashed") +
  scale_color_manual(name = "Methods",
                     values = c("Direct" = "gold",
                                "OLS" = "#FFA500",
                                "WLS" = "#FF7F7F",
                                "MinT" = "#9C1111",
                                "Combi" = "red"),
                     guide ="legend") +
  scale_shape_manual(
    name = "Methods",
    values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "Combi" = 7),
    guide ="legend") +
  scale_x_continuous(limits = c(89.95, 90.05),
                     breaks = seq(89.95, 90.05, by = 0.025),
                     labels = function(x) sprintf("%.3f", x)) +
  labs(title = "", x = "Coverage", y = "Length") +
  guides(color = guide_legend(title = "Methods"), shape = guide_legend(title = "Methods")) +  # Fusionne visuellement
  theme_minimal() +
  theme(axis.text = element_text(size = 15),
        axis.title =  element_text(size = 20),
        legend.background =  element_rect(linetype= "solid"))+
  theme(legend.position = "bottom")
p_legend
ggsave(filename = paste0("./Plot/legend.pdf"), plot = p_legend,width = 4, height = 4)


# Les données pour la légende seulement
legend_data <- data.frame(
  Method = c("Direct", "OLS", "WLS", "MinT", "Combi"),
  x = 1:5,
  y = rep(1, 5)  # tous alignés
)

# Créer une "légende" personnalisée avec un graphique bidon
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Method, shape = Method)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "Direct" = "gold",
    "OLS" = "#FFA500",
    "WLS" = "#FF7F7F",
    "MinT" = "#9C1111",
    "Combi" = "red"
  )) +
  scale_shape_manual(values = c(
    "Direct" = 17,
    "OLS" = 15,
    "WLS" = 18,
    "MinT" = 16,
    "Combi" = 7
  )) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.background =  element_rect(linetype= "solid")
  ) +
  guides(
    color = guide_legend("Methods", override.aes = list(shape = c(17, 15, 18, 16, 7))),
    shape = "none"
  )
legend_plot
legend_bottom <- cowplot::get_legend(legend_plot)

pdf("./Plot/legend.pdf", width = 7, height = 1.5)
grid::grid.newpage()
grid::grid.draw(legend_bottom)
dev.off()
