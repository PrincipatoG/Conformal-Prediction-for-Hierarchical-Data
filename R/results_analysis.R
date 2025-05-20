rm(list = ls())

library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(lubridate)
library(mgcv)
library(ggplot2)
library(stringr)

args=commandArgs(TRUE)
n <- as.numeric(args[1])
n_obs <- as.numeric(args[2])
num_simulations <- as.numeric(args[3])
config <- args[4]

# model <- "gam"
# n <- 5
# n_obs <- 100000
# num_simulations <- 1000
# config <- "Simu_easy_many"
if(str_detect(config, "hard") ){
  difficulty <- "hard"
}else{
  difficulty <- "easy"
}
if (config == "Simu_hard_many"){
  scaler_x <- 1
  scaler_y <- 2.5
}
if (config == "Simu_easy_many"){
  scaler_x <- 1
  scaler_y <- 5
}
if (config == "Simu_hard_few"){
  scaler_x <- 1
  scaler_y <- 2.5
}
if (config == "Simu_easy_few"){
  scaler_x <- 1 
  scaler_y <- 5
}
# n <- 5
# m <- n + 3
# n_obs <- 1000
# num_simulations <- 100
# model <- "gam"
methods <- c("Direct", "OLS","WLS","MinT","Combi")
alpha <- 0.1
result <- readRDS(file = paste0("./Outputs/results_",model,"_",n,"_",n_obs,"_",num_simulations,".rds")) #%>% mutate(Method = factor(method, levels = c("Direct","OLS","WLS","Combi","MinT")))
# result <- readRDS(file = paste0("./Outputs/results_new_",difficulty,"_",model,"_",n,"_",n_obs,"_",num_simulations,".rds")) #%>% mutate(Method = factor(method, levels = c("Direct","OLS","WLS","Combi","MinT")))
# result_nodes <- readRDS(file = "./Output/Simulation_1/Result_nodes.rds")
summarised_results <- result %>% group_by(simulation, method, métrique) %>% summarise(Averaged_values = mean(value), .groups = "drop")
Coverage <- summarised_results %>% filter(métrique == "Coverage")
Length <- summarised_results %>% filter(métrique == "Length") %>% mutate(Averaged_values = m*Averaged_values)# Je vais *m car j'ai pris la moyenne et non la somme avant

result_global <- data.frame(Method = Length$method, Length = Length$Averaged_values, Coverage = Coverage$Averaged_values, Target = rep(0.9, length(methods)) )

# saveRDS(result_global,"../../Plot Bachir le boss/SYNTETHIC DATA/result_global_artificial.rds")
# saveRDS(result,"../../Plot Bachir le boss/SYNTETHIC DATA/result_nodes_artificial.rds")

mean_result <- result_global %>%
  group_by(Method) %>%
  summarise(
    Mean_Coverage = mean(Coverage, na.rm = TRUE),
    Mean_Length = mean(Length, na.rm = TRUE),
    sqrt_Mean_Length = sqrt(mean(Length, na.rm = TRUE)),
    SE_Coverage = 1.96*sd(Coverage, na.rm = TRUE) / sqrt(n()),
    SE_Length = 1.96*sd(Length, na.rm = TRUE) / sqrt(n())
  )
mean_result <- mean_result %>% mutate(Method = factor(Method, levels = c("Direct","OLS","WLS","Combi","MinT")))
# C'est une modification de mean_result que je veux pour le tableau LaTeX

df <- data.frame(
  A = c(2.3, 4.5, 1.2),
  B = c(1.1, 3.2, 1.2),
  C = c(3.4, 2.1, 0.9)
)
uncertainty <- df <- data.frame(
  A = c(0.3, 0.5, 0.2),
  B = c(0.1, 0.2, 0.2),
  C = c(0.4, 0.1, 0.3)
)

df <- t(mean_result[c("Method","Mean_Length")])
colnames(df) <- as.character(unlist(df[1, ]))
df <- df[-1, ]  
uncertainty <- t(mean_result[c("Method","SE_Length")])
uncertainty <- uncertainty[-1, ]  

get_table <- function(data, uncertainty) {
  latex_rows <- sapply(seq_len(nrow(data)),  function(i) {
    row <- data[i, ]
    row_unc <- uncertainty[i, ]
    min_val <- min(row, na.rm = TRUE)
    
    formatted_row <- mapply(function(x, u) {
      value_str <- paste0(format(x, digits = 3), " $\\pm$ ", format(u, digits = 2))
      if (x == min_val) {
        return(paste0("\\textbf{", value_str, "}"))
      } else {
        return(value_str)
      }
    }, row, row_unc)
    
    paste(formatted_row, collapse = " & ")
  })
  
  # En-têtes
  col_names <- paste(colnames(data), collapse = " & ")
  
  # Corps
  body <- paste(latex_rows, collapse = " \\\\\n")
  
  # Tableau complet
  table_latex <- paste0(
    "\\begin{tabular}{", paste(rep("c", ncol(data)), collapse = ""), "}\n",
    col_names, " \\\\\n\\hline\n",
    body, " \\\\\n",
    "\\end{tabular}"
  )
  
  return(table_latex)
}
cat(get_table(df,uncertainty))




p1 <- ggplot(mean_result, aes(x = 100 * Mean_Coverage, y = Mean_Length, color = Method, shape = Method)) +
  geom_point(size = 7) + # Points pour les moyennes
  geom_errorbar(aes(ymin = Mean_Length - SE_Length, ymax = Mean_Length + SE_Length), width = 0.02 * scaler_x) +
  geom_errorbarh(aes(xmin = 100 * (Mean_Coverage - SE_Coverage), xmax = 100 * (Mean_Coverage + SE_Coverage)), height = 0.2 * scaler_y) +
  geom_vline(xintercept = 100 * (1 - alpha), col = "gray", linetype = "dashed") +
  # geom_text(aes(label = pval_text), vjust = -1.5, size = 5, color = "black") + # Ajouter les p-values
  scale_color_manual(name = "Methods",
                     values = c("Direct" = "gold",
                                "OLS" = "#FFA500",
                                "WLS" = "#FF7F7F",
                                "MinT" = "#9C1111",
                                "Combi" = "red")) +
  scale_shape_manual(
    name = "Methods",
    values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "Combi" = 7)) +
  scale_x_continuous(limits = c(89.5, 90.5),
                     breaks = seq(89.5, 90.5, by = 0.1),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))
p1
# Sauvegarder le graphique
ggsave(paste0("./Plot_config/",config,"/Global_Coverage_Length_Averages_",model,"_",n,"_",n_obs,"_",num_simulations,".png"), p1)
print(p1)
b1 <- ggplot(mean_result, aes(x = Method, y = Mean_Length, color = Method, shape = Method)) +
  geom_point(size = 6) + 
  labs(title = "",
       x = "",
       y = "Length") + 
  theme_minimal() + 
  geom_errorbar(aes(ymin = Mean_Length - SE_Length, ymax = Mean_Length + SE_Length), width = 0.25) +
  scale_shape_manual(
    name = "Methods",
    values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "MinT_oracle" = 16, "Combi" = 7),
    labels = c("Direct" = "Direct",
               "OLS" = "OLS",
               "WLS" = "WLS",
               "MinT" = "MinT",
               "MinT_oracle" = "oracle MinT",
               "Combi" = "Combi")) +
  scale_color_manual(name = "Methods",
                     values = c("Direct" = "gold",  # gold
                                "OLS" = "#FFA500",   # Orange
                                "WLS" = "#FF7F7F", # Rouge clair ajusté
                                "Combi" = "red",       # Rouge pour l'étoile
                                "MinT" = "#9C1111",
                                "MinT_oracle" = "#482076"),
                     labels = c("Direct" = "Direct",
                                "OLS" = "OLS",
                                "WLS" = "WLS",
                                "MinT" = "MinT",
                                "MinT_oracle" = "oracle MinT",
                                "Combi" = "Combi"))  +  # Rouge foncé  
  theme_void() +
  theme(axis.line.x = element_blank(),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.line.y = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        # axis.text.x = element_text(size = 11),
        axis.text.y = element_text(color = "black", size = 12, margin = margin(r = 5)),
        axis.title.y = element_text(size = 13, margin = margin(r = 15), angle = 90))+
  guides(fill = guide_legend(override.aes = list(size = c(2, 2, 2, 2, 2, 2)))) # Réduire la taille de la moyenne dans la légende
b1
ggsave(paste0("./Plot_config_final/",config,"/boxplot_artificial.pdf"), b1)

# nodes <- unique(result$variable)
# # node <- nodes[[1]]
# for (node in nodes){
#   
#   if(node %in% c("V2","V4")){
#     scaler_y <- scaler_y/2
#     print(scaler_y)
#   }
#   data_node <- result %>% filter(variable == node)
#   Coverage_node <- data_node %>% filter(métrique == "Coverage")
#   Length_node <- data_node %>% filter(métrique == "Length") %>% mutate(value = sqrt(value))
#   data_node <- data.frame(Method = Length_node$method, Length = Length_node$value, Coverage = Coverage_node$value, Target = rep(0.9, length(methods)) )
#   data_node <- data_node %>% mutate(Method = factor(Method, levels = c("Direct","OLS","WLS", "Combi","MinT")))
#   
#   mean_result_node <- data_node %>%
#     group_by(Method) %>%
#     summarise(
#       Mean_Coverage = mean(Coverage, na.rm = TRUE),
#       Mean_Length = mean(Length, na.rm = TRUE),
#       SE_Coverage = 1.96*sd(Coverage, na.rm = TRUE) / sqrt(n()),
#       SE_Length = 1.96*sd(Length, na.rm = TRUE) / sqrt(n())
#     )
#   
#   p_node <- ggplot(mean_result_node, aes(x = 100 * Mean_Coverage, y = Mean_Length, color = Method, shape = Method)) +
#     geom_point(size = 7) + # Points pour les moyennes
#     geom_errorbar(aes(ymin = Mean_Length - SE_Length, ymax = Mean_Length + SE_Length), width = 0.02 * scaler_x) +
#     geom_errorbarh(aes(xmin = 100 * (Mean_Coverage - SE_Coverage), xmax = 100 * (Mean_Coverage + SE_Coverage)), height = 0.2 * scaler_y) +
#     geom_vline(xintercept = 100 * (1 - alpha), col = "gray", linetype = "dashed") +
#     # geom_text(aes(label = pval_text), vjust = -1.5, size = 5, color = "black") + # Ajouter les p-values
#     scale_color_manual(name = "Methods",
#                        values = c("Direct" = "gold",
#                                   "OLS" = "#FFA500",
#                                   "WLS" = "#FF7F7F",
#                                   "MinT" = "#9C1111",
#                                   "Combi" = "red")) +
#     scale_shape_manual(
#       name = "Methods",
#       values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "Combi" = 7)) +
#     scale_x_continuous(limits = c(89.5, 90.5),
#                        breaks = seq(89.5, 90.5, by = 0.1),
#                        labels = function(x) sprintf("%.1f", x)) +
#     labs(title = "", x = "", y = "") +
#     theme_minimal() +
#     theme(legend.position = c(0.9, 0.9),
#           legend.justification = c(1, 1),
#           legend.title = element_text(size = 18),
#           legend.text = element_text(size = 14))
#   p_node
#   ggsave(paste0("./Plot_config/",config,"/Nodes_Coverage_Length_",node,"_",model,"_",n,"_",n_obs,"_",num_simulations,".png"), p_node)
#   # print(p_node)
# }

node <- "V1"
if(node %in% c("V2","V4")){
  scaler_y <- scaler_y/2
  print(scaler_y)
}
data_node <- result %>% filter(variable == node)
Coverage_node <- data_node %>% filter(métrique == "Coverage")
Length_node <- data_node %>% filter(métrique == "Length") %>% mutate(value = value)
data_node <- data.frame(Method = Length_node$method, Length = Length_node$value, Coverage = Coverage_node$value, Target = rep(0.9, length(methods)) )
data_node <- data_node %>% mutate(Method = factor(Method, levels = c("Direct","OLS","WLS", "Combi","MinT")))

mean_result_node <- data_node %>%
  group_by(Method) %>%
  summarise(
    Mean_Coverage = mean(Coverage, na.rm = TRUE),
    Mean_Length = mean(Length, na.rm = TRUE),
    sqrt_Mean_Length = sqrt(mean(Length, na.rm = TRUE)),
    SE_Coverage = 1.96*sd(Coverage, na.rm = TRUE) / sqrt(n()),
    SE_Length = 1.96*sd(Length, na.rm = TRUE) / sqrt(n())
  )

p_node <- ggplot(mean_result_node, aes(x = 100 * Mean_Coverage, y = Mean_Length, color = Method, shape = Method)) +
  geom_point(size = 5) + 
  geom_errorbar(aes(ymin = Mean_Length - SE_Length, ymax = Mean_Length + SE_Length), width = 0.02 * scaler_x) +
  geom_errorbarh(aes(xmin = 100 * (Mean_Coverage - SE_Coverage), xmax = 100 * (Mean_Coverage + SE_Coverage)), height = 0.2 * scaler_y) +
  geom_vline(xintercept = 100 * (1 - alpha), col = "gray", linetype = "dashed") +
  scale_color_manual(name = "Methods",
                     values = c("Direct" = "gold",
                                "OLS" = "#FFA500",
                                "WLS" = "#FF7F7F",
                                "MinT" = "#9C1111",
                                "Combi" = "red")) +
  scale_shape_manual(
    name = "Methods",
    values = c("Direct" = 17, "OLS" = 15, "WLS" = 18, "MinT" = 16, "Combi" = 7)) +
  scale_x_continuous(limits = c(89.5, 90.5),
                     breaks = seq(89.5, 90.5, by = 0.1),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme_void()+
  theme(legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.line.y = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        # axis.text.x = element_text(size = 11),
        axis.text.y = element_text(color = "black", size = 12, margin = margin(r = 5)),
        axis.title.y = element_text(size = 13, margin = margin(r = 15), angle = 90))
p_node
ggsave(paste0("./Plot_config/",config,"/Nodes_Coverage_Length_",node,"_",model,"_",n,"_",n_obs,"_",num_simulations,".png"), p_node)
