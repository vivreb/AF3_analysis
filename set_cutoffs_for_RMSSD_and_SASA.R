library(readr)
library(dplyr)
library(tidyr)
library(mixtools)
library(AdaptGauss)
library(ggplot2)


all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff <- read_csv("all_kinases_pae_weighted_rmssd_and_surface_stu.csv") %>% 
  dplyr::select(-c(...1))


sorted_rmssd_vector <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>% 
  arrange(pae_weighted_rmssd) %>% 
  pull(pae_weighted_rmssd)

set.seed(124)
fit_RMSSD <- normalmixEM(sorted_rmssd_vector, mean.constr = c(0.9, 8),
                         sd.constr = c(0.5, NA))

Chi2testMixtures(Data = sorted_rmssd_vector,
                 Means = fit_RMSSD$mu,
                 SDs = fit_RMSSD$sigma,
                 Weights = fit_RMSSD$lambda,
                 PlotIt = TRUE)

plot(fit_RMSSD, whichplots = 2, marginal = T, alpha = 0.5)

plotFDR(fit_RMSSD$posterior, alpha = 0.05, compH0 = 1, pctfdr = 0.05)


sorted_SASA_vector <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>% 
  mutate(absolute_SASA_change = abs(change_in_residue_surface)) %>% 
  arrange(absolute_SASA_change) %>% 
  pull(absolute_SASA_change)

set.seed(124)
fit_SASA <- normalmixEM(sorted_SASA_vector, mean.constr = c(15, 150),
                        sd.constr = c(40, NA))

Chi2testMixtures(Data = sorted_SASA_vector,
                 Means = fit_SASA$mu,
                 SDs = fit_SASA$sigma,
                 Weights = fit_SASA$lambda,
                 PlotIt = TRUE)

plot(fit_SASA, whichplots = 2)

plotFDR(fit_SASA$posterior, alpha = 0.1, compH0 = 2)


# Modified function to extract FDR values without plotting
extractFDRValues <- function(post1, post2 = NULL, compH0 = 1, pctfdr = 0.3) {
  # Number of data points
  n <- dim(post1)[1]
  
  # Compute cumulative sums and FDR for the first posterior probabilities (post1)
  cs1 <- cumsum(post1[, compH0])
  fdr1 <- cs1 / (1:n)
  
  # Initialize a list to store the results
  result <- list(fdr1 = fdr1)
  
  # If there is a second posterior probability set (post2)
  if (!is.null(post2)) {
    cs2 <- cumsum(post2[, compH0])
    fdr2 <- cs2 / (1:n)
    result$fdr2 <- fdr2
  }
  
  # Return the FDR values
  return(result)
}

fdr_results_RMSSD <- extractFDRValues(fit_RMSSD$posterior, compH0 = 1)
min(fdr_results_RMSSD[["fdr1"]])

index_for_RMSSD_cutoff <- which(fdr_results_RMSSD[["fdr1"]] == min(fdr_results_RMSSD[["fdr1"]]))

sorted_rmssd_vector[index_for_RMSSD_cutoff]

sorted_rmssd_vector[241] # 50% 
sorted_rmssd_vector[247] # 95%


fdr_results_SASA <- extractFDRValues(fit_SASA$posterior, compH0 = 2)
min(fdr_results_SASA[["fdr1"]])

index_for_SASA_cutoff <- which(fdr_results_SASA[["fdr1"]] == min(fdr_results_SASA[["fdr1"]]))
sorted_SASA_vector[index_for_SASA_cutoff]

sorted_SASA_vector[285] # 50%
sorted_SASA_vector[313] # 95%



#################### Nice plot #################### 

set.seed(124)

density_RMSSD_background <- rnorm(10000, mean = fit_RMSSD$mu[1], sd = fit_RMSSD$sigma[1]) #* fit_RMSSD$lambda[1]
density_RMSSD_changing <- rnorm(10000, mean = fit_RMSSD$mu[2], sd = fit_RMSSD$sigma[2]) #* fit_RMSSD$lambda[2]

densities <- c(density_RMSSD_background, density_RMSSD_changing) 
densities <- data.frame(densities,id = rep(c(1, 2), each=10000), 
                        weights = rep(c(fit_RMSSD$lambda[1], fit_RMSSD$lambda[2]), each=10000))

scale_factor = 1370

rmssd_plot <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>%  
  ggplot(aes(x = pae_weighted_rmssd)) +
  geom_histogram(binwidth = 0.2) +
  stat_density(data = densities %>% filter(id == 1), aes(x = densities, y = ..density.. * scale_factor * fit_RMSSD$lambda[1] / sum(..density..)),  color = "#CBBBBB", position = "identity", geom = "line", linewidth = 0.5) +
  stat_density(data = densities %>% filter(id == 2), aes(x = densities, y = ..density.. * scale_factor * fit_RMSSD$lambda[2] / sum(..density..)),  color = "#95ACD3", position = "identity", geom = "line", linewidth = 0.5) +
  scale_x_continuous(limits = c(-1, 55)) +
  geom_vline(xintercept = sorted_rmssd_vector[237], linetype = "dashed") +
  #geom_vline(xintercept = 0.8, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "RMSSD (Å)", y = "Count")
rmssd_plot

pdf("Z:/Viviane/Computational/AF3/AF3_pae_weighted_RMSSD_histogram_with_densities_new_cutoff_.pdf",
    width = 5,
    height = 2,
    pointsize = 15)
rmssd_plot
dev.off()

set.seed(124)

density_SASA_background <- rnorm(1000, mean = fit_SASA$mu[1], sd = fit_SASA$sigma[1]) #* fit_SASA$lambda[1]
density_SASA_changing <- rnorm(1000, mean = fit_SASA$mu[2], sd = fit_SASA$sigma[2]) #* fit_SASA$lambda[2]

densities <- c(density_SASA_background, density_SASA_changing) 
densities <- data.frame(densities,id = rep(c(1, 2), each=1000), 
                        weights = rep(c(fit_SASA$lambda[1], fit_SASA$lambda[2]), each=1000))

scale_factor = 1450

sasa_plot <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>%  
  ggplot(aes(x = abs(change_in_residue_surface))) +
  geom_histogram(binwidth = 3) +
  stat_density(data = densities %>% filter(id == 1), aes(x = densities, y = ..density.. * scale_factor * fit_SASA$lambda[1] / sum(..density..)),  color = "#CBBBBB", position = "identity", geom = "line", linewidth = 0.5) +
  stat_density(data = densities %>% filter(id == 2), aes(x = densities, y = ..density.. * scale_factor * fit_SASA$lambda[2] / sum(..density..)),  color = "#95ACD3", position = "identity", geom = "line", linewidth = 0.5) +
  scale_x_continuous(limits = c(-1, 750)) +
  geom_vline(xintercept = sorted_SASA_vector[280], linetype = "dashed") +
  #geom_vline(xintercept = 15, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = expression("Absolute change in SASA at ATP-binding site (Å"^2*")"), y = "Count")
sasa_plot

pdf("Z:/Viviane/Computational/AF3/AF3_sasa_histogram_with_densities_new_cutoff_.pdf",
    width = 5,
    height = 2,
    pointsize = 15)
sasa_plot
dev.off()


all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>% 
  mutate(changing = ifelse(pae_weighted_rmssd > sorted_rmssd_vector[237], "Yes", "No")) %>% 
  mutate(changing = ifelse(abs(change_in_residue_surface) > sorted_SASA_vector[280], "Yes", changing))


all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>%  write.csv("all_kinases_pae_weighted_rmssd_and_surface_stu_gmm_based_cutoff.csv")


test <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>% 
  mutate(changing_rmssd = ifelse(pae_weighted_rmssd > sorted_rmssd_vector[237], "yes", "no")) %>% 
  mutate(changing_sasa = ifelse(abs(change_in_residue_surface) > sorted_SASA_vector[280], "yes", "no")) %>% 
  mutate(change_in_both = ifelse(changing_rmssd == "yes" & changing_sasa == "yes", "yes", "no")) %>% 
  mutate(any_change = ifelse(changing_rmssd == "yes" | changing_sasa == "yes", "yes", "no")) %>% 
  mutate(no_change = ifelse(changing_rmssd == "no" & changing_sasa == "no", "yes", "no")) %>% 
  dplyr::select(c(uniprot_id, no_change, any_change, changing_rmssd, changing_sasa, change_in_both)) %>% 
  pivot_longer(cols = no_change:change_in_both, names_to = "type_of_change", values_to = "change") %>% 
  filter(change == "yes") %>% 
  group_by(type_of_change) %>% 
  mutate(change_count = n()) %>% 
  ungroup() %>% 
  dplyr::select(-c(change, uniprot_id)) %>% 
  distinct() %>% 
  mutate(type_of_change = factor(type_of_change, levels = c("no_change", "any_change", "changing_rmssd", "changing_sasa", "change_in_both"))) %>% 
  ggplot(aes(x = type_of_change, y = change_count, fill = type_of_change,)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(labels=c("no_change" = "No change", 
                            "any_change" = "Predicted structural change",                            
                            "changing_rmssd" = "Change in RMSSD",
                            "changing_sasa" = "Change in SASA",
                            "change_in_both" = "Change in RMSSD and SASA")) +
  labs(x = "", y = "Count") +
  theme(axis.text.x = element_text(angle=50, hjust = 1)) +
  scale_fill_manual(values = c("#CBBBBB", "#95ACD3", "#95ACD3", "#95ACD3", "#95ACD3"))


test

pdf("Z:/Viviane/Computational/AF3/AF3_changes_summary_statistics_updated_gmm_based_cutoff_.pdf",
    width = 3,
    height = 5,
    pointsize = 8)
test
dev.off()


############################ Different plot strategy ############################ 

scale_factor = 160
set.seed(124)

rmssd_plot <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>%  
  ggplot(aes(x = pae_weighted_rmssd)) +
  geom_histogram(binwidth = 0.2) +  
  stat_function(fun = function(x) dnorm(x, mean = fit_RMSSD$mu[1], sd = fit_RMSSD$sigma[1]) * scale_factor * fit_RMSSD$lambda[1], color = "#CBBBBB", geom = "line", linewidth = 0.5) +
  stat_function(fun = function(x) dnorm(x, mean = fit_RMSSD$mu[2], sd = fit_RMSSD$sigma[2]) * scale_factor * fit_RMSSD$lambda[2], color = "#95ACD3", geom = "line", linewidth = 0.5) +
  scale_x_continuous(limits = c(-1, 55)) +
  geom_vline(xintercept = sorted_rmssd_vector[237], linetype = "dashed") +
  #geom_vline(xintercept = 0.8, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "RMSSD (Å)", y = "Count")

pdf("Z:/Viviane/Computational/AF3/AF3_pae_weighted_RMSSD_histogram_with_densities_new_cutoff_.pdf",
    width = 5,
    height = 2,
    pointsize = 15)
rmssd_plot
dev.off()

scale_factor = 2940
set.seed(124)

sasa_plot <- all_kinases_pae_weighted_rmssd_and_surface_stu_new_cutoff %>%  
  ggplot(aes(x = abs(change_in_residue_surface))) +
  geom_histogram(binwidth = 3) +
  stat_function(fun = function(x) dnorm(x, mean = fit_SASA$mu[1], sd = fit_SASA$sigma[1]) * scale_factor * fit_SASA$lambda[1], color = "#CBBBBB", geom = "line", linewidth = 0.5) +
  stat_function(fun = function(x) dnorm(x, mean = fit_SASA$mu[2], sd = fit_SASA$sigma[2]) * scale_factor * fit_SASA$lambda[2], color = "#95ACD3", geom = "line", linewidth = 0.5) +
  scale_x_continuous(limits = c(-1, 750)) +
  geom_vline(xintercept = sorted_SASA_vector[280], linetype = "dashed") +
  #geom_vline(xintercept = 15, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = expression("Absolute change in SASA at ATP-binding site (Å"^2*")"), y = "Count")

pdf("Z:/Viviane/Computational/AF3/AF3_sasa_histogram_with_densities_new_cutoff_.pdf",
    width = 5,
    height = 2,
    pointsize = 15)
sasa_plot
dev.off()

