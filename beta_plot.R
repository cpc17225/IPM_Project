library(tidyverse)

ipm_parms <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_comp.RDS")
ipm_parm <- do.call(rbind, ipm_parms)
ipm_parm <- data.frame(Vital_rate_data = rownames(ipm_parm), ipm_parm, row.names = NULL)

betas <- ipm_parm %>% 
  mutate(
    vital_rate = case_when(
      str_starts(Vital_rate_data, "g") ~ "Growth",
      str_starts(Vital_rate_data, "s") ~ "Survival",
      str_starts(Vital_rate_data, "r") ~ "Reproduction",
      str_starts(Vital_rate_data, "f") ~ "Flowering",
      TRUE ~ "Other"
    )
  )

betas[24:26, "vital_rate"] <- "Recruit"
betas <- betas[-c(27:33),]

beta <- betas %>% 
  mutate(
    climate = case_when(
      str_detect(Vital_rate_data, "sum") ~ "Summer temperature",
      str_detect(Vital_rate_data, "spr") ~ "Spring temperature",
      str_detect(Vital_rate_data, "snw") ~ "Snowpack",
      str_detect(Vital_rate_data, "sm") ~ "Snowmelt",
      TRUE ~ "Other"
    )
  ) %>% 
  filter(vital_rate != "Other") %>% 
  filter(climate != "Other") %>% 
  rename("Estimate" = ipm_parm) %>% 
  mutate(
    lag = str_sub(Vital_rate_data, -1)
  )

beta[10:11, "lag"] = "0,1"
beta[1:2, "lag"] = "2,3"


ggplot(beta, aes(x = climate, y = Estimate, fill = vital_rate, group = Vital_rate_data)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c(
    "Growth"       = "darkseagreen4",
    "Recruit"      = "firebrick",
    "Reproduction" = "steelblue",
    "Survival"     = "chocolate"
  )) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


beta_plot <- beta %>%
  mutate(
    effect_dir = case_when(
      Estimate > 0 ~ "Positive (+)",
      Estimate < 0 ~ "Negative (-)",
      TRUE ~ "No Effect (0)"
    ),
    climate_clean = str_replace_all(climate, "_", " ")
  )

ggplot(beta_plot, aes(x = climate_clean, y = vital_rate, fill = effect_dir)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = paste0("Lag ", lag)), 
    color = "white", 
    size = 4.5,           
    lineheight = 0.8,      
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Positive (+)" = "#2E7D32", # Green
      "Negative (-)" = "#C62828"  # Red
    )
  ) +
  labs(
    x = "Climate Variable",
    y = "Vital Rate",
    fill = "Climate Effect"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    legend.position = "top"
  )

ggsave("Beta_plot.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)



survival_beta <- beta %>% 
  filter(vital_rate == "Survival")

ggplot(survival_beta, aes(x = climate, y = Estimate, fill = climate)) +
  geom_col(width = 0.7) +
  labs(x = "Climate Variable", y = "Slope Estimate for Survival Rate") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c(
    "Snowpack" = "lightblue",
    "Spring temperature" = "springgreen4",
    "Summer temperature" = "red3"
  )) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
  
