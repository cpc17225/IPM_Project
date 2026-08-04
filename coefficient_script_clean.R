### Packages and loading data----

library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(readr)
library(performance)


Tg_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")
Tg_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_rep.rds")
Climate_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Climate_data_temp.rds")
Climate_data <- Climate_data %>% 
  rename("Year" = year)
Tg_Climate <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_clean_scale_full_data.rds")
Tg_Climate_Flower <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_clean_scale_flower_data.rds")

#climate data for ipmr
ipm_climate_data <- Climate_data %>% 
  filter(Year >=1979) %>% 
  mutate(scaled_snowmelt = scale(snowmelt)) %>% 
  mutate(scaled_snowpack = scale(snowpack)) %>% 
  mutate(scaled_summer_temp = scale(summer.mean.temp)) %>% 
  mutate(scaled_spring_temp = scale(spring.mean.temp))



### Survival analysis----
#model selection from dredging
model_surv_comp <- glmmTMB(survival ~ comp_size + I(comp_size^2) +
                           snowpack_lag3 + spring.mean.temp_lag4 + summer.mean.temp_lag1,
                           family = binomial,
                           data = Tg_Climate,
                           na.action = na.omit)
summary(model_surv_comp)
plot(simulateResiduals(model_surv_comp))
check_collinearity(model_surv_comp)


ggplot(Tg_data) +
  geom_smooth(aes(x = comp_size, y = survival, color = "Survival"),
              method = "glm",
              method.args = list(family = "binomial"),
              formula = y ~ poly(x,2),
              fullrange = FALSE,
              linewidth = 1.2) +
  geom_smooth(aes(x = comp_size_prev, y = rep, color = "Reproduction"),
              method = "glm",
              method.args = list(family = "binomial")) +
  theme_classic(base_size = 11) +
  labs(y = "Probability", x = "Composite size", color = "Vital Rate") +
  scale_color_manual(values = c(
    Survival = "firebrick",
    Reproduction = "steelblue"
  )) +
  theme(legend.position = c(0.3, 0.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
  

ggsave("Survival_vital_rate.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)



### Growth analysis----
model_growth_comp <- glmmTMB(comp_size ~ comp_size_prev +
                               spring.mean.temp_lag1 + spring.mean.temp_lag2 + summer.mean.temp_lag2,
                             dispformula = ~ comp_size_prev,
                             family = t_family(link = "identity"),
                             data = Tg_Climate,
                             na.action = na.omit)
summary(model_growth_comp)
plot(simulateResiduals(model_growth_comp))
check_collinearity(model_growth_comp)

nu <- as.numeric(family_params(model_growth_comp))
print(nu)
#df = 6.460343 < 30  => t_distribution is a better fit and
#data has leptokurtosis
#have to implement t-distribution into growth kernel

ggplot(Tg_Climate, aes(x = comp_size_prev, y = comp_size)) +
  geom_point()




### Reproductive analysis----
model_rep_comp <- glmmTMB(rep ~ comp_size_prev +
                            snowpack_lag4 + summer.mean.temp_lag0 + spring.mean.temp_lag0,
                          family = binomial,
                          data = Tg_Climate,
                          na.action = na.omit)
summary(model_rep_comp)
plot(simulateResiduals(model_rep_comp))
check_collinearity(model_rep_comp)


ggplot(Tg_data, aes(x = comp_size_prev, y = rep)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"))



### Flower analysis----
model_flower_comp <- glmmTMB(NumFlowers ~ comp_size_prev,
                             family = poisson,
                             data = Tg_Climate_Flower,
                             na.action = na.omit)
summary(model_flower_comp)
plot(simulateResiduals(model_flower_comp))

ggplot(Tg_rep, aes(x = comp_size_prev, y = NumFlowers))+
  geom_point()



### Recruit distribution----

#Getting first-year plants and their size
first_appearance <- Tg_data%>%
  group_by(Tag) %>% 
  mutate(Tag = Tag) %>% 
  summarize(first_year = min(Year),
            first_size_comp = comp_size[Year == min(Year)]) %>% 
  filter(first_year > 1980) %>% 
  na.omit()

#four individuals identified  as true seedlings in demographic dataset
#tags: 706, 1031A, 345D, 1082
#sizes of these seedling
z0_known <- c(-3.1606430, -1.941617, -2.69576, -2.695760)
mu_z0 <- mean(z0_known)
s_z0  <- sd(z0_known)
n_z0  <- length(z0_known)

#coefficients from growth model
grow_cond <- fixef(model_growth_comp)$cond
grow_disp <- fixef(model_growth_comp)$disp

#growth model parms
beta_0  <- as.numeric(grow_cond["(Intercept)"])    # Intercept
beta_1  <- as.numeric(grow_cond["comp_size_prev"])    # Slope
gamma_0 <- as.numeric(grow_disp["(Intercept)"])   # Disp Intercept
gamma_1 <- as.numeric(grow_disp["comp_size_prev"])    # Disp Slope

#degrees of freedom of t-distrbution for growth model
nu_growth = nu

#percentile of growth
p <- 0.95

#quantiles
Q_z0     <- qt(p, df = n_z0 - 1)
Q_growth <- qt(p, df = nu_growth)

#growth functions
calc_mu    <- function(z) beta_0 + beta_1 * z
calc_sigma <- function(z) exp(gamma_0 + gamma_1 * z)

#recursive calculations for max growth after 2 years
z0_max <- mu_z0 + Q_z0 * s_z0
z1_max <- calc_mu(z0_max) + Q_growth * calc_sigma(z0_max)

cat("Upper bound starting size (z0_max):", z0_max, "\n")
cat("Maximum size after Year 1 (z1_max):", z1_max, "\n")
#max size for recruits entering population:
#z1 = -0.03082991

#filtering to only recruits under max growth
recruits <- first_appearance %>% 
  filter(first_size_comp<z1_max)

full_recruit_data <- first_appearance %>% 
  mutate(Classification = if_else(first_size_comp<z1_max, "Recruit", "Juvenile/Adult"))

#histograms comparing first-detected and recruits
ggplot(full_recruit_data, aes(x = first_size_comp, color = Classification, fill = Classification)) +
  geom_histogram(binwidth = 0.2) +
  scale_color_manual(values = c(
    "Recruit" = "grey70",
    "Juvenile/Adult" = "black"
  )) + 
  scale_fill_manual(values = c(
    "Recruit" = "steelblue",
    "Juvenile/Adult" = "grey50"
  )) +
  theme_minimal(base_size = 11) +
  labs(x = "Composite Size", y = "Count") +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave("Recruits_distribution.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)


#adding climate data
#lagging climate data for a max of four years
climate_lagged <- Climate_data %>% 
  arrange(Year)
for (i in 0:4) {
  climate_lagged[[paste0("snowmelt_lag", i)]] <- dplyr::lag(climate_lagged$snowmelt, i)
  climate_lagged[[paste0("snowpack_lag", i)]] <- dplyr::lag(climate_lagged$snowpack, i)
}

#left joining
recruit_climate <- recruits %>% 
  rename("Year" = first_year) %>% 
  left_join(climate_lagged, by = "Year")

#recruit size distribution model
model_recruit_size <- lm(first_size_comp ~ scale(snowmelt_lag0) + scale(snowmelt_lag1),
                              data = recruit_climate)
summary(model_recruit_size)
#checked for non-constant variance and normality of residuals, no issues found
#correlation between snowmelt and snowmelt with lag of t=1 is low
cor(climate_lagged$snowmelt_lag0, climate_lagged$snowmelt_lag1, use = "complete.obs")
plot(model_recruit_size)

#coefficients for recruit distribution
recruits_mean <- mean(recruits$first_size_comp)
recruits_sd <- sd(recruits$first_size_comp)

#coefficients for recruit distribution with climate
recruits_int_c <- coef(model_recruit_size)["(Intercept)"]
recruits_slope_sm0 <- coef(model_recruit_size)["scale(snowmelt_lag0)"]
recruits_slope_sm1 <- coef(model_recruit_size)["scale(snowmelt_lag1)"]
recruits_sd_c <- summary(model_recruit_size)$sigma



### Survival coefficients----

#cannot estimate seedling survival
#can estimate survival from seedling -> recruit bank
#using true seedlings from recruit analysis

recruit_survival <- predict(
  model_surv_comp,
  newdata = data.frame(comp_size = mu_z0,
                       snowpack_lag3 = 0,
                       spring.mean.temp_lag4 = 0,
                       summer.mean.temp_lag1 = 0),
  type = "response"
)




### Extracting coefficients----

#growth Parameters
grow_cond <- fixef(model_growth_comp)$cond
grow_disp <- fixef(model_growth_comp)$disp
grow_df   <- nu

#survival parameters
surv_cond <- fixef(model_surv_comp)$cond

#reproduction parameters
rep_cond  <- fixef(model_rep_comp)$cond

#flower parameters
flow_cond <- fixef(model_flower_comp)$cond



#compile into the master list for ipmr
ipmr_params_comp <- list(
  # Growth (Student-t with size-dependent dispersion)
  g_int        = as.numeric(grow_cond["(Intercept)"]),
  g_slope      = as.numeric(grow_cond["comp_size_prev"]),
  g_spr_lag1   = as.numeric(grow_cond["spring.mean.temp_lag1"]),
  g_spr_lag2   = as.numeric(grow_cond["spring.mean.temp_lag2"]),
  g_sum_lag2   = as.numeric(grow_cond["summer.mean.temp_lag2"]),
  
  d_int        = as.numeric(grow_disp["(Intercept)"]),
  d_slope      = as.numeric(grow_disp["comp_size_prev"]),
  t_df         = as.numeric(grow_df),
  
  # Survival
  s_int        = as.numeric(surv_cond["(Intercept)"]),
  s_slope      = as.numeric(surv_cond["comp_size"]),
  s_quad       = as.numeric(surv_cond["I(comp_size^2)"]),
  s_snw_lag3   = as.numeric(surv_cond["snowpack_lag3"]),
  s_spr_lag4   = as.numeric(surv_cond["spring.mean.temp_lag4"]),
  s_sum_lag1   = as.numeric(surv_cond["summer.mean.temp_lag1"]),
  
  # Probability of Reproduction
  r_int        = as.numeric(rep_cond["(Intercept)"]),
  r_slope      = as.numeric(rep_cond["comp_size_prev"]),
  r_snw_lag4   = as.numeric(rep_cond["snowpack_lag4"]),
  r_sum_lag0   = as.numeric(rep_cond["summer.mean.temp_lag0"]),
  r_spr_lag0   = as.numeric(rep_cond["spring.mean.temp_lag0"]),
  
  # Number of Flowers
  f_int        = as.numeric(flow_cond["(Intercept)"]),
  f_slope      = as.numeric(flow_cond["comp_size_prev"]),
  
  # Germination (estimation)
  g_est        = 0.01,
  
  # Size distribution for recruits
  rec_int_c =    as.numeric(coef(model_recruit_size)["(Intercept)"]),
  rec_slope_sm0 = as.numeric(coef(model_recruit_size)["scale(snowmelt_lag0)"]),
  rec_slope_sm1 = as.numeric(coef(model_recruit_size)["scale(snowmelt_lag1)"]),
  rec_sd_c = as.numeric(summary(model_recruit_size)$sigma),
  
  # Survival for seedlings/juveniles
  s_SB         = 0.4,
  s1           = recruit_survival,
  
  # Number of seeds
  num_seeds    = 300,
  
  # Mean climate (0 since variables scaled)
  snowmelt_mean  = 0,
  snowpack_mean = 0,
  summer_temp_mean = 0,
  spring_temp_mean = 0
)

saveRDS(ipmr_params_comp, file = "ipmr_parms_comp.RDS")

### Size variables----

#omega variables

L <- min(Tg_Climate$comp_size, na.rm = TRUE)
#min size = -3.990192
U <- max(Tg_Climate$comp_size, na.rm = TRUE)
#max size = 6.025928




#starting size variables
Tg_data_1979 <- Tg_Climate %>% 
  filter(Year == 1979) %>% 
  pull(comp_size)

#total population size in 1979 = 70 (should be hidden plants in recruit bank)

n_mesh <- 100

mesh_breaks <- seq(L, U, length.out = n_mesh + 1)

emp_hist <- hist(Tg_data_1979, breaks = mesh_breaks, plot = FALSE)
initial_size_vector <- emp_hist$counts

saveRDS(initial_size_vector, file = "initial_size_vector.rds")



#starting variables for 1996
#pull the actual plant sizes from 1996
sizes_1996 <- Tg_Climate %>% 
  filter(Year == 1996) %>% 
  pull(comp_size_prev) 

#bin into your 100 IPM mesh points (using your existing mesh_breaks)
emp_hist_1996 <- hist(sizes_1996, breaks = mesh_breaks, plot = FALSE)
initial_size_vector_1996 <- emp_hist_1996$counts

saveRDS(initial_size_vector_1996, file = "initial_size_vector_1996.rds")

Tg_Climate %>% 
  filter(Year == 1995) %>% 
  pull(comp_size)
