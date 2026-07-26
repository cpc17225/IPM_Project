### Packages and loading data----

library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(readr)

Tg_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")
Tg_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_rep.rds")
Climate_data <- read_csv("C:/Users/Owner/Downloads/barr_snowmelt_date_2022.csv")
Climate_data <- Climate_data %>% 
  rename("Year" = year)


### Survival analysis----
model_surv_comp <- glmmTMB(survival ~ comp_size + I(comp_size^2),
                           family = binomial,
                           data = Tg_data)
summary(model_surv_comp)
plot(simulateResiduals(model_surv_comp))

ggplot(Tg_data, aes(x = comp_size, y = survival)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = y ~ poly(x,2))



### Growth analysis----
model_growth_comp <- glmmTMB(comp_size ~ comp_size_prev,
                             dispformula = ~ comp_size_prev,
                             family = t_family(link = "identity"),
                             data = Tg_data)
summary(model_growth_comp)
plot(simulateResiduals(model_growth_comp))

nu <- as.numeric(family_params(model_growth_comp))
print(nu)
#df = 6.460343 < 30  => that t_distribution is a better fit and
#data has leptokurtosis
#have to implement t-disribution into growth kernel



### Reproductive analysis----
model_rep_comp <- glmmTMB(rep ~ comp_size_prev,
                          family = binomial,
                          data = Tg_data)
summary(model_rep_comp)

ggplot(Tg_data, aes(x = comp_size_prev, y = rep)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"))



### Flower analysis----
model_flower_comp <- glmmTMB(NumFlowers ~ comp_size_prev,
                             family = poisson,
                             data = Tg_rep)
summary(model_flower_comp)
plot(simulateResiduals(model_flower_comp))
diagnose(model_flower_comp)

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
model_recruit_size <- lm(first_size_comp ~ snowmelt_lag0 + snowmelt_lag1,
                              data = recruit_climate)
summary(model_recruit_size)
#checked for nonconstant variance and normality of residuals, no issues found
#correlation between snowmelt and snowmelt with lag of t=1 is low
cor(climate_lagged$snowmelt_lag0, climate_lagged$snowmelt_lag1, use = "complete.obs")
plot(model_recruit_size)

#coefficients for recruit distribution
recruits_mean <- mean(recruits$first_size_comp)
recruits_sd <- sd(recruits$first_size_comp)

#cofficients for recruit distribution with climate
f1_int <- coef(model_recruit_size)["(Intercept)"]
f1_sm0 <- coef(model_recruit_size)["snowmelt_lag0"]
f1_sm1 <- coef(model_recruit_size)["snowmelt_lag1"]
f1_sd <- summary(model_recruit_size)$sigma



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
  # Growth
  g_int      = as.numeric(grow_cond["(Intercept)"]),
  g_slope    = as.numeric(grow_cond["comp_size_prev"]),
  d_int      = as.numeric(grow_disp["(Intercept)"]),
  d_slope    = as.numeric(grow_disp["comp_size_prev"]),
  t_df       = as.numeric(grow_df),
  
  # Survival
  s_int      = as.numeric(surv_cond["(Intercept)"]),
  s_slope    = as.numeric(surv_cond["comp_size"]),
  s_quad     = as.numeric(surv_cond["I(comp_size^2)"]),
  
  # Probability of Reproduction
  r_int      = as.numeric(rep_cond["(Intercept)"]),
  r_slope    = as.numeric(rep_cond["comp_size_prev"]),
  
  # Number of Flowers
  f_int      = as.numeric(flow_cond["(Intercept)"]),
  f_slope    = as.numeric(flow_cond["comp_size_prev"]),
  
  #Germination (estimation),
  g_est = 0.02,
  
  #Size distribution for recruits
  f1_mean = recruits_mean,
  f1_sd = recruits_sd,
  
  #Survival for seedlings/juveniles
  s_SB = 0.75,
  s1 = 0.75*0.75,
  
  #Number of seeds
  num_seeds = 300
)

saveRDS(ipmr_params, file = "ipmr_parms.RDS")
