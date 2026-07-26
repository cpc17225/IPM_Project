### Packages and loading data----

library(glmmTMB)
library(tidyverse)
library(DHARMa)

full_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/full_data.rds")
full_data_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/full_data_rep.rds")
Tg_climate <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_rds.rds")
Tg_climate_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_rep_rds.rds")


### Survival analysis----

model_surv_LLL_c <- glmmTMB(survival ~ log(LLL) + snowmelt + snowpack,
                            family = binomial,
                            data = full_data)
summary(model_surv_LLL_c)



### Growth analysis----
model_growth_LLL_c <- glmmTMB(log(LLL) ~ log(LLLprev) + snowmelt + snowpack,
                              family = t_family(link = "identity"),
                              dispformula = ~log(LLLprev),
                              data = full_data)
summary(model_growth_LLL_c)
plot(simulateResiduals(model_growth_LLL_c))


# Data visualization with climate
ggplot(full_data, aes(x = snowmelt, y = log(LLL) - log(LLLprev))) +
  geom_point() + 
  geom_smooth(method = "lm")
ggplot(full_data, aes(x = snowpack, y = log(LLL) - log(LLLprev)))+
  geom_point()+
  geom_smooth(method = "lm")


### Reproductive analysis----

model_rep_LLL_c <- glmmTMB(rep ~ log(LLLprev) + snowpack + snowmelt,
                           family = binomial,
                           data = full_data)
summary(model_rep_LLL_c)


### Number of flowers of size z, given reproduction

model_flower_LLL_c <- glmmTMB(NumFlowers ~ log(LLLprev) + snowmelt,
                              family = nbinom1(),
                              data = full_data_rep)
summary(model_flower_LLL_c)



### Lagged  effects----


#Survival
surv_models <- vector("list", 5)
names(surv_models) <- paste0("lag_", 0:4)

for(i in 0:4){
  
  # Build climate variable names dynamically
    snow_pack <- paste0("snowpack_lag", i)
    snow_melt <- paste0("snowmelt_lag", i)
  
  # Build formula dynamically
  form <- as.formula(paste("survival ~ log(LLL) + I(log(LLL)^2) + ", snow_pack, "+", snow_melt))
  
  # Fit model and store
  surv_models[[paste0("lag_", i)]] <- glmmTMB(form,
                                              family = binomial,
                                              data = Tg_climate)
}

# Extract AIC for each model
aic_values <- sapply(surv_models, AIC)
print(round(aic_values, 2))

# Identify best model
best_lag  <- names(which.min(aic_values))
best_model <- surv_models[[best_lag]]

cat("\nBest model:", best_lag, "\n")
cat("AIC:", round(min(aic_values), 2), "\n\n")
#Best model is lag 3
#However, model with lag 2 has all significant covariates (in lag 3, snow melt isn't significant)
#Sign flips on snowpack beta value
#Snowmelt beta value is >> snowpack in both models
#DHARMa plots look ver similar
summary(best_model)
summary(surv_models$lag_2)
plot(simulateResiduals(best_model))
plot(simulateResiduals(surv_models$lag_2))
plot(jitter(Tg_climate$survival) ~ Tg_climate$snowpack_lag3)


## Flower model
flower_models <- vector("list", 5)
names(flower_models) <- paste0("lag_", 0:4)

for(i in 0:4){
  snow_pack <- paste0("snowpack_lag", i)
  snow_melt <- paste0("snowmelt_lag", i)
  form <- as.formula(paste("rep ~ log(LLLprev) +", snow_melt,  "+" ,snow_pack))
  flower_models[[paste0("lag_", i)]] <- glmmTMB(form,
                                                family = binomial,
                                                data = Tg_climate)
}

flower_aic <- sapply(flower_models, AIC)
print(round(flower_aic, 2))
cat("Best flowering lag:", names(which.min(flower_aic)), "\n")
#Best model is lag_4
#In lag 4, snowpack is significant
#In lag 2, snowmelt is significant
#Different lags for different climate variables
summary(flower_models$lag_4)
plot(simulateResiduals(flower_models$lag_4))
summary(flower_models$lag_1)

#Visualization of climatic effects on flowering
Tg_climate %>% 
  mutate(size_bin = cut(snowpack4_z, breaks = 12)) %>% 
  group_by(size_bin) %>% 
  summarise(rep_rate = mean(rep), 
            rep_var = var(rep),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = rep_rate, size = n)) +
  geom_point(aes(size = n, color = rep_var)) +
  ylim(0, .15)

ggplot(Tg_climate, aes(snowpack4_z, rep)) +
  geom_jitter(height = 0.03, alpha = 0.15) +
  geom_smooth(method = "gam",
              method.args = list(family = binomial),
              formula = y ~ s(x),
              color = "red") +
  coord_cartesian(ylim = c(0,0.2))

# Growth
growth_models <- vector("list", 5)
names(growth_models) <- paste0("lag_", 0:4)

for(i in 0:4){
  snow_pack <- paste0("snowpack_lag", i)
  snow_melt <- paste0("snowmelt_lag", i)
  form <- as.formula(paste("log(LLL) ~ log(LLLprev) +",
                           snow_pack, "+", snow_melt))
  growth_models[[paste0("lag_", i)]] <- glmmTMB(form,
                                                family  = t_family(link = "identity"),
                                                dispformula = ~ log(LLLprev),
                                                data    = Tg_climate)
}

growth_aic <- sapply(growth_models, AIC)
print(round(growth_aic, 2))
cat("Best growth lag:", names(which.min(growth_aic)), "\n")
#Best model is lag_0
#No other AIC value is close
#Beta values quite small for growth
summary(growth_models$lag_0)


# Number of flowers
num_flower_models <- vector("list", 5)
names(num_flower_models) <- paste0("lag_", 0:4)

for(i in 0:4){
  snow_pack <- paste0("snowpack_lag", i)
  snow_melt <- paste0("snowmelt_lag", i)
  form <- as.formula(paste("NumFlowers ~ log(LLLprev) +",  snow_melt))
  num_flower_models[[paste0("lag_", i)]] <- glmmTMB(form,
                                                family  = nbinom1(),
                                                data    = Tg_climate_rep)
}
num_flower_aic <- sapply(num_flower_models, AIC)
print(round(num_flower_aic, 2))
cat("Best growth lag:", names(which.min(num_flower_aic)), "\n")
#Best model if lag_4, but climatic effects not significant
#Model is unstable when including snowpack
summary(num_flower_models$lag_4)




### Lagged models (snowpack)----
#Survival
surv_models_p <- vector("list", 5)
names(surv_models_p) <- paste0("lag_", 0:4)

for(i in 0:4){
  
  # Build climate variable names dynamically
  snow_pack <- paste0("snowpack_lag", i)
  
  # Build formula dynamically
  form <- as.formula(paste("survival ~ log(LLL) + I(log(LLL)^2) +", snow_pack))
  
  # Fit model and store
  surv_models_p[[paste0("lag_", i)]] <- glmmTMB(form,
                                              family = binomial,
                                              data = Tg_climate)
}
surv_aic <- sapply(surv_models_p, AIC)
print(surv_aic)
min(surv_aic)
#Still best lag is 3
summary(surv_models_p$lag_3)


#Growth models with only snowpack
growth_models_p <- vector("list", 5)
names(growth_models_p) <- paste0("lag_", 0:4)

for(i in 0:4){
  snow_pack <- paste0("snowpack_lag", i)
  form <- as.formula(paste("log(LLL) ~ log(LLLprev) +",
                           snow_pack))
  growth_models_p[[paste0("lag_", i)]] <- glmmTMB(form,
                                                family  = t_family(link = "identity"),
                                                dispformula = ~ log(LLLprev),
                                                data    = Tg_climate)
}
growth_aic_p <- sapply(growth_models_p, AIC)
print(growth_aic_p)
summary(growth_models_p$lag_4)



#Survival of recruits
#Visualization since distribution and range is changing
Tg_data %>% 
  filter(!is.na(LLL)) %>%
  mutate(lgLLL = log(LLL)) %>% 
  filter(lgLLL < 4.3) %>% 
  mutate(size_bin = cut(lgLLL, breaks = 12)) %>% 
  group_by(size_bin) %>%
  summarise(surv_rate = mean(survival), 
            surv_var = var(survival),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = surv_rate, size = n)) +
  geom_point(aes(size = n, color = surv_var)) +
  ylim(0, 1)


surv_models_recruit <- vector("list", 5)
names(surv_models_recruit) <- paste0("lag_", 0:4)

#No parabolic term from visualization
for(i in 0:4){
  
  # Build climate variable names dynamically
  snow_pack <- paste0("snowpack_lag", i)
  
  # Build formula dynamically
  form <- as.formula(paste("survival ~ log(LLL) + ", snow_pack))
  
  # Fit model and store
  surv_models_recruit[[paste0("lag_", i)]] <- glmmTMB(form,
                                              family = binomial,
                                              data = Tg_climate %>% filter(log(LLL)<3.13))
}
surv_recruit_AIC = sapply(surv_models_recruit, AIC)
print(surv_recruit_AIC)
summary(surv_models_recruit$lag_4)
#Still lag 3 is best model
#Use recruit distributions from GMM and do survival analysis on those
#So filtering isn't arbitrary




### Explicit lagged models----

model_surv_LLL_c3 <- glmmTMB(survival ~ log(LLL) + I(log(LLL)^2) + snowpack3_z,
                             family = binomial,
                             data = Tg_climate)
summary(model_surv_LLL_c3)

model_rep_LLL_c4 <- glmmTMB(rep ~ log(LLLprev) + snowpack4_z,
                           family = binomial,
                           data = Tg_climate)
summary(model_rep_LLL_c4)

#Not significant (despite visualization that it looks significant)
#model_rep_LLL_c4_poly <- glmmTMB(
#  rep ~ log(LLLprev) + snowpack4_z + I(snowpack4_z^2),
#  family = binomial,
#  data = Tg_climate)
#summary(model_rep_LLL_c4_poly)

#Growth model. After taking out correlated covariates, snowpack is not
#Significant. Same as null regression
model_growth_LLL <- glmmTMB(log(LLL) ~ log(LLLprev),
                               family = t_family(link = "identity"),
                               dispformula = ~ log(LLLprev),
                               data = Tg_climate)
summary(model_growth_LLL)


#Number of flowers model
#No climate data significant
model_flower_LLL <- glmmTMB(NumFlowers ~ log(LLLprev),
                            family = nbinom1(),
                            data = Tg_climate_rep)
summary(model_flower_LLL)



### Extracting coefficients----

#Growth Parameters
grow_cond_c <- fixef(model_growth_LLL)$cond
grow_disp_c <- fixef(model_growth_LLL)$disp
grow_df_c   <- exp(family_params(model_growth_LLL))

#Survival parameters
surv_cond_c <- fixef(model_surv_LLL_c3)$cond

#Reproduction parameters
rep_cond_c  <- fixef(model_rep_LLL_c4)$cond

#Flower parameters
flow_cond_c <- fixef(model_flower_LLL)$cond

# 5. Compile into the master list for ipmr
ipmr_params_c <- list(
  # Growth
  g_int_c      = as.numeric(grow_cond_c["(Intercept)"]),
  g_slope_c    = as.numeric(grow_cond_c["log(LLLprev)"]),
  d_int_c      = as.numeric(grow_disp_c["(Intercept)"]),
  d_slope_c    = as.numeric(grow_disp_c["log(LLLprev)"]),
  t_df_c       = as.numeric(grow_df_c),
  
  # Survival
  s_int_c      = as.numeric(surv_cond_c["(Intercept)"]),
  s_slope_c    = as.numeric(surv_cond_c["log(LLL)"]),
  s_quad_c     = as.numeric(surv_cond_c["I(log(LLL)^2)"]),
  s_slope_sp    = as.numeric(surv_cond_c["snowpack3_z"]),
  
  # Probability of Reproduction
  r_int_c      = as.numeric(rep_cond_c["(Intercept)"]),
  r_slope_c    = as.numeric(rep_cond_c["log(LLLprev)"]),
  r_slope_sp   = as.numeric(rep_cond_c["snowpack4_z"]),
  
  # Number of Flowers
  f_int_c      = as.numeric(flow_cond_c["(Intercept)"]),
  f_slope_c    = as.numeric(flow_cond_c["log(LLLprev)"]),
  
  #Germination (estimation),
  g_est = 0.02,
  
  #Detection probability (semi-estimated, value from GMM)
  det_p = 0.45,
  
  #Size distribution for pathway 1 recruits
  f1_mean = 3.125133,
  f1_sd = 0.6612659,
  
  #Size distribution for pathway 2 recruits
  f2_mean = 4.216923,
  f2_sd = 0.3048252,
  
  #Survival for seedlings/juveniles
  s_SB = 0.6*0.99,
  s1 = 0.5,
  s2 = 0.4,
  
  #Number of seeds
  num_seeds = 300
)

saveRDS(ipmr_params_c, file = "ipmr_parms_c.rds")
