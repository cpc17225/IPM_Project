### Packages and loading data----

library(glmmTMB)
library(tidyverse)
library(DHARMa)

full_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/full_data.rds")
full_data_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/full_data_rep.rds")
Tg_climate <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate.rds")
Tg_climate_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_rep.rds")


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
