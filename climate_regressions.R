### Packages and loading data----

library(MuMIn)
library(car)
library(glmmTMB)

Tg_data_climate <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_clean_scale.rds")

options(na.action = "na.fail")

global_surv <- glmmTMB(
  survival ~ comp_size + I(comp_size^2) + 
    spring.mean.temp_lag0 + spring.mean.temp_lag1 + spring.mean.temp_lag2 + spring.mean.temp_lag3 + spring.mean.temp_lag4 +
    summer.mean.temp_lag0 + summer.mean.temp_lag1 + summer.mean.temp_lag2 + summer.mean.temp_lag3 + summer.mean.temp_lag4 +
    snowpack_lag0 + snowpack_lag1 + snowpack_lag2 + snowpack_lag3 + snowpack_lag4,
  family = binomial,
  data = Tg_data_climate
)

surv_dredge <- dredge(
  global_surv,
  rank = "AICc",
  fixed = c("cond(comp_size)", "cond(I(comp_size^2))"),
  m.lim = c(2,5)
)
