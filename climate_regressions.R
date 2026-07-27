### Packages and loading data----

library(MuMIn)
library(performance)
library(glmmTMB)

Tg_data_climate <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_clean_scale.rds")
Tg_data_climate_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_climate_clean_scale_rep.rds")

options(na.action = "na.fail")



### Survival----

#full survival model
global_surv <- glmmTMB(
  survival ~ comp_size + I(comp_size^2) + 
    spring.mean.temp_lag0 + spring.mean.temp_lag1 + spring.mean.temp_lag2 +
    spring.mean.temp_lag3 + spring.mean.temp_lag4 +
    summer.mean.temp_lag0 + summer.mean.temp_lag1 +
    summer.mean.temp_lag2 + summer.mean.temp_lag3 + summer.mean.temp_lag4 +
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

#model selection
head(surv_dredge, 10)
best_surv_model <- get.models(surv_dredge, subset = 1)[[1]]
summary(best_surv_model)
check_collinearity(best_surv_mod)



### Reproduction----

#full reproduction model
global_rep <- glmmTMB(
  rep ~ comp_size_prev + 
    spring.mean.temp_lag0 + spring.mean.temp_lag1 + spring.mean.temp_lag2 +
    spring.mean.temp_lag3 + spring.mean.temp_lag4 +
    summer.mean.temp_lag0 + summer.mean.temp_lag1 + summer.mean.temp_lag2 +
    summer.mean.temp_lag3 + summer.mean.temp_lag4 +
    snowpack_lag0 + snowpack_lag1 + snowpack_lag2 + snowpack_lag3 + snowpack_lag4,
  family = binomial,
  data = Tg_data_climate_rep
)

rep_dredge <- dredge(
  global_rep,
  rank = "AICc",
  fixed = c("cond(comp_size_prev)"),
  m.lim = c(1,4)
)
#model selection
#three models within delta AIC = 2
#"second best" model is most biologically realistic (uses current summer and spring temps)
head(rep_dredge, 10)
best_rep_model <- get.models(rep_dredge, subset = 1)[[1]]
summary(best_rep_model)
second_rep_model <- get.models(rep_dredge, subset = 2)[[1]]
summary(second_rep_model)
check_collinearity(second_rep_model)
third_rep_model <- get.models(rep_dredge, subset = 3)[[1]]
summary(third_rep_model)



### Growth----

#only observations with size at t and t-1
Tg_growth <- Tg_data_climate %>% 
  drop_na(
    comp_size_prev
  )
#full growth model
global_growth <- glmmTMB(
  comp_size ~ comp_size_prev + 
    spring.mean.temp_lag0 + spring.mean.temp_lag1 + spring.mean.temp_lag2 + spring.mean.temp_lag3 + spring.mean.temp_lag4 +
    summer.mean.temp_lag0 + summer.mean.temp_lag1 + summer.mean.temp_lag2 + summer.mean.temp_lag3 + summer.mean.temp_lag4 +
    snowpack_lag0 + snowpack_lag1 + snowpack_lag2 + snowpack_lag3 + snowpack_lag4,
  family = t_family(link = "identity"),
  dispformula = ~ comp_size_prev,
  data = Tg_growth
)

growth_dredge <- dredge(
  global_growth,
  rank = "AICc",
  fixed = c("cond(comp_size_prev)", "disp(comp_size_prev)"),
  m.lim = c(2,5)
)
#model selection
head(growth_dredge, 10)
best_growth_model <- get.models(growth_dredge, subset = 1)[[1]]
summary(best_growth_model)




### Flowering----

#only reproductive individuals
Tg_flower <- Tg_data_climate_rep %>% 
  filter(rep==1)

#full flower model
global_flower <- glmmTMB(
  NumFlowers ~ comp_size_prev + 
    spring.mean.temp_lag0 + spring.mean.temp_lag1 + spring.mean.temp_lag2 + spring.mean.temp_lag3 + spring.mean.temp_lag4 +
    summer.mean.temp_lag0 + summer.mean.temp_lag1 + summer.mean.temp_lag2 + summer.mean.temp_lag3 + summer.mean.temp_lag4 +
    snowpack_lag0 + snowpack_lag1 + snowpack_lag2 + snowpack_lag3 + snowpack_lag4,
  family = poisson(),
  data = Tg_flower
)

flower_dredge <- dredge(
  global_flower,
  rank = "AICc",
  fixed = c("cond(comp_size_prev)"),
  m.lim = c(1,4)
)

#model selection
head(flower_dredge, 10)


### Supported models----

#Survival model
# survival ~ snowpack_lag3 + spring.mean.temp_lag4 + summer.mean.temp_lag1 +
#            comp_size + (comp_size)^2
# family = binomial(logit)

#delta AIC = 5.39


#Reproductive model
# rep ~ snowpack_lag4 + spring.mean.temp_lag0 + summer.mean.temp_lag0 +
#       comp_size_prev
# family = binomial(logit)

#three models withing delta AIC = 2. chose most biologically realistic model


#Growth model
# comp_size ~ spring.mean.temp_lag1 + spring.mean.temp_lag2 + summer.mean.temp_lag2 +
#             comp_size_prev
# dispformula = ~ comp_size_prev
# family = t_family(identity)

#delta AIC = 7.89


#Flower model
# NumFlowers ~ comp_size_prev
# family = poission(log)

#many models within delta AIC = 2
#null model of no climate within this range => best model
