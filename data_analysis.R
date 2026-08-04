## Packages and loading data----

library(glmmTMB)
library(tidyverse)
library(goftest)
library(moments)
library(mgcv)
library(DHARMa)

Tg_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")
Tg_rep <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_rep.rds")


### Survival analysis----

# GLMM for survival analysis using the length of longest leaf (LLL)
model_surv_LLL <- glmmTMB(survival ~ log(LLL),
                          data = Tg_data,
                          family = binomial)
summary(model_surv_LLL)
plot(simulateResiduals(model_surv_LLL))
# Trying a GAM
model_surv_LLL_gam <- gam(survival ~ s(log(LLL)),
                              data = Tg_data,
                              family = binomial)
summary(model_surv_LLL_gam)

# Adding a polynomial term
# This is the best fitting model. p value very small and DHARMa plot looks great
model_surv_LLL_poly <- glmmTMB(survival ~ log(LLL) + I(log(LLL)^2),
                               data = Tg_data,
                               family = binomial)
summary(model_surv_LLL_poly)
plot(simulateResiduals(model_surv_LLL_poly))

# GLMM for survival analysis using the number of leaves
model_surv_num <- glmmTMB(survival ~ numleaves + (1|Tag),
                          data = Tg_data,
                          family = binomial)
summary(model_surv_num)
plot(simulateResiduals(model_surv_num))


## Data visualization with LLL and numleaves for survival

#Using LLL
Tg_data %>% 
  filter(!is.na(LLL)) %>%
  mutate(size_bin = cut(LLL, breaks = 12)) %>% 
  group_by(size_bin) %>% 
  summarise(surv_rate = mean(survival), 
            surv_var = var(survival),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = surv_rate, size = n)) +
  geom_point(aes(size = n, color = surv_var)) +
  ylim(0, 1)
#Using log(LLL)
Tg_data %>% 
  filter(!is.na(LLL)) %>%
  mutate(lgLLL = log(LLL)) %>% 
  mutate(size_bin = cut(lgLLL, breaks = 12)) %>% 
  group_by(size_bin) %>% 
  summarise(surv_rate = mean(survival), 
            surv_var = var(survival),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = surv_rate, size = n)) +
  geom_point(aes(size = n, color = surv_var)) +
  ylim(0, 1)
#Using a gam-style smooth survival curve
ggplot(filter(Tg_data, !is.na(LLL)),
       aes(log(LLL), survival)) +
  geom_jitter(height = 0.03, alpha = 0.15) +
  geom_smooth(method = "gam",
              method.args = list(family = binomial),
              formula = y ~ s(x),
              color = "red")


#Using numleaves
Tg_data %>% 
  mutate(size_bin = cut(numleaves, breaks = 7)) %>% 
  group_by(size_bin) %>% 
  summarise(surv_rate = mean(survival), 
            surv_var = var(survival),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = surv_rate, size = n)) +
  geom_point(aes(size = n, color = surv_var)) +
  ylim(0, 1)



### Growth analysis----

#These models are acting a little odd. Could be outliers, I'm not sure?

#For LLL
model_growth_LLL <- glmmTMB(log(LLL) ~ log(LLLprev),
                            family = gaussian(),
                            dispformula = ~ log(LLLprev),
                            data = Tg_data)
summary(model_growth_LLL)
plot(simulateResiduals(model_growth_LLL))

#For numleaves
model_growth_num <- glmmTMB(log(numleaves) ~ log(numleavesprev),
                            family = gaussian(),
                            dispformula = ~ log(numleavesprev),
                            data = Tg_data)
summary(model_growth_num)
plot(simulateResiduals(model_growth_num))



##Data visualization of growth

#For LLL
ggplot(Tg_data, aes(x = LLLprev, y = LLL)) +
  geom_point()

#For LLL (using log(LLL))
ggplot(Tg_data, aes(x = log(LLLprev), y = log(LLL))) +
  geom_point() +
  xlim(1.5, 5) + 
  ylim(1.5, 5.25) + 
  geom_smooth(method = "lm")

hist(Tg_data$LLLprev, breaks = 20)
hist(log(Tg_data$LLLprev), breaks = 20)
hist(Tg_data$LLL, breaks = 20)
hist(log(Tg_data$LLL), breaks = 20)
  

#For numleaves
ggplot(Tg_data, aes(x = numleavesprev, y = numleaves)) + 
  geom_point() +
  xlim(0,65) + 
  ylim(0,65)

hist(Tg_data$numleavesprev, breaks = 20)
hist(log(Tg_data$numleavesprev), breaks = 20)

hist(Tg_data$numleaves, breaks = 20)
hist(log(Tg_data$numleaves), breaks = 20)


### Reproductive analysis----

## Probability of reproducing at size z

#Using LLL
model_rep_LLL <- glmmTMB(rep ~ log(LLLprev),
                         family = binomial,
                         data = Tg_data)
summary(model_rep_LLL)
plot(simulateResiduals(model_rep_LLL))

#Using numleaves
model_rep_num <- glmmTMB(rep ~ numleavesprev + (1|Tag),
                         family = binomial,
                         data = Tg_data)
summary(model_rep_num)
plot(simulateResiduals(model_rep_num))


## Data visualization of probability of reproduction

#Using LLL
Tg_data %>% 
  filter(LLLprev>60) %>% 
  mutate(size_bin = cut(LLLprev, breaks = 6)) %>% 
  group_by(size_bin) %>% 
  summarise(rep_rate = mean(rep), 
            rep_var = var(rep),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = rep_rate, size = n)) +
  geom_point(aes(size = n, color = rep_var)) +
  ylim(0, 1)

#Using numleaves
Tg_data %>% 
  filter(numleavesprev>8) %>% 
  mutate(size_bin = cut(numleavesprev, breaks = 6)) %>% 
  group_by(size_bin) %>% 
  summarise(rep_rate = mean(rep), 
            rep_var = var(rep),
            n= n()) %>% 
  ggplot(aes(x = size_bin, y = rep_rate, size = n)) +
  geom_point(aes(size = n, color = rep_var)) +
  ylim(0, 1)


## Number of flowers of size z, given reproduction

#Using LLL
model_flower_LLL <- glmmTMB(NumFlowers ~ log(LLLprev),
                            family = nbinom1(),
                            data = Tg_rep)
summary(model_flower_LLL)
plot(simulateResiduals(model_flower_LLL))

#Using numleaves
#This model has some issues, may be from numleaves outlier?
model_flower_num <- glmmTMB(NumFlowers ~ scale(numleavesprev),
                            family = nbinom1(),
                            data = Tg_rep)
summary(model_flower_num)
plot(simulateResiduals(model_flower_num))


## Data visualization of the number of flowers

#Using LLL
ggplot(Tg_rep, aes(x = LLLprev, y = NumFlowers)) +
  geom_point()
#Using log(LLL)
ggplot(Tg_rep, aes(x = log(LLLprev), y = NumFlowers)) +
  geom_point() +
  xlim(4,5.3)

#Using numleaves
ggplot(Tg_rep, aes(x = numleavesprev, y = NumFlowers)) +
  geom_point()






### Diagnostic stuff----


#Cramer-von Mises test for uniform residual distribution in growth model using numleaves
#Significant deviation from uniformity, seen with the bump in the QQ plot
res_num <- simulateResiduals(model_growth_num)
scaled_res_num <- res_num$scaledResiduals
scaled_res_num1 <- na.omit(scaled_res_num)
cvm.test(scaled_res_num1, null = "punif", min = 0, max = 1)


#Skewness and kurtosis measures
resids1 <- residuals(model_growth_num)
skewness(resids1)
kurtosis(resids1)

resids2 <- residuals(model_growth_LLL)
skewness(resids2)
kurtosis(resids2)



#Trying a t distribution for growth

#For LLL
#Best fitting model so far
model_growth_LLL_t <- glmmTMB(log(LLL) ~ log(LLLprev),
                              family = t_family(link = "identity"),
                              dispformula = ~ log(LLLprev),
                              data = Tg_data)
summary(model_growth_LLL_t)
plot(simulateResiduals(model_growth_LLL_t))

#For numleaves
model_growth_num_t <- glmmTMB(numleaves ~ numleavesprev,
                              family = t_family(link = "identity"),
                              dispformula = ~ numleavesprev,
                              data = Tg_data)
summary(model_growth_num_t)
plot(simulateResiduals(model_growth_num_t))



### Population size over time ----

#Getting population size per year
pop_size <- Tg_data %>% 
  group_by(Year) %>% 
  summarize(pop = n()) %>% 
  mutate(logpop = log(pop))

pop_size$Year <- as.numeric(as.character(pop_size$Year))

# Plotting N(t)
ggplot(pop_size, aes(x = Year, y = pop)) +
  geom_point() +
  scale_x_continuous(breaks = seq(1979, 2025, by = 5))
# Plotting log[N(t)]
ggplot(pop_size, aes(x = Year, y = logpop)) +
  geom_point() +
  scale_x_continuous(breaks = seq(1979, 2025, by = 5))

#Getting lambda
model_lambda <- lm(logpop ~ Year,
                  data = pop_size[pop_size$Year>1995, ])
summary(model_lambda)
bounds <- confint(model_lambda, level = 0.95)
lambda <- exp(coef(model_lambda)["Year"])
lambda_bounds <- exp(bounds[2,])

model_lambda_inc <-lm(logpop ~ Year,
                      data = pop_size[pop_size$Year<1996, ])
summary(model_lambda_inc)
lambda_inc <- exp(coef(model_lambda_inc)["Year"])
bounds <- confint(model_lambda_inc, level = 0.95)
lambda_inc_bounds <- exp(bounds[2,])

#Lambda = 0.9527436 
plot(model_lambda)




### Extracting coefficients----

#Growth Parameters
grow_cond <- fixef(model_growth_LLL_t)$cond
grow_disp <- fixef(model_growth_LLL_t)$disp
grow_df   <- exp(family_params(model_growth_LLL_t))

#Survival parameters
surv_cond <- fixef(model_surv_LLL_poly)$cond

#Reproduction parameters
rep_cond  <- fixef(model_rep_LLL)$cond

#Flower parameters
flow_cond <- fixef(model_flower_LLL)$cond

# 5. Compile into the master list for ipmr
ipmr_params <- list(
  # Growth
  g_int      = as.numeric(grow_cond["(Intercept)"]),
  g_slope    = as.numeric(grow_cond["log(LLLprev)"]),
  d_int      = as.numeric(grow_disp["(Intercept)"]),
  d_slope    = as.numeric(grow_disp["log(LLLprev)"]),
  t_df       = as.numeric(grow_df),
  
  # Survival
  s_int      = as.numeric(surv_cond["(Intercept)"]),
  s_slope    = as.numeric(surv_cond["log(LLL)"]),
  s_quad     = as.numeric(surv_cond["I(log(LLL)^2)"]),
  
  # Probability of Reproduction
  r_int      = as.numeric(rep_cond["(Intercept)"]),
  r_slope    = as.numeric(rep_cond["log(LLLprev)"]),
  
  # Number of Flowers
  f_int      = as.numeric(flow_cond["(Intercept)"]),
  f_slope    = as.numeric(flow_cond["log(LLLprev)"]),
  
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

saveRDS(ipmr_params, file = "ipmr_parms.RDS")

ggplot(Tg_data, aes(x = sqrt(numleaves), y = LLL)) + 
  geom_point() +
  ylim(0,175)

plot(sqrt(Tg_data$numleaves), Tg_data$LLL)
points(sqrt(21), 44, col = 2, cex = 2)

     