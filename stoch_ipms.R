### Packages and loading data----

library(ipmr)
library(tidyverse)
library(Metrics)
library(patchwork)
library(mgcv)
library(gratia)

ipm_parms <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_comp.RDS")
initial_size <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/initial_size_vector.RDS")

ipm_parms$s_SB = 0.32

### Helper functions/variables----

#used to help manage t-distributed growth
#these are used in growth distribution parts of the survival kernels
dt_scaled <- function(x, mean, sd, df) {
  dt((x - mean) / sd, df = df) / sd
}

pt_scaled <- function(q, mean, sd, df) {
  pt((q - mean) / sd, df = df)
}

#omega bounds for size variables
#lower bound for size (-0.5 to account for size eviction)
L = -3.990192 - 0.5
#upper bound for size (+0.5 to account for size eviction)
U = 6.025928 + 0.5

#these +-0.5 are arbitrary and don't necessarily follow the 20% rule
#but +-0.5 standard deviations was convenient and seems to work in the ipms



### Null IPM----



### Building kernels----

#defining the structure of the ipms
#general ipm because we have the recruit bank
#we are using a density independent ipm
#stochastic to feed climate data into the ipm
#we are using kernel resampling because we are feeding historical data
#each kernel is resamapled each year

null_stoch_ipm <- init_ipm(sim_gen = "general",
                 di_dd = "di",
                 det_stoch = "stoch",
                 kern_param = "kern")

#growth/survival kernel (P -> P)
null_stoch_ipm <- define_kernel(
  proto_ipm = null_stoch_ipm,
  name = "P_yr",
  family = "CC",
  formula = s * G * (1-pf) * d_compsize,
  #d_compsize is just an integration factor so ipmr treats the formula correctly
  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
               s_snw_lag3 * snowpack_mean +
               s_spr_lag4 * spring_temp_mean +
               s_sum_lag1 * summer_temp_mean),
  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
  mu_g = g_int + g_slope*compsize_1 +
    g_spr_lag1 * spring_temp_mean +
    g_spr_lag2 * spring_temp_mean +
    g_sum_lag2 * summer_temp_mean,
  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * snowpack_mean +
                r_sum_lag0 * summer_temp_mean +
                r_spr_lag0 * spring_temp_mean),
  data_list = ipm_parms,
  states = list(c("compsize")),
  evict_cor = TRUE,
  #uses t-distribution for eviction correction function
  evict_fun = truncated_distributions(fun = "t_scaled",
                                      target = "G"),
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
  
)

#fecundity kernel (Fe -> RB)
null_stoch_ipm <- define_kernel(
  proto_ipm = null_stoch_ipm,
  name = "Fe_yr",
  family = "CD",
  formula = pf * fs * num_seeds * g_est * s_SB,
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * snowpack_mean +
                r_sum_lag0 * summer_temp_mean +
                r_spr_lag0 * spring_temp_mean),
  fs = exp(f_int + f_slope*compsize_1),
  data_list = ipm_parms,
  states = list(c("compsize", "RB")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)

#recruit bank kernel (RB -> P)
null_stoch_ipm <- define_kernel(
  proto_ipm = null_stoch_ipm,
  name = "RB_yr",
  family = "DC",
  formula = s1 * f1_dist * d_compsize,
  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
  mu_f1 = rec_int_c +
    rec_slope_sm0 * snowmelt_mean +
    rec_slope_sm1 * snowmelt_mean,
  sd_f1 = rec_sd_c,
  data_list = ipm_parms,
  states = list(c("RB", "compsize")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)



### Defining states----

#defining where the kernels go
null_stoch_ipm <- define_impl(
  proto_ipm = null_stoch_ipm,
  make_impl_args_list(
    kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
    #midpoint integration rule over 3 distinct states
    int_rule = c(rep("midpoint", 3)),
    state_start = c("compsize", "compsize", "RB"),
    state_end = c("compsize", "RB", "compsize")
  )
)


#defining omega
null_stoch_ipm <- define_domains(
  proto_ipm = null_stoch_ipm,
  # 100 for number of meshpoints to integrate over
  # 100 is the standard, "default" number of meshpoints
  compsize = c(L, U, 100)
)

#defining initial population state
null_stoch_ipm <- define_pop_state(
  proto_ipm = null_stoch_ipm,
  pop_vectors = list(
    #initial size distribution in 1979
    n_compsize = initial_size,
    #recruit bank number is rough estimation
    n_RB = 10
    #empirical size distribution does not empirically account for individuals
    #"hidden" in the recruit bank
  )
)

#generating IPM
null_stoch_ipm <- make_ipm(
  proto_ipm = null_stoch_ipm,
  iterations = 47,
  kernel_seq = 1:47,
  #uses the helper functions
  usr_funs = list(
    dt_scaled = dt_scaled,
    pt_scaled = pt_scaled
  ),
  return_all_envs = TRUE
)

lambda(null_stoch_ipm, log = FALSE)
#IS convergent to asymptotic
is_conv_to_asymptotic(null_stoch_ipm)

#annual lambdas by year
#besides transient dynamics, these approach a constant value
#because we are using a constant climate for this ipm
annual_lambdas_null <- null_stoch_ipm$pop_state$lambda
lambda_df <- data.frame(
  Year = 1979:2025,
  lamdba = annual_lambdas_null,
  log_lambda = log(annual_lambdas_null)
)

lambda_df


#set actual 1979 population size
N_1979 <- 70

#calculate absolute population size over time using cumulative products
absolute_pop <- N_1979 * c(1, cumprod(annual_lambdas_null))

#clean data frame
absolute_pop_df <- data.frame(
  Year = 1979:2026,
  Total_N = absolute_pop
)

#after transient dynamics wash out, we steady out into approximately a linear curve
ggplot(absolute_pop_df, aes(x = Year, y = Total_N)) +
  geom_line()









### Historic climate IPM----



### Additional data needed----

#these parameters include the time-lag corrected climate inputs "baked into"
#each vital rate, so can call each climate variable without having to explicitly
#define the time lags within the ipm
ipm_parms_clim <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipm_parms_climate.RDS")

ipm_parms_clim$s_SB = 0.32

### Building kernels----

clim_stoch_ipm <- init_ipm(sim_gen = "general",
                           di_dd = "di",
                           det_stoch = "stoch",
                           kern_param = "kern")

#growth/survival kernel (P -> P)
clim_stoch_ipm <- define_kernel(
  proto_ipm = clim_stoch_ipm,
  name = "P_yr",
  family = "CC",
  formula = s * G * (1-pf) * d_compsize,
  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
               s_snw_lag3 * sp_surv_yr +
               s_spr_lag4 * spt_surv_yr +
               s_sum_lag1 * sut_surv_yr),
  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
  mu_g = g_int + g_slope*compsize_1 +
    g_spr_lag1 * spt1_grow_yr +
    g_spr_lag2 * spt2_grow_yr +
    g_sum_lag2 * sut_grow_yr,
  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * sp_rep_yr +
                r_sum_lag0 * sut_rep_yr +
                r_spr_lag0 * spt_rep_yr),
  data_list = ipm_parms_clim,
  states = list(c("compsize")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "t_scaled",
                                      target = "G"),
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
  
)

#fecundity kernel (Fe -> RB)
clim_stoch_ipm <- define_kernel(
  proto_ipm = clim_stoch_ipm,
  name = "Fe_yr",
  family = "CD",
  formula = pf * fs * num_seeds * g_est * s_SB,
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * sp_rep_yr +
                r_sum_lag0 * sut_rep_yr +
                r_spr_lag0 * spt_rep_yr),
  fs = exp(f_int + f_slope*compsize_1),
  data_list = ipm_parms_clim,
  states = list(c("compsize", "RB")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)

#recruit bank kernel (RB -> P)
clim_stoch_ipm <- define_kernel(
  proto_ipm = clim_stoch_ipm,
  name = "RB_yr",
  family = "DC",
  formula = s1 * f1_dist * d_compsize,
  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
  mu_f1 = rec_int_c +
    rec_slope_sm0 * sm0_rec_yr +
    rec_slope_sm1 * sm1_rec_yr,
  sd_f1 = rec_sd_c,
  data_list = ipm_parms_clim,
  states = list(c("RB", "compsize")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)



### Defining states----

#defining where the kernels go
clim_stoch_ipm <- define_impl(
  proto_ipm = clim_stoch_ipm,
  make_impl_args_list(
    kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
    int_rule = c(rep("midpoint", 3)),
    state_start = c("compsize", "compsize", "RB"),
    state_end = c("compsize", "RB", "compsize")
  )
)


#defining omega
clim_stoch_ipm <- define_domains(
  proto_ipm = clim_stoch_ipm,
  # 100 for number of meshpoints to integrate over
  compsize = c(L, U, 100)
)

#defining initial population state
clim_stoch_ipm <- define_pop_state(
  proto_ipm = clim_stoch_ipm,
  pop_vectors = list(
    #initial size distribution in 1979
    n_compsize = initial_size,
    #recruit bank number is rough estimation
    n_RB = 10
  )
)

#generating IPM
clim_stoch_ipm <- make_ipm(
  proto_ipm = clim_stoch_ipm,
  iterations = 47,
  kernel_seq = 1:47,
  usr_funs = list(
    dt_scaled = dt_scaled,
    pt_scaled = pt_scaled
  ),
  return_all_envs = TRUE
)

#lambda is not convergent to asymptotic because we are feeding yearly-differing
#climate variables
lambda(clim_stoch_ipm, log = FALSE)


annual_lambda_clim <- clim_stoch_ipm$pop_state$lambda
lambda_df_clim <- lambda_df %>% 
  mutate(lamdba_clim = annual_lambda_clim) %>% 
  mutate("log(lambda_clim)" = log(annual_lambda_clim))
lambda_df_clim

#finding lambda for just 1994-2026 (post-crash)
lambda_df_post_crash <- lambda_df_clim %>% 
  filter(Year>1995)
lambda_post_crash <- exp(mean(log(lambda_df_post_crash$lamdba_clim)))
#after population crash, stochastic lambda = 0.9272315
lambda_post_crash
#from empirical data, should be 0.9527436


### Analysis----
Tg_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")

#set actual 1979 population size
N_1979 <- 70

#calculate absolute population size over time using cumulative products
absolute_pop_clim <- N_1979 * c(1, cumprod(annual_lambda_clim))

#clean data frame
absolute_pop_clim_df <- data.frame(
  Year = 1979:2026,
  Total_N = absolute_pop_clim,
  logN = log(absolute_pop_clim)
)

ggplot(absolute_pop_clim_df, aes(x = Year, y = Total_N)) +
  geom_point()
ggplot(absolute_pop_clim_df, aes(x = Year, y = logN)) +
  geom_point() +
  labs(y = "log[N(t)] from IPM")

pop_size <- Tg_data %>% 
  group_by(Year) %>% 
  summarize(pop = n()) %>% 
  mutate(logpop = log(pop))

pop_size$Year <- as.numeric(as.character(pop_size$Year))

emp_lambda <- lead(pop_size$pop)/pop_size$pop

#plotting N(t) from empirical data
ggplot(pop_size, aes(x = Year, y = pop)) +
  geom_point() +
  scale_x_continuous(breaks = seq(1979, 2025, by = 5))
#plotting log[N(t)] from empirical data
pop_size <- pop_size %>%
  mutate(Regime = ifelse(Year < 1996, "Historical Boom (< 1996)", "Declining Period (1996-2025)"))
ggplot(pop_size, aes(x = Year, y = logpop, color = Regime)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Historical Boom (< 1996)" = "gray60", 
                                "Declining Period (1996-2025)" = "darkred"))+
  scale_x_continuous(breaks = seq(1979, 2025, by = 5)) +
  labs(y = "log[N(t)] from empirical data",
       color = "Demographic Regime")+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.8, 0.85))

ggsave("Figure1_Empirircal_Population.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)



#combining all populations into one dataframe
combined_pop <- data.frame(
  emp_pop = pop_size$pop,
  clim_pop = absolute_pop_clim_df$Total_N,
  null_pop = absolute_pop_df$Total_N,
  Year = 1979:2026
)




### 1996-2026 IPM----

initial_size_1996 <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/initial_size_vector_1996.RDS")


clim_stoch_ipm_1996 <- init_ipm(sim_gen = "general",
                                di_dd = "di",
                                det_stoch = "stoch",
                                kern_param = "kern")


### Building kernels----


#growth/survival kernel (P -> P)
clim_stoch_ipm_1996 <- define_kernel(
  proto_ipm = clim_stoch_ipm_1996,
  name = "P_yr",
  family = "CC",
  formula = s * G * (1-pf) * d_compsize,
  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
               s_snw_lag3 * sp_surv_yr +
               s_spr_lag4 * spt_surv_yr +
               s_sum_lag1 * sut_surv_yr),
  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
  mu_g = g_int + g_slope*compsize_1 +
    g_spr_lag1 * spt1_grow_yr +
    g_spr_lag2 * spt2_grow_yr +
    g_sum_lag2 * sut_grow_yr,
  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * sp_rep_yr +
                r_sum_lag0 * sut_rep_yr +
                r_spr_lag0 * spt_rep_yr),
  data_list = ipm_parms_clim,
  states = list(c("compsize")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "t_scaled",
                                      target = "G"),
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
  
)

#fecundity kernel (Fe -> RB)
clim_stoch_ipm_1996 <- define_kernel(
  proto_ipm = clim_stoch_ipm_1996,
  name = "Fe_yr",
  family = "CD",
  formula = pf * fs * num_seeds * g_est * s_SB,
  pf = plogis(r_int + r_slope*compsize_1 +
                r_snw_lag4 * sp_rep_yr +
                r_sum_lag0 * sut_rep_yr +
                r_spr_lag0 * spt_rep_yr),
  fs = exp(f_int + f_slope*compsize_1),
  data_list = ipm_parms_clim,
  states = list(c("compsize", "RB")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)

#recruit bank kernel (RB -> P)
clim_stoch_ipm_1996 <- define_kernel(
  proto_ipm = clim_stoch_ipm_1996,
  name = "RB_yr",
  family = "DC",
  formula = s1 * f1_dist * d_compsize,
  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
  mu_f1 = rec_int_c +
    rec_slope_sm0 * sm0_rec_yr +
    rec_slope_sm1 * sm1_rec_yr,
  sd_f1 = rec_sd_c,
  data_list = ipm_parms_clim,
  states = list(c("RB", "compsize")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)



### Defining states----

#defining where the kernels go
clim_stoch_ipm_1996 <- define_impl(
  proto_ipm = clim_stoch_ipm_1996,
  make_impl_args_list(
    kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
    int_rule = c(rep("midpoint", 3)),
    state_start = c("compsize", "compsize", "RB"),
    state_end = c("compsize", "RB", "compsize")
  )
)


#defining omega
clim_stoch_ipm_1996 <- define_domains(
  proto_ipm = clim_stoch_ipm_1996,
  # 100 for number of mesh points to integrate over
  compsize = c(L, U, 100)
)

#defining population state (starting in 1996)
clim_stoch_ipm_1996 <- define_pop_state(
  proto_ipm = clim_stoch_ipm_1996,
  pop_vectors = list(
    n_compsize = initial_size_1996, 
    n_RB       = 10 #estimation
  )
)

#generating IPM
clim_stoch_ipm_1996 <- make_ipm(
  proto_ipm = clim_stoch_ipm_1996,
  iterations = 30,
  kernel_seq = 18:47,
  usr_funs = list(
    dt_scaled = dt_scaled,
    pt_scaled = pt_scaled
  ),
  return_all_envs = TRUE
)

lambda(clim_stoch_ipm_1996, "log" = FALSE)

lambda_df_1996 <- lambda_df_clim %>% 
  filter(Year >= 1996) %>% 
  mutate(crash_lambda = clim_stoch_ipm_1996$pop_state$lambda)

N_1996 = 125

absolute_pop_clim_1996<- N_1996 * c(1, cumprod(clim_stoch_ipm_1996$pop_state$lambda))

combined_pop_1996 <- combined_pop %>% 
  filter(Year >= 1996) %>% 
  mutate(crash_pop = absolute_pop_clim_1996)

df_all <- bind_rows(
  combined_pop %>% mutate(Panel = "A"),
  combined_pop_1996 %>% mutate(Panel = "B")
)

ggplot(df_all, aes(x = Year)) +
  geom_point(aes(y = log(emp_pop), color = "Empirical"), size = 2) +
  geom_line(aes(y = log(clim_pop), color = "Climate"), linewidth = 1) +
  geom_line(aes(y = log(null_pop), color = "Null"), linewidth = 1, linetype = "dashed") +
  # This line will automatically only render in Panel B where crash_pop exists
  geom_line(aes(y = log(crash_pop), color = "Climate (post 1996)"), linewidth = 1, na.rm = TRUE) +
  facet_wrap(~ Panel, scales = "free_x") +
  scale_color_manual(
    name   = "Model Type",
    values = c("Empirical"           = "black",
               "Climate"             = "#4393C3",
               "Climate (post 1996)" = "#B2182B",
               "Null"                = "#AAAAAA"),
    breaks = c("Empirical", "Climate", "Climate (post 1996)", "Null")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape    = c(16, NA, NA, NA),
        linetype = c("blank", "solid", "solid", "dashed")
      )
    )
  ) +
  labs(y = "log[N(t)]") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12, hjust = 0))

ggsave("Figure2_Population_Hindcasts.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)


ggplot(absolute_pop_clim_df, aes(x = Year, y = logN)) +
  geom_line(color = "steelblue", linewidth = 1) +
  theme_classic(base_size = 11) +
  labs(y = "log[N(t)]") +
  geom_vline(xintercept = 1996, linewidth = 0.8)

ggsave("IPM_hindcast_climate.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)






all_lambda_df <- lambda_df_clim %>% 
  mutate(emp_lambda = emp_lambda[-48])


all_lambda_1996_df <- lambda_df_1996 %>% 
  mutate(emp_lambda = emp_lambda[-c(1:17, 48)])



### Hindcast analysis----


#RMSE shows climate IPM better predicts empirical population dynamics over null IPM
#root mean square error for full data
rmse_null <- Metrics::rmse(log(combined_pop$emp_pop), log(combined_pop$null_pop))
rmse_clim <- Metrics::rmse(log(combined_pop$emp_pop), log(combined_pop$clim_pop))
print(paste("Null RMSE:", rmse_null))
print(paste("Climate RMSE:", rmse_clim))

#correlation of lambda for full data
#correlation results were not consistent with the rest of the results
#null IPM might just be steadier and that's why it has better correlation with empirical lambda?
cor_null = cor(all_lambda_df$emp_lambda, all_lambda_df$lamdba, method = "pearson")
cor_clim = cor(all_lambda_df$emp_lambda, all_lambda_df$lamdba_clim, method = "pearson")
print(paste("Climate Correlation:", cor_clim))
print(paste("Null Correlation:", cor_null))

#root mean square error for 1996 data
#RMSE shows 1996-2026 IPM better predicts declining period
#(just a shift down, all movements are the same)
rmse_null_1996 <- Metrics::rmse(log(combined_pop_1996$emp_pop), log(combined_pop_1996$null_pop))
rmse_clim_1996 <- Metrics::rmse(log(combined_pop_1996$emp_pop), log(combined_pop_1996$clim_pop))
rmse_crash_1996 <- Metrics::rmse(log(combined_pop_1996$emp_pop), log(combined_pop_1996$crash_pop))
print(paste("Null RMSE (1996-):", rmse_null_1996))
print(paste("Climate RMSE (1996-):", rmse_clim_1996))
print(paste("Climate Crash RMSE (1996-):", rmse_crash_1996))

#correlation of lambda for 1996 data
cor_null_1996 = cor(all_lambda_1996_df$emp_lambda, all_lambda_1996_df$lamdba, method = "pearson")
cor_crash_1996 = cor(all_lambda_1996_df$emp_lambda, all_lambda_1996_df$crash_lambda, method = "pearson")
cor_clim_1996 = cor(all_lambda_1996_df$emp_lambda, all_lambda_1996_df$lamdba_clim, method = "pearson")

print(paste("Climate Correlation (1996-):", cor_clim_1996))
print(paste("Null Correlation (1996-):", cor_null_1996))
print(paste("Climate Crash Correlation (1996-):", cor_crash_1996))



### Elasticity----


#null climate elasticities
run_null_ipm <- function(modified_parms) {
  
  temp_ipm <- init_ipm(sim_gen = "general", di_dd = "di", 
                       det_stoch = "stoch", kern_param = "kern") %>%
    
    # P Kernel (Survival & Growth)
    define_kernel(name = "P_yr", family = "CC",
                  formula = s * G * (1-pf) * d_compsize,
                  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
                               s_snw_lag3 * snowpack_mean +
                               s_spr_lag4 * spring_temp_mean +
                               s_sum_lag1 * summer_temp_mean),
                  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
                  mu_g = g_int + g_slope*compsize_1 +
                    g_spr_lag1 * spring_temp_mean +
                    g_spr_lag2 * spring_temp_mean +
                    g_sum_lag2 * summer_temp_mean,
                  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * snowpack_mean +
                                r_sum_lag0 * summer_temp_mean +
                                r_spr_lag0 * spring_temp_mean),
                  data_list = modified_parms,
                  states = list(c("compsize")),
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions(fun = "t_scaled", target = "G"),
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Fecundity Kernel
    define_kernel(name = "Fe_yr", family = "CD",
                  formula = pf * fs * num_seeds * g_est * s_SB,
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * snowpack_mean +
                                r_sum_lag0 * summer_temp_mean +
                                r_spr_lag0 * spring_temp_mean),
                  fs = exp(f_int + f_slope*compsize_1),
                  data_list = modified_parms,
                  states = list(c("compsize", "RB")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Recruit Bank Kernel
    define_kernel(name = "RB_yr", family = "DC",
                  formula = s1 * f1_dist * d_compsize,
                  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
                  mu_f1 = rec_int_c + rec_slope_sm0 * snowmelt_mean + rec_slope_sm1 * snowmelt_mean,
                  sd_f1 = rec_sd_c,
                  data_list = modified_parms,
                  states = list(c("RB", "compsize")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Where kernels go
    define_impl(make_impl_args_list(
      kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
      int_rule = rep("midpoint", 3),
      state_start = c("compsize", "compsize", "RB"),
      state_end = c("compsize", "RB", "compsize")
    )) %>%
    define_domains(compsize = c(L, U, 100)) %>%
    define_pop_state(pop_vectors = list(n_compsize = runif(100), n_RB = 15)) %>%
    
    # Run the model
    make_ipm(iterations = 150, # 150 is plenty for the null to reach asymptotic lambda
             kernel_seq = sample(1:47, 150, replace = TRUE),
             usr_funs = list(dt_scaled = dt_scaled, pt_scaled = pt_scaled))
  
  #return the asymptotic lambda
  return(temp_ipm$pop_state$lambda[150]) 
}

#baseline lambda to calculate elasticity
baseline_lambda <- run_null_ipm(ipm_parms)

#tested parameters
params_to_test <- c("s_int", "g_int", "r_int", "s_SB", "g_est", "num_seeds")
elasticity_results <- data.frame(Parameter = character(), 
                                 New_Lambda = numeric(), 
                                 Elasticity = numeric())

#loop to calculate elasticity
for (p in params_to_test) {
  
  #baseline parameters
  test_parms <- ipm_parms
  
  #perturb parameter by 5%
  perturbation <- test_parms[[p]] * 0.05
  test_parms[[p]] <- test_parms[[p]] + perturbation
  
  #run the model with the tweaked parameter
  new_lambda <- run_null_ipm(test_parms)
  
  #calculate Elasticity: proportional change in lambda / proportional change in parameter (0.05)
  elas <- ((new_lambda - baseline_lambda) / baseline_lambda) / 0.05
  
  #store results
  elasticity_results <- rbind(elasticity_results, 
                              data.frame(Parameter = p, 
                                         New_Lambda = new_lambda, 
                                         Elasticity = elas))
}

print(elasticity_results)

elasticity_results$Biological_Process <- c("Adult Survival", "Growth", "Probability of Flowering", 
                                           "Seedling Survival", "Germination Rate", "Number of seeds")


#1996-2025 climate resampled IPM
#does not sequentially run climate years
#resamples climate variables from modern climate regime for 200 iterations
run_clim_ipm <- function(modified_parms) {
  
  set.seed(42)
  
  temp_ipm <- init_ipm(sim_gen = "general", di_dd = "di", 
                       det_stoch = "stoch", kern_param = "kern") %>%
    
    # P Kernel (Survival & Growth)
    define_kernel(name = "P_yr", family = "CC",
                  formula = s * G * (1-pf) * d_compsize,
                  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
                               s_snw_lag3 * sp_surv_yr +
                               s_spr_lag4 * spt_surv_yr +
                               s_sum_lag1 * sut_surv_yr),
                  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
                  mu_g = g_int + g_slope*compsize_1 +
                    g_spr_lag1 * spt1_grow_yr +
                    g_spr_lag2 * spt2_grow_yr +
                    g_sum_lag2 * sut_grow_yr,
                  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * sp_rep_yr +
                                r_sum_lag0 * sut_rep_yr +
                                r_spr_lag0 * spt_rep_yr),
                  data_list = modified_parms,
                  states = list(c("compsize")),
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions(fun = "t_scaled", target = "G"),
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Fecundity Kernel
    define_kernel(name = "Fe_yr", family = "CD",
                  formula = pf * fs * num_seeds * g_est * s_SB,
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * sp_rep_yr +
                                r_sum_lag0 * sut_rep_yr +
                                r_spr_lag0 * spt_rep_yr),
                  fs = exp(f_int + f_slope*compsize_1),
                  data_list = modified_parms,
                  states = list(c("compsize", "RB")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Recruit Bank Kernel
    define_kernel(name = "RB_yr", family = "DC",
                  formula = s1 * f1_dist * d_compsize,
                  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
                  mu_f1 = rec_int_c +
                    rec_slope_sm0 * sm0_rec_yr +
                    rec_slope_sm1 * sm1_rec_yr,
                  sd_f1 = rec_sd_c,
                  data_list = modified_parms,
                  states = list(c("RB", "compsize")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Where kernels go
    define_impl(make_impl_args_list(
      kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
      int_rule = rep("midpoint", 3),
      state_start = c("compsize", "compsize", "RB"),
      state_end = c("compsize", "RB", "compsize")
    )) %>%
    define_domains(compsize = c(L, U, 100)) %>%
    define_pop_state(pop_vectors = list(n_compsize = initial_size, n_RB = 10)) %>%
    
    # Run the model
    make_ipm(iterations = 200,
             kernel_seq = sample(18:47, size = 200, replace = TRUE),
             usr_funs = list(dt_scaled = dt_scaled,
                             pt_scaled = pt_scaled))
  
  #return the asymptotic lambda
  annual_lambdas <- temp_ipm$pop_state$lambda
  stoch_lambda <- exp(mean(log(annual_lambdas)))
  
  return(stoch_lambda)
}

baseline_lambda_climate <- run_clim_ipm(ipm_parms_clim)

#tested parameters
params_to_test <- c("s_int", "g_int", "r_int", "s_SB", "g_est", "num_seeds")
elasticity_results_climate <- data.frame(Parameter = character(), 
                                 New_Lambda = numeric(), 
                                 Elasticity = numeric())


#loop to calculate elasticity
for (p in params_to_test) {
  
  #baseline parameters
  test_parms <- ipm_parms_clim
  
  #perturb parameter by 5%
  perturbation <- test_parms[[p]] * 0.05
  test_parms[[p]] <- test_parms[[p]] + perturbation
  
  #run the model with the tweaked parameter
  new_lambda <- run_clim_ipm(test_parms)
  
  #calculate Elasticity: proportional change in lambda / proportional change in parameter (0.05)
  elas <- ((new_lambda - baseline_lambda_climate) / baseline_lambda_climate) / 0.05
  
  #store results
  elasticity_results_climate <- rbind(elasticity_results_climate, 
                              data.frame(Parameter = p, 
                                         New_Lambda = new_lambda, 
                                         Elasticity = elas))
}


elasticity_results_climate$Biological_Process <- c("Adult Survival", "Growth", "Probability of Flowering", 
                                           "Seedling Survival", "Germination Rate", "Number of seeds")




#1979-1995 climate regime resampled IPM
run_clim_pre_ipm <- function(modified_parms) {
  
  set.seed(42)
  
  temp_ipm <- init_ipm(sim_gen = "general", di_dd = "di", 
                       det_stoch = "stoch", kern_param = "kern") %>%
    
    # P Kernel (Survival & Growth)
    define_kernel(name = "P_yr", family = "CC",
                  formula = s * G * (1-pf) * d_compsize,
                  s = plogis(s_int + s_slope*compsize_1 + s_quad*(compsize_1)^2 +
                               s_snw_lag3 * sp_surv_yr +
                               s_spr_lag4 * spt_surv_yr +
                               s_sum_lag1 * sut_surv_yr),
                  G = dt_scaled(compsize_2, mean = mu_g, sd = g_sd, df = t_df),
                  mu_g = g_int + g_slope*compsize_1 +
                    g_spr_lag1 * spt1_grow_yr +
                    g_spr_lag2 * spt2_grow_yr +
                    g_sum_lag2 * sut_grow_yr,
                  g_sd = exp(0.5 * (d_int + d_slope*compsize_1)),
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * sp_rep_yr +
                                r_sum_lag0 * sut_rep_yr +
                                r_spr_lag0 * spt_rep_yr),
                  data_list = modified_parms,
                  states = list(c("compsize")),
                  evict_cor = TRUE,
                  evict_fun = truncated_distributions(fun = "t_scaled", target = "G"),
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Fecundity Kernel
    define_kernel(name = "Fe_yr", family = "CD",
                  formula = pf * fs * num_seeds * g_est * s_SB,
                  pf = plogis(r_int + r_slope*compsize_1 +
                                r_snw_lag4 * sp_rep_yr +
                                r_sum_lag0 * sut_rep_yr +
                                r_spr_lag0 * spt_rep_yr),
                  fs = exp(f_int + f_slope*compsize_1),
                  data_list = modified_parms,
                  states = list(c("compsize", "RB")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Recruit Bank Kernel
    define_kernel(name = "RB_yr", family = "DC",
                  formula = s1 * f1_dist * d_compsize,
                  f1_dist = dnorm(compsize_2, mu_f1, sd_f1),
                  mu_f1 = rec_int_c +
                    rec_slope_sm0 * sm0_rec_yr +
                    rec_slope_sm1 * sm1_rec_yr,
                  sd_f1 = rec_sd_c,
                  data_list = modified_parms,
                  states = list(c("RB", "compsize")),
                  evict_cor = FALSE,
                  uses_par_sets = TRUE,
                  par_set_indices = list(yr = 1:47)) %>%
    
    # Where kernels go
    define_impl(make_impl_args_list(
      kernel_names = c("P_yr", "Fe_yr", "RB_yr"),
      int_rule = rep("midpoint", 3),
      state_start = c("compsize", "compsize", "RB"),
      state_end = c("compsize", "RB", "compsize")
    )) %>%
    define_domains(compsize = c(L, U, 100)) %>%
    define_pop_state(pop_vectors = list(n_compsize = initial_size, n_RB = 10)) %>%
    
    # Run the model
    make_ipm(iterations = 200,
             kernel_seq = sample(1:17, size = 200, replace = TRUE),
             usr_funs = list(dt_scaled = dt_scaled,
                             pt_scaled = pt_scaled))
  
  #return the asymptotic lambda
  annual_lambdas <- temp_ipm$pop_state$lambda
  stoch_lambda <- exp(mean(log(annual_lambdas)))
  
  return(stoch_lambda)
}


baseline_lambda_pre_climate <- run_clim_pre_ipm(ipm_parms_clim)

#tested parameters
params_to_test <- c("s_int", "g_int", "r_int", "s_SB", "g_est", "num_seeds")
elasticity_results_pre_climate <- data.frame(Parameter = character(), 
                                         New_Lambda = numeric(), 
                                         Elasticity = numeric())


#loop to calculate elasticity
for (p in params_to_test) {
  
  #baseline parameters
  test_parms <- ipm_parms_clim
  
  #perturb parameter by 5%
  perturbation <- test_parms[[p]] * 0.05
  test_parms[[p]] <- test_parms[[p]] + perturbation
  
  #run the model with the tweaked parameter
  new_lambda <- run_clim_pre_ipm(test_parms)
  
  #calculate Elasticity: proportional change in lambda / proportional change in parameter (0.05)
  elas <- ((new_lambda - baseline_lambda_pre_climate) / baseline_lambda_pre_climate) / 0.05
  
  #store results
  elasticity_results_pre_climate <- rbind(elasticity_results_pre_climate, 
                                      data.frame(Parameter = p, 
                                                 New_Lambda = new_lambda, 
                                                 Elasticity = elas))
}


elasticity_results_pre_climate$Biological_Process <- c("Adult Survival", "Growth", "Probability of Flowering", 
                                                   "Seedling Survival", "Germination Rate", "Number of seeds")
df_null <- elasticity_results %>%
  mutate(Model = "Mean Climate (Null)")

df_clim <- elasticity_results_climate %>%
  mutate(Model = "Stochastic Climate (1996-2025)")

df_pre_clim <- elasticity_results_pre_climate %>% 
  mutate(Model = "Stochastic Climate (1979-1995)")

combined_elasticity <- bind_rows(df_null, df_clim, df_pre_clim)
combined_elasticity$Model <- factor(combined_elasticity$Model, 
                                    levels = c("Mean Climate (Null)",
                                               "Stochastic Climate (1996-2025)",
                                               "Stochastic Climate (1979-1995)"))

climate_elasticities <- bind_rows(df_clim, df_pre_clim)


#elasticity plot
ggplot(combined_elasticity, aes(x = reorder(Biological_Process, Elasticity), y = Elasticity, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(
    aes(
      label = sprintf("%.3f", Elasticity),
      hjust = ifelse(Elasticity >= 0, -0.2, 1.2)
    ),
    position = position_dodge(width = 0.8),
    size = 3.2
  ) +
  
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  scale_fill_manual(values = c("Mean Climate (Null)" = "gray70", 
                               "Stochastic Climate (1996-2025)" = "darkred",
                               "Stochastic Climate (1979-1995)" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  theme_classic(base_size = 15) +
  labs(x = "Biological Process",
       y = "Elasticity",
       fill = "Model Type") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )

ggsave("Elasticity_Figure_stochastic.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)
  
ggsave("Survival_Elasticites.png", 
       plot = last_plot(), 
       width = 8.5, 
       height = 5.5, 
       units = "in", 
       dpi = 300)





### Monocarpy analysis----

#null df
r_int_surface <- data.frame()
baseline_r_int  <- ipm_parms_clim$r_int

#additive perturbation of reproduction intercept
#baseline r_int is -5.751604
#creates a linear relationship (for this domain of perturbation)
#modern climate (1996-2025)
for (delta in seq(from = -1.5, to = 1.5, by = 0.15)) {
  
  test_parms <- ipm_parms_clim
  
  test_parms$r_int <- baseline_r_int + delta
  
  new_lambda <- run_clim_ipm(test_parms)
  
  sens <- (new_lambda - baseline_lambda_climate) / delta
  if (delta == 0) sens <- 0
  
  r_int_surface <- rbind(r_int_surface, data.frame(
    Delta_Shift = delta,               # Negative = lower repro, Positive = higher repro
    r_int_value = test_parms$r_int,   # Absolute intercept value
    New_Lambda  = new_lambda,
    Sensitivity = sens
  ))
  
  

}

print(r_int_surface)

ggplot(r_int_surface, aes(x = Delta_Shift, y = New_Lambda)) +
  geom_point()

#null df
r_int_surface_pre <- data.frame()
#historical climate (1979-2025)
for (delta in seq(from = -1.5, to = 1.5, by = 0.15)) {
  
  test_parms <- ipm_parms_clim
  
  test_parms$r_int <- baseline_r_int + delta
  
  new_lambda <- run_clim_pre_ipm(test_parms)
  
  sens <- (new_lambda - baseline_lambda_pre_climate) / delta
  if (delta == 0) sens <- 0
  
  r_int_surface_pre <- rbind(r_int_surface_pre, data.frame(
    Delta_Shift = delta,               # Negative = lower repro, Positive = higher repro
    r_int_value = test_parms$r_int,   # Absolute intercept value
    New_Lambda  = new_lambda,
    Sensitivity = sens
  ))
}

ggplot(r_int_surface_pre, aes(x = Delta_Shift, y = New_Lambda)) +
  geom_point()

r_int_surface <- r_int_surface %>% 
  mutate(Model = as.factor("Modern Climate"))

r_int_surface_pre <- r_int_surface_pre %>% 
  mutate(Model = as.factor("Historical Climate"))

r_surface <- bind_rows(r_int_surface, r_int_surface_pre)

ggplot(r_surface, aes(x = Delta_Shift, y = New_Lambda, color = Model)) +
  geom_point() +
  scale_color_manual(values = c(
    "Modern Climate" = "red",
    "Historical Climate" = "green"
  )) +
  theme_classic(base_size = 11) +
  labs(y = "New Lambda", x = "Delta Shift")


#slopes differ between historical climate and modern climate
#obviously modern climate is shifted down, but the slopes shift too

mod_lm <- lm(r_int_surface$New_Lambda ~ r_int_surface$Delta_Shift)
summary(mod_lm)

his_lm <- lm(r_int_surface_pre$New_Lambda ~ r_int_surface_pre$Delta_Shift)
summary(his_lm)

r_slope_surface <- data.frame()
baseline_r_slope <- ipm_parms_clim$r_slope

for (delta in seq(from = -1.5, to = 1.5, by = 0.1)) {
  
  test_parms <- ipm_parms_clim
  test_parms$r_slope <- baseline_r_slope + delta
  
  new_lambda    <- run_clim_ipm(test_parms)
  new_lambda_pre <- run_clim_pre_ipm(test_parms)
  
  sens_mod <- (new_lambda - baseline_lambda_climate) / delta
  sens_his <- (new_lambda_pre - baseline_lambda_pre_climate) / delta
  if (delta == 0) sens_mod <- 0
  if (delta == 0) sens_his <- 0
  
  r_slope_surface <- rbind(r_slope_surface, data.frame(
    Delta_Shift   = delta,
    r_slope_value = test_parms$r_slope,
    New_Lambda    = new_lambda,
    New_Lambda_Pre = new_lambda_pre,
    Sensitivity_Mod   = sens_mod,
    Sensitivity_his = sens_his
  ))
}

#nonlinear so ignore these
slope_lm_mod  <- lm(New_Lambda     ~ Delta_Shift, data = r_slope_surface)
slope_lm_pre  <- lm(New_Lambda_Pre ~ Delta_Shift, data = r_slope_surface)
summary(slope_lm_mod)
summary(slope_lm_pre)



ggplot(r_slope_surface, aes(x = Delta_Shift)) +
  geom_point(aes(y = New_Lambda), color = "green") +
  geom_point(aes(y = New_Lambda_Pre), color = "red") +
  theme_classic(base_size = 11) +
  labs(y = "Lambda", x = "Delta Shift")

ggplot(r_slope_surface, aes(x = Delta_Shift)) +
  geom_point(aes(y = Sensitivity_Mod), color = "green") +
  geom_point(aes(y = Sensitivity_his), color = "red") +
  labs(y = "Sensitivity", x = "Delta Shift") +
  theme_classic(base_size = 11)

#lambda gams
#slopes differ under climate regimes
non_mod <- gam(New_Lambda ~ s(Delta_Shift),
                   family = gaussian,
                   data = r_slope_surface)
summary(non_mod)
non_his <- gam(New_Lambda_Pre ~ s(Delta_Shift),
        family = gaussian,
        data = r_slope_surface)
summary(non_his)

#sensitivity gams
#slopes differ under climate regimes
non_sens_mod <- gam(Sensitivity_Mod ~ s(Delta_Shift),
               family = gaussian,
               data = r_slope_surface)
summary(non_sens_mod)

non_sens_his <- gam(Sensitivity_his ~ s(Delta_Shift),
                    family = gaussian,
                    data = r_slope_surface)
summary(non_sens_his)


#long-format data for analysis
long_lambda <- pivot_longer(
  data = r_slope_surface,
  cols = c(New_Lambda, New_Lambda_Pre),
  names_to = "Model",
  values_to = "Lambda"
) %>% 
  mutate(Model = factor(Model))


#better sensitivity analysis
combined_sens <- gam(Lambda ~ Model + s(Delta_Shift, by = Model),
                     data = long_lambda)

#modern climate derivative curve (under perturbation set)
d_modern <- gratia::derivatives(combined_sens, select = "s(Delta_Shift):ModelNew_Lambda",
                        data = data.frame(Delta_Shift = seq(-1.5, 1.5, length.out = 300),
                                          Model = "New_Lambda"))
#historical climate derivative curve
d_hist   <- gratia::derivatives(combined_sens, select = "s(Delta_Shift):ModelNew_Lambda_Pre",
                        data = data.frame(Delta_Shift = seq(-1.5, 1.5, length.out = 300),
                                          Model = "New_Lambda_Pre"))

deriv_diff <- data.frame(
  Delta_Shift = d_modern$Delta_Shift,
  diff = d_modern$.derivative - d_hist$.derivative,
  se   = sqrt(d_modern$.se^2 + d_hist$.se^2)  # valid if the two group smooths are fit independently/orthogonally, which they are with by=Model
) %>%
  mutate(lower = diff - 1.96*se, upper = diff + 1.96*se)

#different additive perturbations yield different derivatives under
#historical and modern climates
ggplot(deriv_diff, aes(Delta_Shift, diff)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed")



