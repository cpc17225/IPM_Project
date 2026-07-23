### Packages and loading data----


library(ipmr)
library(tidyverse)

## Adding climate data to IPM parameters

ipm_parms_base <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_c.rds")
climate_data_lag <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/climate_data_with_lag.rds")

sp_surv_list <- as.list(climate_data_lag$snowpack_lag3)
names(sp_surv_list) <- paste0("sp_surv_", 1:47)

sp_rep_list <- as.list(climate_data_lag$snowpack_lag4)
names(sp_rep_list) <- paste0("sp_rep_", 1:47)

# All combined into one list
ipm_parms_c <- c(ipm_parms_base, sp_surv_list, sp_rep_list)


### Creating IPM Kernels----

climate_ipm <- init_ipm(sim_gen = "general",
                        di_dd = "di",
                        det_stoch = "stoch",
                        kern_param = "kern")
#Defining growth kernel (P -> P)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "P_yr",
  family = "CC",
  formula = s * G * (1-pf) * d_logLLL,
  s = plogis(s_int_c + s_slope_c*logLLL_1 + s_quad_c*(logLLL_1)^2 +
               s_slope_sp*sp_surv_yr),
  G = dnorm(logLLL_2, mu_g, sqrt(g_var)),
  g_var = exp(d_int_c + d_slope_c*logLLL_1),
  mu_g = g_int_c + g_slope_c * logLLL_1,
  pf = plogis(r_int_c + r_slope_c*logLLL_1 + 
                r_slope_sp * sp_rep_yr),
  data_list = ipm_parms_c,
  states = list(c("logLLL")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G"),
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)

#Defining fecundity kernel (Fe -> SB)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "Fe_yr",
  family ="CD",
  formula = pf * fs * num_seeds * g_est * s_SB,
  pf = plogis(r_int_c + r_slope_c * logLLL_1 +
                r_slope_sp * sp_rep_yr),
  fs = exp(f_int_c + f_slope_c * logLLL_1),
  data_list = ipm_parms_c,
  states = list(c("logLLL", "SB")),
  evict_cor = FALSE,
  uses_par_sets = TRUE,
  par_set_indices = list(yr = 1:47)
)

#Defining seedling kernel (SB -> RB1)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "S",
  family = "DD",
  formula = s1,
  data_list = ipm_parms_c,
  states = list(c("SB", "RB1")),
  evict_cor = FALSE,
  uses_par_sets = FALSE
)

#Defining Recruit Bank 1 to Recruit Bank 2 kernel (RB1 -> RB2)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "T_RB1_RB2",
  family = "DD",
  formula = s2 * (1-det_p),
  data_list = ipm_parms_c,
  states = list(c("RB1", "RB2")),
  uses_par_sets = FALSE,
  evict_cor = FALSE
)

#Defining Recruit Bank 1 to Population (RB1 -> P)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "R1",
  family = "DC",
  formula = s2 * det_p * f1_dist * d_logLLL,
  f1_dist = dnorm(logLLL_2, f1_mean, f1_sd),
  states = list(c("RB1", "logLLL")),
  data_list = ipm_parms_c,
  evict_cor = FALSE,
  uses_par_sets = FALSE
)

#Defining Recruit Bank 2 to Population (RB2 -> P)
climate_ipm <- define_kernel(
  proto_ipm = climate_ipm,
  name = "R2",
  family = "DC",
  formula = det_p * f2_dist * s2 * d_logLLL,
  f2_dist = dnorm(logLLL_2, f2_mean, f2_sd),
  states = list(c("RB2", "logLLL")),
  data_list = ipm_parms_c,
  evict_cor = FALSE,
  uses_par_sets = FALSE
)

#Defining where the kernels go
climate_ipm <- define_impl(
  proto_ipm = climate_ipm,
  make_impl_args_list(
    kernel_names = c("P_yr", "Fe_yr", "S", "T_RB1_RB2", "R1", "R2"),
    int_rule = c(rep("midpoint", 6)),
    state_start = c("logLLL", "logLLL", "SB", "RB1", "RB1", "RB2"),
    state_end = c("logLLL", "SB", "RB1", "RB2", "logLLL", "logLLL")
  )
)

#Defining omega
climate_ipm <- define_domains(
  proto_ipm = climate_ipm,
  logLLL = c(0.65,
             5.6,
             100)
)

#Defining initial population state given average population sizes and estimations
climate_ipm <- define_pop_state(
  proto_ipm = climate_ipm,
  pop_vectors = list(
    n_logLLL = runif(100),
    n_SB = 25,
    n_RB1 = 20,
    n_RB2 = 10
  )
)

#Generating IPM
climate_ipm <- make_ipm(
  proto_ipm = climate_ipm,
  iterations = 47,
  kernel_seq = 1:47,
  usr_funs = list(dt = dt),
  return_all_envs = TRUE
)


### Diagnostics and data----


lambda(climate_ipm, log = FALSE)
is_conv_to_asymptotic(climate_ipm)

annual_lambdas <- climate_ipm$pop_state$lambda
lambda_df <- data.frame(
  Year = 1979:2025,
  lambda = annual_lambdas,
  log_lambda = log(annual_lambdas)
)
