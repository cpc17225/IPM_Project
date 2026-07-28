### Packages and loading data----

library(ipmr)
library(tidyverse)

ipm_parms <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_comp.RDS")
initial_size <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/initial_size_vector.RDS")


### Helper functions/variables----
dt_scaled <- function(x, mean, sd, df) {
  dt((x - mean) / sd, df = df) / sd
}

pt_scaled <- function(q, mean, sd, df) {
  pt((q - mean) / sd, df = df)
}

#omega bounds for size variables
#lower bound for size (-0.5 to account for size eviction)
L = -4.011008 - 0.5
#upper bound for size (+0.5 to account for size eviction)
U = 6.060386 + 0.5

### Building kernels----

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
    int_rule = c(rep("midpoint", 3)),
    state_start = c("compsize", "compsize", "RB"),
    state_end = c("compsize", "RB", "compsize")
  )
)


#defining omega
null_stoch_ipm <- define_domains(
  proto_ipm = null_stoch_ipm,
  # 100 for number of meshpoints to integrate over
  compsize = c(L, U, 100)
)

#defining initial population state
null_stoch_ipm <- define_pop_state(
  proto_ipm = null_stoch_ipm,
  pop_vectors = list(
    n_compsize = initial_size,
    #recruit bank number is rough estimation
    n_RB = 5
  )
)

#generating IPM
null_stoch_ipm <- make_ipm(
  proto_ipm = null_stoch_ipm,
  iterations = 47,
  kernel_seq = 1:47,
  usr_funs = list(
    dt_scaled = dt_scaled,
    pt_scaled = pt_scaled
  ),
  return_all_envs = TRUE
)

lambda(null_stoch_ipm, log = FALSE)
is_conv_to_asymptotic(null_stoch_ipm)

annual_lambdas <- null_stoch_ipm$pop_state$lambda
lambda_df <- data.frame(
  Year = 1979:2025,
  lamdba = annual_lambdas,
  log_lambda = log(annual_lambdas)
)

lambda_df

