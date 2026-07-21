### Packages and loading data----

library(tidyverse)
library(ipmr)
library(popdemo)

ipm_parms <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms.rds")

### Creating IPM Kernels----

null_ipm <- init_ipm(sim_gen = "general",
                     di_dd = "di",
                     det_stoch = "det")

#Growth and survival kernel (P -> P)
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "P",
  family = "CC",
  formula = s * G * (1-pf) * d_logLLL,
  s = plogis(s_int + s_slope*logLLL_1 + s_quad*(logLLL_1)^2),
  G = dnorm(logLLL_2, mu_g, sqrt(g_var)),
  mu_g = g_int + g_slope * logLLL_1,
  g_var = exp(d_int + d_slope * logLLL_1),
  pf = plogis(r_int + r_slope*logLLL_1),
  data_list = ipm_parms,
  states = list(c("logLLL")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
  
)

#Fecundity kernel (F -> SB)
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "Fe",
  family = "CD",
  formula = pf * fs * num_seeds * g_est * s_SB,
  pf = plogis(r_int + r_slope*logLLL_1),
  fs = exp(f_int + f_slope*logLLL_1),
  data_list = ipm_parms,
  states = list(c("logLLL", "SB")),
  evict_cor = FALSE
)

#Seedling Kernel (SB -> RB1)
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "S",
  family = "DD",
  formula = s1,
  data_list = ipm_parms,
  states = list(c("SB", "RB1")),
  evict_cor = FALSE
)

#Recruit Bank 1 to Recruit Bank 2 Kernel (RB1 -> RB2)
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "T_RB1_RB2",
  family = "DD",
  formula = s2 * (1-det_p),
  data_list = ipm_parms,
  states = list(c("RB1", "RB2")),
  evict_cor = FALSE
)

#Recruit Bank 1 to Population (RB1 -> P)
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "R1",
  family = "DC",
  formula = s2 * det_p * f1_dist * d_logLLL,
  f1_dist = dnorm(logLLL_2, f1_mean, f1_sd),
  states = list(c("RB1", "logLLL")),
  data_list = ipm_parms,
  evict_cor = FALSE
)

#Recruit bank 2 to Populations
null_ipm <- define_kernel(
  proto_ipm = null_ipm,
  name = "R2",
  family = "DC",
  formula = det_p * f2_dist * s2 * d_logLLL,
  f2_dist = dnorm(logLLL_2, f2_mean, f2_sd),
  states = list(c("RB2", "logLLL")),
  data_list = ipm_parms,
  evict_cor = FALSE
)

#Defining where the kernels go
null_ipm <- define_impl(
  proto_ipm = null_ipm,
  make_impl_args_list(
    kernel_names = c("P", "Fe", "S", "T_RB1_RB2", "R1", "R2"),
    int_rule = c(rep("midpoint", 6)),
    state_start = c("logLLL", "logLLL", "SB", "RB1", "RB1", "RB2"),
    state_end = c("logLLL", "SB", "RB1", "RB2", "logLLL", "logLLL")
  )
)

### States and initialization----

#Defining omega
null_ipm <- define_domains(
  proto_ipm = null_ipm,
  logLLL = c(0.65,
             5.6,
             100)
)

#Defining initial population state given average population sizes and estimations
# n_logLLL follows a random uniform distribution for starting population vector
# 100 is number of mesh points
# The other constants are abritrary - do not affect lambda
null_ipm <- define_pop_state(
  proto_ipm = null_ipm,
  pop_vectors = list(
    n_logLLL = runif(100),
    n_SB = 25,
    n_RB1 = 20,
    n_RB2 = 10
  )
)

#Generating IPM
null_ipm <- make_ipm(proto_ipm = null_ipm,
                     iterations = 100)

### Diagnostics----
lambda(null_ipm)

#Converges = TRUE
is_conv_to_asymptotic(null_ipm)

#Calculating left and right eigenvectors
right <- right_ev(null_ipm)
left <- left_ev(null_ipm)
#Dot product of left EV and right EV
dot <- sum(right$logLLL_w * left$logLLL_v)
#Sensitivity of size (given in mesh rank) to lambda
sensitivity <- (right$logLLL_w * left$logLLL_v)/dot
plot(sensitivity)

# Make this explicit - overlay sensitivity with size distribution
par(mfrow = c(1,1))

# Convert mesh index to actual log(LLL) values
mesh_points <- seq(0.65, 5.6, length.out = length(sensitivity))

plot(mesh_points, sensitivity,
     xlab = "log(LLL)", ylab = "Sensitivity",
     main = "Size-specific demographic importance",
     type = "l", lwd = 2)

# Also plot the stable size distribution alone
plot(mesh_points, right$logLLL_w / sum(right$logLLL_w),
     xlab = "log(LLL)", ylab = "Proportion",
     main = "Stable size distribution",
     type = "l", lwd = 2)

# Full sensitivity matrix (outer product)
n <- length(right$logLLL_w)
S_matrix <- outer(left$logLLL_v, right$logLLL_w) / dot

K_mat <- matrix(0, nrow = 103, ncol = 103)

# 2. Extract your list of 6 computed kernels from the model
kernels <- null_ipm$sub_kernels

# 3. Sew the quilt together (matching your EQ1 - EQ5)

# N -> N (EQ 5: Adult survival & growth)
K_mat[1:100, 1:100] <- kernels$P

# N -> SB (EQ 1: Adults producing seeds into the bank)
K_mat[101, 1:100] <- kernels$Fe

# SB -> RB1 (EQ 2: Seeds surviving into first-year recruits)
K_mat[102, 101] <- kernels$S

# RB1 -> RB2 (EQ 3: First-year recruits surviving to second year)
K_mat[103, 102] <- kernels$T_RB1_RB2

# RB1 -> N (EQ 4 part 1: First-year recruits emerging into continuous sizes)
K_mat[1:100, 102] <- kernels$R1

# RB2 -> N (EQ 4 part 2: Second-year recruits emerging into continuous sizes)
K_mat[1:100, 103] <- kernels$R2
w_full <- c(right$logLLL_w, right$SB_w, right$RB1_w, right$RB2_w)
v_full <- c(left$logLLL_v, left$SB_v, left$RB1_v, left$RB2_v)

# 2. Calculate the global dot product (the inner product)
dot_full <- sum(v_full * w_full)

# 3. Calculate the 103x103 Sensitivity Matrix via the outer product
S_matrix <- outer(v_full, w_full) / dot_full

# 4. Calculate the Elasticity Matrix (Hadamard/element-wise product)
#    Note: Make sure your lambda is a single numeric scalar here
E_matrix <- (S_matrix * K_mat) / Re(eigen(K_mat)$values[1])

# Sum of elasticities should = 1
cat("Sum of elasticities:", sum(E_matrix), "\n")

# Plot as surface
mesh_points <- seq(0.65, 5.6, length.out = n)
E_matrix_cont <- E_matrix[-c(101,102,103),-c(101,102,103)]
image(mesh_points, mesh_points, t(E_matrix_cont),
      xlab = "Size at t (log LLL)",
      ylab = "Size at t+1 (log LLL)",
      main = "Elasticity surface K(z',z)")

# Marginal elasticity - collapse over future size
# Tells you which CURRENT sizes drive lambda
elas_marginal <- colSums(E_matrix_cont) * (5.6 - 0.65) / n
plot(mesh_points, elas_marginal,
     type = "l", lwd = 2,
     xlab = "log(LLL)",
     ylab = "Elasticity",
     main = "Marginal elasticity by current size")

E_matrix[101,101]



### ipmr Analyses----

#Examining average lifespan for a given size
#Code from ipmr vignette
#Confused how we can get average lifespan?
#Also seems really low?
make_N <- function(ipm) {
  P <- ipm$sub_kernel$P
  I <- diag(nrow(P))
  N <- solve(I-P)
  return(N)
}

eta_bar <- function(ipm){
  N <- make_N(ipm)
  out = colSums(N)
  return(as.vector(out))
  
}

sigma_eta <- function(ipm){
  N <- make_N(ipm)
  out <- colSums(2*(N %^% 2L)-N) - colSums(N)^2
  return(as.vector(out))
}

mean_1 <- eta_bar(null_ipm)
var_1 <- sigma_eta(null_ipm)
mesh_ps <- int_mesh(null_ipm)$logLLL_1 %>%
  unique()
plot(mesh_ps, mean_1, type = "l", xlab = expression( "Initial size z"[0]))
plot(mesh_ps, var_1, type = "l", xlab = expression( "Initial size z"[0]))


#From ipmr vignette
#Does not work, not sure why
#Trying to make mega-kernel for plotting and further analysis
sub_k_nms <- names(null_ipm$sub_kernels)
mega_mat_text <- c(sub_k_nms[1], sub_k_nms[2], sub_k_nms[3],
                   sub_k_nms[4], sub_k_nms[5], sub_k_nms[6])
mega_mat_2 <- make_iter_kernel(null_ipm,
                               mega_mat = mega_mat_text)