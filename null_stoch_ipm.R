### Packages and loading data----

library(ipmr)
library(tidyverse)

ipm_parms <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_comp.RDS")


### Building kernels

null_ipm <- init_ipm(sim_gen = "general",
                 di_dd = "di",
                 det_stoch = "stoch",
                 kern_param = "kern")

#growth/survival kernel (P -> P)
null_ipm <- define_kernel()