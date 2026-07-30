library(tidyverse)

ipm_parms_base <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/ipmr_parms_comp.RDS")

### Lagged climate dataframe----
Climate_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Climate_data_temp.rds")
Climate_data <- Climate_data %>% 
  rename("Year" = year) %>% 
  mutate(scaled_snowmelt = scale(snowmelt)) %>% 
  mutate(scaled_snowpack = scale(snowpack)) %>% 
  mutate(scaled_summer_temp = scale(summer.mean.temp)) %>% 
  mutate(scaled_spring_temp = scale(spring.mean.temp))
climate_lagged <- Climate_data %>% 
  arrange(Year)
for (i in 0:4) {
  climate_lagged[[paste0("snowmelt_lag", i)]] <- dplyr::lag(climate_lagged$scaled_snowmelt, i)
  climate_lagged[[paste0("snowpack_lag", i)]] <- dplyr::lag(climate_lagged$scaled_snowpack, i)
  climate_lagged[[paste0("summer_temp_lag", i)]] <- dplyr::lag(climate_lagged$scaled_summer_temp, i)
  climate_lagged[[paste0("spring_temp_lag", i)]] <- dplyr::lag(climate_lagged$scaled_spring_temp, i)
}

climate_lagged <- climate_lagged %>% 
  filter(Year >= 1979 & Year < 2026)


### Survival climate variables----

#snowpack lag 3
sp_surv_list <- as.list(climate_lagged$snowpack_lag3)
names(sp_surv_list) <- paste0("sp_surv_", 1:47)

#spring temp lag 4
spt_surv_list <- as.list(climate_lagged$spring_temp_lag4)
names(spt_surv_list) <- paste0("spt_surv_", 1:47)

#summer temp lag 1
sut_surv_list <- as.list(climate_lagged$summer_temp_lag1)
names(sut_surv_list) <- paste0("sut_surv_", 1:47)



### Reproduction climate variables----

#snowpack lag 4
sp_rep_list <- as.list(climate_lagged$snowpack_lag4)
names(sp_rep_list) <- paste0("sp_rep_", 1:47)

#spring temp lag 0
spt_rep_list <- as.list(climate_lagged$spring_temp_lag0)
names(spt_rep_list) <- paste0("spt_rep_", 1:47)

#summer temp lag 0
sut_rep_list <- as.list(climate_lagged$summer_temp_lag0)
names(sut_rep_list) <- paste0("sut_rep_", 1:47)



### Growth climate variables----

#spring temp lag 1
spt1_grow_list <- as.list(climate_lagged$spring_temp_lag1)
names(spt1_grow_list) <- paste0("spt1_grow_", 1:47)

#spring temp lag 1
spt2_grow_list <- as.list(climate_lagged$spring_temp_lag2)
names(spt2_grow_list) <- paste0("spt2_grow_", 1:47)

#summer temp lag 2
sut_grow_list <- as.list(climate_lagged$summer_temp_lag2)
names(sut_grow_list) <- paste0("sut_grow_", 1:47)



### Recruit climate variables

#snowmelt lag 0
sm0_rec_list <- as.list(climate_lagged$snowmelt_lag0)
names(sm0_rec_list) <- paste0("sm0_rec_", 1:47)

#snowmelt lag 1
sm1_rec_list <- as.list(climate_lagged$snowmelt_lag1)
names(sm1_rec_list) <- paste0("sm1_rec_", 1:47)


ipm_parms_climate <- c(ipm_parms_base,
                       sp_surv_list,
                       spt_surv_list,
                       sut_surv_list,
                       sp_rep_list,
                       spt_rep_list,
                       sut_rep_list,
                       spt1_grow_list,
                       spt2_grow_list,
                       sut_grow_list,
                       sm0_rec_list,
                       sm1_rec_list)

saveRDS(ipm_parms_climate, file = "ipm_parms_climate.rds")                       
