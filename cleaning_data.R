### Packages and loading data----

library(tidyverse)
library(readxl)
library(readr)

Grandiflora_Cleaned <- read_excel("Grandiflora_Cleaned.xlsx", na = "NA")
Climate_Data <- read_csv("C:/Users/Owner/Downloads/barr_snowmelt_date_2022.csv")

### Renaming, data structures----

#Main dataset, name changes and conditional columns
clean_data <- Grandiflora_Cleaned %>% 
  rename("survival" = "Alive(t+1)") %>% 
  rename("numleaves" = "NumLeaves(t)") %>% 
  rename("numleavesprev" = "NumLeaves(t-1)") %>% 
  rename("LLL" = "LLL(t)") %>% 
  rename("LLLprev" = "LLL(t-1)") %>% 
  mutate(Gophered = if_else(Comments == "Gophered", 1, 0, missing = 0)) %>% 
  mutate(rep = if_else(NumFlowers>0, 1, 0)) %>% 
  mutate(logLLL = log(LLL))

str(clean_data)

#Data structure changes
clean_data$Year <- as.numeric(clean_data$Year)
clean_data$numleavesprev <- as.numeric(clean_data$numleavesprev)
clean_data$LLLprev <- as.numeric(clean_data$LLLprev)

### Removing outliers----

#Flower on side of rosette. Plant survived. Removing flower
clean_data[4408, "NumFlowers"] <- 0

#Plant survived after reproductive event. Not possible (only survived for one year)
clean_data <- clean_data %>% slice(-3189)

#LLL and numleaves should always be NA when a plant is flowering
clean_data <- clean_data %>% 
  mutate(LLL = replace(LLL, rep == 1, NA)) %>% 
  mutate(numleaves = replace(numleaves, rep == 1, NA))

#Cases when numleaves is too high for a low LLL
#Could be from herbivory or bad measurement
clean_data[3285, "numleaves"] <- NA
clean_data[3285, "LLL"] <- NA
clean_data[3286, "numleavesprev"] <- NA
clean_data[3287, "LLLprev"] <- NA
clean_data[3131, "numleaves"] <- NA
clean_data[3131, "LLL"] <- NA
clean_data[3616, "numleaves"] <- NA
clean_data[3616, "LLL"] <- NA
clean_data[1015, "numleaves"] <- NA
clean_data[1015, "LLL"] <- NA
clean_data[1016, "numleavesprev"] <- NA
clean_data[1016, "LLLprev"] <- NA
clean_data[1484, "numleaves"] <- NA
clean_data[1484, "LLL"] <- NA
clean_data[1485, "numleavesprev"] <- NA
clean_data[1485, "LLLprev"] <- NA
clean_data[2244, "numleaves"] <- NA
clean_data[2244, "LLL"] <- NA
clean_data[2245, "numleavesprev"] <- NA
clean_data[2245, "LLLprev"] <- NA
clean_data[2931, "numleaves"] <- NA
clean_data[2931, "LLL"] <- NA
clean_data[3030, "numleaves"] <- NA
clean_data[3030, "LLL"] <- NA
clean_data[3031, "numleavesprev"] <- NA
clean_data[3031, "LLLprev"] <- NA
clean_data[3143, "numleaves"] <- NA
clean_data[3143, "LLL"] <- NA

#Case where LLL is too high given number of leaves
clean_data[3606, "numleaves"] <- NA
clean_data[3606, "LLL"] <- NA
clean_data[3607, "numleavesprev"] <- NA
clean_data[3607, "LLLPrev"] <- NA



### Combining LLL and numleaves----

#Using the sqrt of numleaves makes relationship more linear
#cor(sqrt(numleaves), LLL) = 0.7867325
#z-scaling both sqrt(numleaves) and LLL then adding together
#Should remove some of the noise and make a better size metric

clean_data <- clean_data %>% 
  mutate(scaled_sqrt_num_leaves = c(scale(sqrt(numleaves)))) %>% 
  mutate(scaled_LLL = c(scale(LLL))) %>% 
  mutate(scaled_sqrt_numleavesprev = c(scale(sqrt(numleavesprev)))) %>% 
  mutate(scaled_LLLprev = c(scale(LLLprev))) %>% 
  mutate(comp_size = scaled_LLL + scaled_sqrt_num_leaves) %>% 
  mutate(comp_size_prev = scaled_LLLprev + scaled_sqrt_numleavesprev)




saveRDS(clean_data, file = "Tg_data.rds")


### Reproductive dataset----

#Reproductive dataset
rep_data <- clean_data %>% 
  filter(NumFlowers>0)

saveRDS(rep_data, file = "Tg_rep.rds")



### Temperature values to climate data----
temp_data <- read_csv("phensyn_weather_2022.csv")
Climate_data_temp <- Climate_Data %>% 
  left_join(temp_data, by = "year")


### Combining demographic and climatic data----

#Fixing case of year and removing unnecessary columns
Climate_Data_fix <- Climate_data_temp %>% 
  rename("Year" = "year") %>% 
  select(-c(snowmelt.doy, snow.year.x, snow.year.y,
            fall.precip.mm.cb, fall.weather.n,
            summer.precip.mm.cb, summer.weather.n, 
            spring.weather.n,
            winterspring.precip, winter.weather.n))

full_data <- clean_data %>% 
  left_join(Climate_Data_fix, by = "Year")

saveRDS(full_data, file = "full_data.rds")

#Making climatic reproductive dataset
full_data_rep <- full_data %>% 
  filter(NumFlowers>0)

saveRDS(full_data_rep, file = "full_data_rep.rds")



add_climate_lag <- function(data, climate, time) {
  # Rename climate columns to indicate lag
  climate_lagged <- climate %>%
    rename_with(~ paste0(.x, "_lag", time), -Year)
  
  # Shift: subtract time from Year so year t in data joins to year t-time climate
  temp <- data %>%
    mutate(Year_shifted = Year - time) %>%
    left_join(climate_lagged, by = c("Year_shifted" = "Year")) %>%
    select(-Year_shifted)
  
  return(temp)
}

# Build up all lags in one dataset
Tg_climate_lag <- clean_data
for(i in 0:4){
  Tg_climate_lag <- add_climate_lag(Tg_climate_lag, Climate_Data_fix, i)
}

# Adding z-scores for snowpack
Tg_climate_lag <- Tg_climate_lag %>% 
  mutate(snowpack0_z = c(scale(snowpack_lag0))) %>% 
  mutate(snowpack1_z = c(scale(snowpack_lag1))) %>% 
  mutate(snowpack2_z = c(scale(snowpack_lag2))) %>% 
  mutate(snowpack3_z = c(scale(snowpack_lag3))) %>% 
  mutate(snowpack4_z = c(scale(snowpack_lag4))) %>%
  mutate(snowmelt0_z = c(scale(snowmelt_lag0))) %>% 
  mutate(snowmelt1_z = c(scale(snowmelt_lag1))) %>% 
  mutate(snowmelt2_z = c(scale(snowmelt_lag2))) %>% 
  mutate(snowmelt3_z = c(scale(snowmelt_lag3))) %>% 
  mutate(snowmelt4_z = c(scale(snowmelt_lag4)))

Tg_climate_clean <- Tg_climate_lag %>% 
  select(-c(Comments, logLLL, scaled_sqrt_num_leaves, scaled_sqrt_numleavesprev,
            scaled_LLL, scaled_LLLprev))


saveRDS(Tg_climate_clean, file = "Tg_climate_rds.rds")

#Reproductive Tg_climate
Tg_climate_rep <- Tg_climate %>% 
  filter(NumFlowers>0)

saveRDS(Tg_climate_rep, file = "Tg_climate_rep_rds.rds")



### Lagged climate years----
add_climate_lag_yr <- function(data, climate, time) {
  # Rename climate columns to indicate lag (excluding the join key 'Year')
  climate_lagged <- climate %>%
    rename_with(~ paste0(.x, "_lag", time), -Year)
  
  # Shift: subtract time from Year so year t in data joins to year t-time climate
  temp <- data %>%
    mutate(Year_shifted = Year - time) %>%
    left_join(climate_lagged, by = c("Year_shifted" = "Year")) %>%
    select(-Year_shifted)
  
  return(temp)
}

# Start with a copy of your base data
climate_data_with_lag <- Climate_Data

# Loop through lags 0 to 4 and update climate_data_with_lag sequentially
for (i in 0:4) {
  climate_data_with_lag <- add_climate_lag(
    data = climate_data_with_lag, 
    climate = Climate_Data, 
    time = i
  )
}
#Exclude Na
climate_data_with_lag <- climate_data_with_lag %>% 
  filter(Year>=1979 & Year < 2026)

saveRDS(climate_data_with_lag, file = "climate_data_with_lag.rds")
