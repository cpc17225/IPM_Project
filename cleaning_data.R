### Packages and loading data----

library(tidyverse)
library(readxl)

Grandiflora_Cleaned <- read_excel("~/Grandiflora_Cleaned.xlsx", na = "NA")
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



saveRDS(clean_data, file = "Tg_data.rds")


### Reproductive dataset----

#Reproductive dataset
rep_data <- clean_data %>% 
  filter(NumFlowers>0)

saveRDS(rep_data, file = "Tg_rep.rds")


### Combining demographic and climatic data----

Climate_Data <- Climate_Data %>% 
  rename("Year" = "year")

full_data <- clean_data %>% 
  left_join(Climate_Data, by = "Year")

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
Tg_climate <- clean_data
for(i in 0:4){
  Tg_climate <- add_climate_lag(Tg_climate, Climate_Data, i)
}

# Adding z-scores for snowpack
Tg_climate <- Tg_climate %>% 
  mutate(snowpack0_z = c(scale(snowpack_lag0))) %>% 
  mutate(snowpack1_z = c(scale(snowpack_lag1))) %>% 
  mutate(snowpack2_z = c(scale(snowpack_lag2))) %>% 
  mutate(snowpack3_z = c(scale(snowpack_lag3))) %>% 
  mutate(snowpack4_z = c(scale(snowpack_lag4)))

saveRDS(Tg_climate, file = "Tg_climate_rds.rds")

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

# 1. Start with a copy of your base data
climate_data_with_lag <- Climate_Data

# 2. Loop through lags 0 to 4 and update climate_data_with_lag sequentially
for (i in 0:4) {
  climate_data_with_lag <- add_climate_lag(
    data = climate_data_with_lag, 
    climate = Climate_Data, 
    time = i
  )
}
saveRDS(climate_data_with_lag, file = "climate_data_with_lag.rds")
