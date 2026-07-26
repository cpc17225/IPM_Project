### Packages and data preparing----

library(tidyverse)
library(mclust)
library(readxl)
library(glmmTMB)


full_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")
Climate_Data <- read_csv("C:/Users/Owner/Downloads/barr_snowmelt_date_2022.csv")
Climate_Data <- Climate_Data %>% 
  rename("Year" = year)


#Find LLL and numleaves for first year of identification
first_appearance <- full_data%>%
  group_by(Tag) %>% 
  mutate(Tag = Tag) %>% 
  summarize(first_year = min(Year),
            first_size_LLL = LLL[Year == min(Year)],
            first_size_num = numleaves[Year == min(Year)])

# Recruits are plants first appearing after the earliest census year
# 1980 stills seems to have a lot of large plants discovered
recruits <- first_appearance %>%
  filter(first_year > 1980)

### Visualization----


# Recruit size distribution - fit a Gaussian or gamma
hist(recruits$first_size_LLL, breaks = 30)
hist(log(recruits$first_size_LLL), breaks = 50)

# Annual recruit counts----
# Filtered past 1980 since 1980 still identified many new plants that did not
#seem to be seedlings/juveniles
recruit_counts <- recruits %>%
  group_by(first_year) %>%
  filter(first_year>1980) %>% 
  summarize(n_recruits = n())

#Do the number of recruits vary by year?
ggplot(recruit_counts, aes(x = first_year, y = n_recruits)) +
  geom_point()

#Take away NA for ease of coding
recruits1 = na.omit(recruits)

### Distribution of recruit sizes----

#Density function instead of histogram
plot(density(log(recruits1$first_size_LLL)), 
     main = "Smooth Density of Log Recruit Size",
     col = "darkblue", lwd = 2)


#Use a Gaussian mixing model to see if the left-side peak is real
mix_model <- Mclust(log(recruits1$first_size_LLL))
summary(mix_model)
plot(mix_model, what = "BIC")

#Extract the statistics about the distributions
weights <- mix_model$parameters$pro
means <- mix_model$parameters$mean
sds <- sqrt(mix_model$parameters$variance$sigmasq)
#Into data frame
data.frame(Weight = weights, Mean_LogLLL = means, SD_LogLLL = sds)


#Distinguish pathways
plant_assignments <- mix_model$classification
length(plant_assignments)
recruits1$assign <- plant_assignments - 1
comp1_data <- log(recruits1$first_size_LLL)[plant_assignments == 1]
comp2_data <- log(recruits1$first_size_LLL)[plant_assignments == 2]

hist(comp1_data, breaks = 20)
hist(comp2_data, breaks = 20)

# Pathway 1 QQ plot
qqnorm(comp1_data, main = "Q-Q Plot: Component 1 (Year 2 Juveniles)")
qqline(comp1_data, col = "red", lwd = 2)

# Pathway 2 QQ Plot
qqnorm(comp2_data, main = "Q-Q Plot: Component 2 (Year 3 Juveniles (or later individuals))")
qqline(comp2_data, col = "red", lwd = 2)

head(mix_model$z)

tapply(recruits1$first_size_LLL,
       mix_model$classification,
       summary)

R1_individuals <- recruits1 %>% 
  filter(plant_assignments == 1)

model <- glmmTMB(assign ~ first_year,
                 family = binomial,
                 data = recruits1)
summary(model)

ggplot(recruits1, aes(x = first_year, y = assign))+
  geom_point()+
  geom_smooth(method = "gam", method.args = list(family = "binomial"))

recruit_climate <- recruits1 %>% 
  rename("Year" = first_year) %>% 
  left_join(Climate_Data, by = "Year")

ggplot(recruit_climate, aes(x = snowpack, y = assign))+
  geom_smooth(method = "gam", method.args = list(family = "binomial"))
ggplot(recruit_climate, aes(x = snowmelt, y = assign))+
  geom_smooth(method = "gam", method.args = list(family = "binomial"))

recruit_climate$P_large <- mix_model$z[,2]
hist(recruit_climate$P_large)

model_large <- lm(log(first_size_LLL) ~ snowmelt,
                  data = recruit_climate)
summary(model_large)
plot(model_large)

model_s <- glmmTMB(log(first_size_LLL) ~ s(snowmelt),
                   family = gaussian,
                   data = recruit_climate)
summary(model_s)
plot(simulateResiduals(model_s))

mu_s <- predict(model_s)

res <- log(recruit_climate$first_size_LLL) - mu_s

hist(res)
qqnorm(res)
density(res)

ggplot(recruit_climate, aes(x = snowmelt, y = log(first_size_LLL)))+
  geom_smooth(method = "gam", method.args = list(family = "gaussian"), formula = y ~ s(x, k = 5))
ggplot(recruit_climate, aes(x = snowpack, y = log(first_size_LLL)))+
  geom_smooth(method = "lm", formula = y ~ poly(x,3))
ggplot(recruit_climate, aes(x = scale(snowmelt * snowpack), y = log(first_size_LLL))) +
  geom_smooth(method = "lm", formula = y ~ poly(x,3))

model_recruit_c <- glmmTMB(log(first_size_LLL) ~ scale(snowmelt),
                           family = gaussian(),
                           dispformula = ~ scale(snowmelt),
                           data = recruit_climate)
summary(model_recruit_c)
plot(simulateResiduals(model_recruit_c))
