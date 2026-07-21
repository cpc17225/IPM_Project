### Packages and data preparing----

library(tidyverse)
library(mclust)

full_data <- readRDS("C:/Users/Owner/OneDrive/Desktop/R/IPM_Project/Tg_data.rds")

#Find LLL and numleaves for first year of identification
first_appearance <- full_data%>%
  group_by(Tag) %>%
  summarize(first_year = min(Year),
            first_size_LLL = LLL[Year == min(Year)],
            first_size_num = numleaves[Year == min(Year)])

# Recruits are plants first appearing after the earliest census year
# 1980 stills seems to have a lot of large plants discovered
recruits <- first_appearance %>%
  filter(first_year > 1980)

### Visualization----


# Recruit size distribution - fit a Gaussian or gamma
hist(recruits$first_size_LLL)
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
mix_model <- Mclust(log(recruits$first_size_LLL))
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
comp1_data <- log(recruits$first_size_LLL)[plant_assignments == 1]
comp2_data <- log(recruits$first_size_LLL)[plant_assignments == 2]

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

