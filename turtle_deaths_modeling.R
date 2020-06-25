rm(list=ls())

library(tidyverse)
#library(ggplot2)
library(lubridate)
library(readr)
library(reshape2)
library(jagsUI)
#library(bayesplot)

col.def <- cols(Date = col_date(format = "%m/%d/%Y"),
                Heat = col_integer(),
                Cliff = col_integer(),
                Total = col_integer())

# Begin and end dates for each season
Day.begin <- "10-15"
Day.end <- "04-30"
# define nesting seasons - not the same as calendar years
# also define day of season, starts on 10-15

data.1 <- read_csv(file = "data/Turtle_deaths_with_zeros.csv",
                   col_types = col.def) %>%
  mutate(Year = year(Date),
         Season = ifelse(Date >= paste0(Year, "-", Day.end), 
                         Year, Year-1),)

# find the number of days between Day1 and Day2
n.days <- as.numeric(as.Date(paste0("2018-", Day.end)) - as.Date(paste0("2017-", Day.begin)))

seasons <- unique(data.1$Season)

data.0 <- data.frame(Season = rep(seasons, each = n.days),
                     DOS = rep(seq(from = 1, to = n.days, by = 1), 
                               times = length(seasons))) %>% 
  mutate(Date = as.Date(paste0(Season, "-", Day.begin)) + days(DOS))

data.0 %>% left_join(data.1, by = "Date") -> data.0.1

obs.y <- array(NA, c(length(seasons), 2, n.days))

for (i in 1:length(seasons)){
  obs.y[i, ,] <- t(filter(data.0.1, Season.x == seasons[i]) %>% 
                     select(c("Heat", "Cliff")))
  
}

jags.data <- list(y = log(obs.y + 1),
                  mu.X0 = c(-2, -2), #array(data = -2, dim = c(dim(obs.y)[1], dim(obs.y)[2])),
                  day = seq(from = 1, to = n.days, by = 1),
                  n.days = length(seq(from = 1, to = n.days, by = 1)),
                  n.years = dim(obs.y)[1],
                  pi = pi)

jags.params <- c("beta.cos", "beta.sin",
                 #"sigma.X", 
                 "sigma.heat.y", "sigma.cliff.y", "rho",
                 "y", "X", "deviance")

MCMC.params <- list(n.chains = 4,
                    n.samples = 50000,
                    n.burnin = 10000,
                    n.thin = 5)


if (!file.exists(paste0("RData/DFS_X0_imputation_mvnorm_", Day.begin, "_", Day.end, ".rds"))){
  jm <- jags(jags.data,
             inits = NULL,
             parameters.to.save= jags.params,
             model.file = 'models/model_DFS_imputation_X0_mvnorm.txt',
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.thin = MCMC.params$n.thin,
             n.iter = MCMC.params$n.samples,
             DIC = T, parallel=T)
  
  saveRDS(jm, paste0("RData/DFS_X0_imputation_mvnorm_", Day.begin, "_", Day.end, ".rds"))
  
} else {
  jm <- readRDS(paste0("RData/DFS_X0_imputation_mvnorm_", Day.begin, "_", Day.end, ".rds"))
}
