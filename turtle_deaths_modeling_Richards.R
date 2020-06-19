
rm(list=ls())

library(tidyverse)
library(ggplot2)
library(lubridate)
library(readr)
library(reshape2)
library(jagsUI)
library(bayesplot)
library(loo)

source("Dunstan_functions.R")

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
         DOY = yday(Date),
         Season = ifelse(Date >= paste0(Year, "-", Day.end), 
                         Year, Year-1),
         Season.f = as.factor(Season), 
         DOS = as.numeric(Date - as.Date(paste0(Season, "-", Day.begin)))) %>%  
  group_by(Season.f) %>%
  mutate(Cum_Heat = cumsum(Heat),
         Cum_Cliff = cumsum(Cliff))

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

n.seasons <- length(unique(data.1$Season))
max.days <- data.1 %>% 
  group_by(Season.f) %>% 
  summarise(total.days = max(DOS),
            n.days = n()) 

obs.y.1 <- array(data = NA, dim = c(n.seasons, max(max.days$n.days), 2))
obs.t <- matrix(NA, nrow = n.seasons, ncol = max(max.days$n.days))
n.vec <- P <- vector(mode = "numeric", length = n.seasons)
#season.sector.idx <- matrix(NA, nrow = length(Sectors) * length(seasons), ncol = 2)
#max.vec <- vector(mode = "numeric", length = length(seasons))
max.mat <- matrix(NA, nrow = length(seasons), ncol = 2)

season_idx <- cause_idx <- 1
for (season_idx in 1:n.seasons){
  data.1 %>% filter(Season == seasons[season_idx]) -> tmp
  n.vec[season_idx] <- nrow(tmp)
  P[season_idx] <- floor(max(tmp$DOS)/2)
  
  for (cause_idx in 1:2){
    
    if (nrow(tmp) > 0){
      obs.t[season_idx, 1:nrow(tmp)] <- tmp$DOS
      tmp2 <- ifelse(cause_idx == 1,
                     tmp[, "Heat"],
                     tmp[, "Cliff"])
      
      obs.y.1[season_idx, 1:nrow(tmp), cause_idx] <- unlist(tmp2)
      
      max.mat[season_idx, cause_idx] <- ifelse(cause_idx == 1,
                                               max(tmp$Heat, na.rm = T),
                                               max(tmp$Cliff, na.rm = T))
    }
    
  }
}


MCMC.params <- list(n.chains = 4,
                    n.samples = 50000,
                    n.burnin = 10000,
                    n.thin = 5)

MCMC.params.2 <- list(n.chains = 4,
                    n.samples = 100000,
                    n.burnin = 60000,
                    n.thin = 5)

n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin

model.names <- c("S_K", "St_K", "S_K1_K2",
                 "St_K1_K2", "S1_S2_K", "S1_S2_K1_K2",
                 "S_K_P", "St_K_P", "S_K1_K2_P",
                 "St_K1_K2_P", "S1_S2_K_P", "S1_S2_K1_K2_P")

jags.parameters <- list(c("sigma.y", "S", "K", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K1", "K2", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K1", "K2", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S1", "S2", "K", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S1", "S2", "K1", "K2", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K", "P", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K", "P","max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K1", "K2", "P", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S", "K1", "K2", "P", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S1", "S2", "K", "P", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"),
                        c("sigma.y", "S1", "S2", "K1", "K2", "P", "max",
                          "y", "X", "sigma.max", "deviance", "loglik"))

model.setting.list <- list(jags.parameters = jags.parameters,
                          MCMC.params = MCMC.params,
                          model.names = model.names)

saveRDS(model.setting.list, file = "RData/model_setting_2020-06-17.rds")




# Heat related deaths


jags.data.Heat <- list(y = (obs.y.1[,,1]),
                       t = obs.t,
                       P = P,
                       n.vec = n.vec,
                       n.seasons = n.seasons,
                       max.vec = max.mat[,1])

jags.data.Heat_P <- list(y = (obs.y.1[,,1]),
                         t = obs.t,
                         n.vec = n.vec,
                         n.seasons = n.seasons,
                         max.vec = max.mat[,1])

# Convert the data into a vector also and match with the log-likelihood vector
# for computing looic.
data.heat.vector <- as.vector(jags.data.Heat$y) %>% 
  rep(each = MCMC.params$n.chains * n.per.chain)

loo.out.heat <- vector(mode = "list", length = length(model.names))
DIC.heat <- vector(mode = "numeric", length = length(model.names))
max.Rhat.heat <- vector(mode = "numeric", length = length(model.names))

for (k in 1:length(model.names)){
  if (!file.exists(paste0("RData/Richards_Heat_", model.names[[k]], "_norm.rds"))){
    
    if (grep("_P", model.names[[k]]) == 1){
      jags.data <- jags.data.Heat_P
    } else {
      jags.data <- jags.data.Heat
    }
    
    M1.Heat <- jags(jags.data,
                    inits = NULL,
                    parameters.to.save= jags.parameters[[k]],
                    model.file = paste0("models/model_Richards_", model.names[[k]], "_norm.txt"),
                    n.chains = MCMC.params$n.chains,
                    n.burnin = MCMC.params$n.burnin,
                    n.thin = MCMC.params$n.thin,
                    n.iter = MCMC.params$n.samples,
                    DIC = T, parallel=T)
    
    saveRDS(M1.Heat, paste0("RData/Richards_Heat_", model.names[[k]], "_norm.rds"))
    
  } else {
    M1.Heat <- readRDS(paste0("RData/Richards_Heat_", model.names[[k]], "_norm.rds"))
  }

  DIC.heat[k] <- M1.Heat$DIC
  max.Rhat.heat[k] <- data.frame(M1.Heat$summary) %>% select(Rhat) %>% max(na.rm = T)
  
  if (!file.exists(paste0("RData/Richards_Heat_", model.names[[k]], "_norm_LOOIC.rds"))){
    loo.out.heat[[k]] <- compute.looic(loglik = M1.Heat$sims.list$loglik,
                                       data.vector = data.heat.vector,
                                       MCMC.params = MCMC.params)
    
    saveRDS(loo.out.heat[[k]], paste0("RData/Richards_Heat_", model.names[[k]], 
                                      "_norm_LOOIC.rds"))
  } else {
    loo.out.heat[[k]] <- readRDS(paste0("RData/Richards_Heat_", 
                                        model.names[[k]], "_norm_LOOIC.rds"))
  }

}

LOOIC.Heat <- data.frame(model = c("M1", "M2", "M3", "M4", "M5", "M6",
                                   "M7", "M8", "M9", "M10", "M11", "M12"),
                         LOOIC = unlist(lapply(loo.out.heat, 
                                               FUN = function(x) x$estimates["looic",
                                                                             "Estimate"])),
                         SE = unlist(lapply(loo.out.heat, 
                                            FUN = function(x) x$estimates["looic", "SE"])),
                         DIC = DIC.heat,
                         max.Rhat = max.Rhat.heat) %>%
  arrange(by = LOOIC)

# A quick look at DICs
#DIC.Heat.df <- data.frame(Model = c("M1", "M2", "M3", "M4", "M5", "M6",
#                                    "M7", "M8", "M9", "M10", "M11", "M12"),
#                          DIC = DIC.heat,
#                          max.Rhat = max.Rhat.heat) %>%
#  arrange(by = DIC)
best.model.heat.ID <- 10


best.model.heat.looic <- readRDS(paste0("RData/Richards_Heat_", 
                                        model.names[[best.model.heat.ID]], 
                                        "_norm_LOOIC.rds"))

pareto_k_table(best.model.heat.looic)
    

best.model.heat <- readRDS(paste0("RData/Richards_Heat_", 
	model.names[[best.model.heat.ID]], "_norm.rds"))

# According to Rhat statistics, convergence reached for all. 

#mcmc_trace(best.model.heat$samples, "P")

#mcmc_dens(best.model.heat$samples, "P")



X_Heat_St_K1_K2_P.df <- compute.summary.stats_St_K1_K2_P(best.model.heat$samples, 
                                                         seasons, 
                                                         max.days$total.days)

ggplot() + 
  geom_path(data = X_Heat_St_K1_K2_P.df,
            aes(x = DOS, y = Med)) +
  geom_ribbon(data = X_Heat_St_K1_K2_P.df,
              aes(x = DOS, ymin = Low, ymax = High),
              alpha = 0.4) + 
  geom_point(data = data.1,
             aes(x = DOS, y = Heat)) + 
  facet_wrap(.~ as.factor(Season))


X_Heat_St_K1_K2_P.df %>% 
  group_by(Season) %>%
  summarise(total.low = sum(ceiling(Low), na.rm = T),
            total.med = sum(ceiling(Med), na.rm = T),
            total.high = sum(ceiling(High), na.rm = T)) -> total.deaths.heat

# first just run one at a time
jags.data.Cliff <- list(y = (obs.y.1[,,2]),
                       t = obs.t,
                       P = P,
                       n.vec = n.vec,
                       n.seasons = n.seasons,
                       max.vec = max.mat[,1])

jags.data.Cliff_P <- list(y = (obs.y.1[,,2]),
                          t = obs.t,
                          n.vec = n.vec,
                          n.seasons = n.seasons,
                          max.vec = max.mat[,1])

# Convert the data into a vector also and match with the log-likelihood vector
# for computing looic.
data.cliff.vector <- as.vector(jags.data.Cliff$y) %>% 
  rep(each = MCMC.params$n.chains * n.per.chain)

loo.out.cliff <- vector(mode = "list", length = length(model.names))
DIC.cliff <- vector(mode = "numeric", length = length(model.names))
max.Rhat.cliff <- vector(mode = "numeric", length = length(model.names))

for (k in 1:length(model.names)){
  
  if (!file.exists(paste0("RData/Richards_Cliff_", model.names[[k]], "_norm.rds"))){
    if (grep("_P", model.names[[k]]) == 1){
      jags.data <- jags.data.Cliff_P
    } else {
      jags.data <- jags.data.Cliff
    }
    
    M1.Cliff <- jags(jags.data,
                     inits = NULL,
                     parameters.to.save= jags.parameters[[k]],
                     model.file = paste0("models/model_Richards_", 
                                         model.names[[k]], "_norm.txt"),
                     n.chains = MCMC.params.2$n.chains,
                     n.burnin = MCMC.params.2$n.burnin,
                     n.thin = MCMC.params.2$n.thin,
                     n.iter = MCMC.params.2$n.samples,
                     DIC = T, parallel=T)
    
    saveRDS(M1.Cliff, paste0("RData/Richards_Cliff_", model.names[[k]], "_norm.rds"))
    
  } else {
    M1.Cliff <- readRDS(paste0("RData/Richards_Cliff_", model.names[[k]], "_norm.rds"))
  }
  
  DIC.cliff[k] <- M1.Cliff$DIC
  max.Rhat.cliff[k] <- data.frame(M1.Cliff$summary) %>% select(Rhat) %>% max(na.rm = TRUE)
  
  if (!file.exists(paste0("RData/Richards_Cliff_", model.names[[k]], "_norm_LOOIC.rds"))){
    #data.cliff.
    
    loo.out.cliff[[k]] <- compute.looic(loglik = M1.Cliff$sims.list$loglik,
                                        data.vector = data.cliff.vector,
                                        MCMC.params = MCMC.params.2)
    
    saveRDS(loo.out.cliff[[k]], paste0("RData/Richards_Cliff_", model.names[[k]], 
                                       "_norm_LOOIC.rds"))
  } else {
    loo.out.cliff[[k]] <- readRDS(paste0("RData/Richards_Cliff_", 
                                         model.names[[k]], "_norm_LOOIC.rds"))
  }

}

LOOIC.Cliff <- data.frame(model = c("M1", "M2", "M3", "M4", "M5", "M6",
                                    "M7", "M8", "M9", "M10", "M11", "M12"),
                          LOOIC = unlist(lapply(loo.out.cliff, 
                                                FUN = function(x) x$estimates["looic",
                                                                              "Estimate"])),

                          SE = unlist(lapply(loo.out.cliff, 
                                             FUN = function(x) x$estimates["looic", "SE"])),
                          DIC = DIC.cliff,
                          Max.Rhat = max.Rhat.cliff) %>%
  arrange(by = LOOIC)

# A quick look at DICs
#DIC.Cliff.df <- data.frame(Model = c("M1", "M2", "M3", "M4", "M5", "M6",
#                                     "M7", "M8", "M9", "M10", "M11", "M12"),
#                           DIC = DIC.cliff,
#                           Max.Rhat = max.Rhat.cliff) %>%
#  arrange(by = DIC)

best.model.cliff.ID <- 8

best.model.cliff.looic <- readRDS(paste0("RData/Richards_Cliff_", 
                                        model.names[[best.model.cliff.ID]], "_norm_LOOIC.rds"))

pareto_k_table(best.model.cliff.looic)
    

best.model.cliff <- readRDS(paste0("RData/Richards_Cliff_", 
	model.names[[best.model.cliff.ID]], "_norm.rds"))

X_Cliff_St_K1_K2_P.df <- compute.summary.stats_St_K1_K2_P(best.model.cliff$samples, 
                                                          seasons, 
                                                          max.days$total.days)

ggplot() + 
  geom_path(data = X_Cliff_St_K1_K2_P.df,
            aes(x = DOS, y = Med)) +
  geom_ribbon(data = X_Cliff_St_K1_K2_P.df,
              aes(x = DOS, ymin = Low, ymax = High),
              alpha = 0.4) + 
  geom_point(data = data.1,
             aes(x = DOS, y = Cliff)) + 
  facet_wrap(.~ as.factor(Season))


X_Cliff_St_K1_K2_P.df %>% 
  group_by(Season) %>%
  summarise(total.low = sum(ceiling(Low), na.rm = T),
            total.med = sum(ceiling(Med), na.rm = T),
            total.high = sum(ceiling(High), na.rm = T)) -> total.deaths.cliff


mcmc_trace(best.model.cliff$samples, "P")
