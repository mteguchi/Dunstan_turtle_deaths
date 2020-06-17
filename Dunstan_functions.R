# why is this script not showing up on Git?


Girondot_fcn <- function(d, S, K, P, min, max){
  K <- abs(K)
  S <- abs(S)
  S1 <- -S
  M1 <- (1 + (2 * exp(K) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K))
  M2 <- (1 + (2 * exp(K) - 1) * exp((1/S) * (P - d))) ^ (-1/exp(K))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}

Girondot_fcn_2 <- function(d, S, K1, K2, P, min, max){
  K1 <- abs(K1)
  S <- abs(S)
  K2 <- abs(K2)
  S1 <- -S
  M1 <- (1 + (2 * exp(K1) - 1) * exp((1/S1) * (P - d))) ^ (-1/exp(K1))
  M2 <- (1 + (2 * exp(K2) - 1) * exp((1/S) * (P - d))) ^ (-1/exp(K2))
  N <- min + (max - min) * (M1 * M2)
  return(N)
}

# Extracting posterior samples of deviance or any other variable from jags output:
extract.samples <- function(varname, zm){
  dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
  return(dev)
}

compute.summary.stats_St_K <- function(samples, seasons, max.days, P){
  n.seasons <- length(seasons)
  K <- extract.samples("K", samples)   # just one of these
  X.low <- X.high <- X.med <- matrix(data = NA, 
                                     nrow = n.seasons, 
                                     ncol = max(max.days))
  i <- j <- 1
  for (i in 1:n.seasons){   # year
    max_i <- extract.samples(paste0("max[", i, "]"), samples)
    S <- extract.samples(paste0("S[", i, "]"), samples)
    for (j in 1:max.days[i]){   # days
      # state
      M1 <- (1 + (2 * exp(K) - 1) * exp((1/(-S)) * (P[i] - j))) ^ (-1/exp(K))
      M2 <- (1 + (2 * exp(K) - 1) * exp((1/S) * (P[i] - j))) ^ (-1/exp(K))
      X <-  max_i * (M1 * M2)
      
      X.low[i,j] <- quantile(X, 0.025)
      X.med[i,j] <- quantile(X, 0.50)      
      X.high[i,j] <- quantile(X, 0.975)
      
    }
  }  
  
  X.low.df <- data.frame(t(X.low)) %>% melt() %>% 
    transmute(Low = value)
  X.med.df <- data.frame(t(X.med)) %>% melt() %>%
    transmute(Med = value) 
  X.high.df <- data.frame(t(X.high)) %>% melt() %>%
    transmute(High = value) 
  
  X_St_K.df <- data.frame(DOS = rep(seq(1, ncol(X.low)), 
                                    length.out = n.seasons * ncol(X.low)),
                          Season = rep(seasons,
                                       each = ncol(X.low)),
                          Low = X.low.df$Low,
                          Med = X.med.df$Med,
                          High = X.high.df$High)  
  return(X_St_K.df)
}

compute.summary.stats_S_K <- function(samples, seasons, max.days, P){
  
  K <- extract.samples("K", samples)   # just one of these
  S <- extract.samples("S", samples)
  
  X.low <- X.high <- X.med <- matrix(data = NA, 
                                     nrow = n.seasons, 
                                     ncol = max(max.days))
  i <- j <- 1
  for (i in 1:n.seasons){   # year
    max_i <- extract.samples(paste0("max[", i, "]"), samples)
    for (j in 1:max.days[i]){   # days
      # state
      M1 <- (1 + (2 * exp(K) - 1) * exp((1/(-S)) * (P[i] - j))) ^ (-1/exp(K))
      M2 <- (1 + (2 * exp(K) - 1) * exp((1/S) * (P[i] - j))) ^ (-1/exp(K))
      X <-  max_i * (M1 * M2)
      
      X.low[i,j] <- quantile(X, 0.025)
      X.med[i,j] <- quantile(X, 0.50)      
      X.high[i,j] <- quantile(X, 0.975)
      
    }
  }  
  
  X.low.df <- data.frame(t(X.low)) %>% melt() %>% 
    transmute(Low = value)
  X.med.df <- data.frame(t(X.med)) %>% melt() %>%
    transmute(Med = value) 
  X.high.df <- data.frame(t(X.high)) %>% melt() %>%
    transmute(High = value) 
  
  X_S_K.df <- data.frame(DOS = rep(seq(1, ncol(X.low)), 
                                   length.out = n.seasons * ncol(X.low)),
                         Season = rep(seasons,
                                      each = ncol(X.low)),
                         Low = X.low.df$Low,
                         Med = X.med.df$Med,
                         High = X.high.df$High)
  
  return(X_S_K.df)
}

compute.summary.stats_St_K1_K2 <- function(samples, seasons, max.days, P){
  n.seasons <- length(seasons)
  K1 <- extract.samples("K1", samples)   # just one of these
  K2 <- extract.samples("K2", samples)   # just one of these
  X.low <- X.high <- X.med <- matrix(data = NA, 
                                     nrow = n.seasons, 
                                     ncol = max(max.days))
  i <- j <- 1
  for (i in 1:n.seasons){   # year
    max_i <- extract.samples(paste0("max[", i, "]"), samples)
    S <- extract.samples(paste0("S[", i, "]"), samples)
    for (j in 1:max.days[i]){   # days
      # state
      M1 <- (1 + (2 * exp(K1) - 1) * exp((1/(-S)) * (P[i] - j))) ^ (-1/exp(K1))
      M2 <- (1 + (2 * exp(K2) - 1) * exp((1/S) * (P[i] - j))) ^ (-1/exp(K2))
      X <-  max_i * (M1 * M2)
      
      X.low[i,j] <- quantile(X, 0.025)
      X.med[i,j] <- quantile(X, 0.50)      
      X.high[i,j] <- quantile(X, 0.975)
      
    }
  }  
  
  X.low.df <- data.frame(t(X.low)) %>% melt() %>% 
    transmute(Low = value)
  X.med.df <- data.frame(t(X.med)) %>% melt() %>%
    transmute(Med = value) 
  X.high.df <- data.frame(t(X.high)) %>% melt() %>%
    transmute(High = value) 
  
  X_St_K1_K2.df <- data.frame(DOS = rep(seq(1, ncol(X.low)), 
                                    length.out = n.seasons * ncol(X.low)),
                          Season = rep(seasons,
                                       each = ncol(X.low)),
                          Low = X.low.df$Low,
                          Med = X.med.df$Med,
                          High = X.high.df$High)  
  return(X_St_K1_K2.df)
}

compute.summary.stats_S_K1_K2 <- function(samples, seasons, max.days, P){
  
  K1 <- extract.samples("K1", samples)   # just one of these
  K2 <- extract.samples("K2", samples)   # just one of these
  S <- extract.samples("S", samples)
  
  X.low <- X.high <- X.med <- matrix(data = NA, 
                                     nrow = n.seasons, 
                                     ncol = max(max.days))
  i <- j <- 1
  for (i in 1:n.seasons){   # year
    max_i <- extract.samples(paste0("max[", i, "]"), samples)
    for (j in 1:max.days[i]){   # days
      # state
      M1 <- (1 + (2 * exp(K1) - 1) * exp((1/(-S)) * (P[i] - j))) ^ (-1/exp(K1))
      M2 <- (1 + (2 * exp(K2) - 1) * exp((1/S) * (P[i] - j))) ^ (-1/exp(K2))
      X <-  max_i * (M1 * M2)
      
      X.low[i,j] <- quantile(X, 0.025)
      X.med[i,j] <- quantile(X, 0.50)      
      X.high[i,j] <- quantile(X, 0.975)
      
    }
  }  
  
  X.low.df <- data.frame(t(X.low)) %>% melt() %>% 
    transmute(Low = value)
  X.med.df <- data.frame(t(X.med)) %>% melt() %>%
    transmute(Med = value) 
  X.high.df <- data.frame(t(X.high)) %>% melt() %>%
    transmute(High = value) 
  
  X_S_K1_K2.df <- data.frame(DOS = rep(seq(1, ncol(X.low)), 
                                   length.out = n.seasons * ncol(X.low)),
                         Season = rep(seasons,
                                      each = ncol(X.low)),
                         Low = X.low.df$Low,
                         Med = X.med.df$Med,
                         High = X.high.df$High)
  
  return(X_S_K1_K2.df)
}

