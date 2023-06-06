####################################################################################################################
####################################################################################################################
#------------------------------------- Packages --------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library(lubridate) ### Time-variables
library(expm) ### Matrix powers
library(MASS)
library(mixR) ### Finite mixture models
library(Rcpp) ### Finite mixture models
library(hmmTMB) ### Fitting
library(splines2) ### Periodic spline
library(qdapRegex) ### Periodic spline
library(loo) ### WAIC for Bayesian models
library(rstan) ### Traceplots


####################################################################################################################
####################################################################################################################
####################################################################################################################

### Access to fitted models and posterior samples saved under names specified as above
load("/Users/carolinelestrom/Documents/KU/Bachelor/Data/Modelfit/sims_02_06.RData")

####################################################################################################################
####################################################################################################################
####################################################################################################################










####################################################################################################################
#------------------------------------- State-dependent densities ---------------------------------------------------
####################################################################################################################



#------------------------------------- Per Frequentist - TOY --------------------------------------------------------

perf_df <- data.frame(Y = 0.5, X = 450)

pfdens_TOY <- ggplot(perf_df, aes(X, Y)) +
  geom_point(size = 1, alpha = 1, color = "#901a1E") +
  geom_abline(aes(intercept = 0, slope = 1/900), linewidth = 2, color = "#901a1E") +
  geom_abline(aes(intercept = 1, slope = -1/900), linewidth = 2, color = "#901a1E") +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Per - TOY (F)") +
  theme(legend.text = element_text(size = 9)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL, limits = c(0,1))

#------------------------------------- Per Bayesian - TOY --------------------------------------------------------------

n_iter <- 200
### Weights for state-specific density functions
w <- table(bayes_30_TOY$viterbi())/nrow(d30)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_30_TOY$iters()[i,1:8]
  dens1 <- obs_30_TOY$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                          sd = par[3])
  dens2 <- obs_30_TOY$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                          sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines
pbdens_TOY <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d30,
                 col = "#901a1E",fill = "white", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Per - TOY (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Gorm Frequentist - TOY -------------------------------------------------------------

### Get state-dependent parameters for a few values of TOY
par <- freq_17_TOY$predict(what = "obspar")

### Weights for state-dependent distributions
w <- table(freq_17_TOY$viterbi())/nrow(d17)
### Grid of depths (x axis)
grid <- seq(0, 800, length = 100)

### For each value of TOY, compute state-dependent distributions
pdf_ls <- lapply(1:dim(par)[3], function(i) {
  p <- par[,,i]
  pdf1 <- w[1] * dgamma(grid, shape = p["Depth.mean", "state 1"]^2/p["Depth.sd", "state 1"]^2,
                        rate = p["Depth.mean", "state 1"]/p["Depth.sd", "state 1"]^2)
  pdf2 <- w[2] * dgamma(grid, shape = p["Depth.mean", "state 2"]^2/p["Depth.sd", "state 2"]^2,
                        rate = p["Depth.mean", "state 2"]/p["Depth.sd", "state 2"]^2)
  res <- data.frame(Depth = grid, pdf = c(pdf1, pdf2, pdf1 + pdf2),
                    State = factor(c(paste0("State ", rep(1:2, each = length(grid))), paste0(rep("Total", each = length(grid))))))
  return(res) })
### Turn into dataframe for ggplot
pdf_df <- do.call(rbind, pdf_ls)

gfdens_TOY <- ggplot(pdf_df, aes(Depth, pdf)) +
  geom_histogram(aes(x = Depth, y=..density..), bins = 30,
                 fill = "white", col = "#901a1E", data = d17) +
  geom_line(aes(col = State), linewidth = 1, alpha = 0.8) +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Gorm - TOY (F)") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.73)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Gorm Bayesian - TOY --------------------------------------------------------------

n_iter <- 2000
### Weights for state-specific density functions
w <- table(bayes_17_TOY$viterbi())/nrow(d17)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_17_TOY$iters()[i,1:8]
  dens1 <- obs_17_TOY$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                          sd = par[3])
  dens2 <- obs_17_TOY$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                          sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines
gbdens_TOY <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d17,
                 col = "#901a1E",fill = "white", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Gorm - TOY (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)


#------------------------------------- Ragnar Frequentist - TOY -------------------------------------------------------------

### Get state-dependent parameters for a few values of TOY
par <- freq_6_TOY$predict(what = "obspar")

### Weights for state-dependent distributions
w <- table(freq_6_TOY$viterbi())/nrow(d63)
### Grid of depths (x axis)
grid <- seq(0, 1000, length = 100)

### For each value of TOY, compute state-dependent distributions
pdf_ls <- lapply(1:dim(par)[3], function(i) {
  p <- par[,,i]
  pdf1 <- w[1] * dgamma(grid, shape = p["Depth.mean", "state 1"]^2/p["Depth.sd", "state 1"]^2,
                        rate = p["Depth.mean", "state 1"]/p["Depth.sd", "state 1"]^2)
  pdf2 <- w[2] * dgamma(grid, shape = p["Depth.mean", "state 2"]^2/p["Depth.sd", "state 2"]^2,
                        rate = p["Depth.mean", "state 2"]/p["Depth.sd", "state 2"]^2)
  res <- data.frame(Depth = grid, pdf = c(pdf1, pdf2, pdf1 + pdf2),
                    State = factor(c(paste0("State ", rep(1:2, each = length(grid))), paste0(rep("Total", each = length(grid))))))
  return(res) })
### Turn into dataframe for ggplot
pdf_df <- do.call(rbind, pdf_ls)

rfdens_TOY <- ggplot(pdf_df, aes(Depth, pdf)) +
  geom_histogram(aes(x = Depth, y=..density..), bins = 30,
                 fill = "white", col = "#901a1E", data = d63) +
  geom_line(aes(col = State), linewidth = 1, alpha = 0.9) +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Ragnar - TOY (F)") + xlim(c(0,800)) +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.73)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Ragnar Bayesian - TOY ------------------------------------------------------------

n_iter <- 2000
### Weights for state-specific density functions
w <- table(bayes_6_TOY$viterbi())/nrow(d63)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_6_TOY$iters()[i,1:8]
  dens1 <- obs_6_TOY$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                         sd = par[3])
  dens2 <- obs_6_TOY$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                         sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines
rbdens_TOY <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d63,
                 col = "#901a1E",fill = "white", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Ragnar - TOY (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)





#------------------------------------- Per Frequentist - TOD -------------------------------------------------------------


### Get state-dependent parameters for a few values of TOD
par <- freq_30_TOD$predict(what = "obspar")

### Weights for state-dependent distributions
w <- table(freq_30_TOD$viterbi())/nrow(d30)
### Grid of depths (x axis)
grid <- seq(0, 800, length = 100)

### For each value of TOD, compute state-dependent distributions
pdf_ls <- lapply(1:dim(par)[3], function(i) {
  p <- par[,,i]
  pdf1 <- w[1] * dgamma(grid, shape = p["Depth.mean", "state 1"]^2/p["Depth.sd", "state 1"]^2,
                        rate = p["Depth.mean", "state 1"]/p["Depth.sd", "state 1"]^2)
  pdf2 <- w[2] * dgamma(grid, shape = p["Depth.mean", "state 2"]^2/p["Depth.sd", "state 2"]^2,
                        rate = p["Depth.mean", "state 2"]/p["Depth.sd", "state 2"]^2)
  res <- data.frame(Depth = grid, pdf = c(pdf1, pdf2, pdf1 + pdf2),
                    State = factor(c(paste0("State ", rep(1:2, each = length(grid))), paste0(rep("Total", each = length(grid))))))
  return(res) })
### Turn into dataframe for ggplot
pdf_df <- do.call(rbind, pdf_ls)

pf_TOD <- ggplot(pdf_df, aes(Depth, pdf)) +
  geom_histogram(aes(x = Depth, y=..density..), bins = 30,
                 fill = "white", col = "#901a1E", data = d30) +
  geom_line(aes(col = State), linewidth = 1, alpha = 0.9) +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Per - TOD (F)") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.73)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Per Bayesian - TOD ------------------------------------------------------------

n_iter <- 2000
### Weights for state-specific density functions
w <- table(bayes_30_TOD$viterbi())/nrow(d30)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_30_TOD$iters()[i,1:8]
  dens1 <- obs_30_TOD$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                          sd = par[3])
  dens2 <- obs_30_TOD$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                          sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines

pbdens_TOD <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d30,
                 fill = "white",col = "#901a1E", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Per - TOD (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Gorm Frequentist - TOD ------------------------------------------------------------

### Get state-dependent parameters for a few values of TOD
par <- freq_17_TOD$predict(what = "obspar")

### Weights for state-dependent distributions
w <- table(freq_17_TOD$viterbi())/nrow(d17)
### Grid of depths (x axis)
grid <- seq(0, 800, length = 100)

### For each value of TOD, compute state-dependent distributions
pdf_ls <- lapply(1:dim(par)[3], function(i) {
  p <- par[,,i]
  pdf1 <- w[1] * dgamma(grid, shape = p["Depth.mean", "state 1"]^2/p["Depth.sd", "state 1"]^2,
                        rate = p["Depth.mean", "state 1"]/p["Depth.sd", "state 1"]^2)
  pdf2 <- w[2] * dgamma(grid, shape = p["Depth.mean", "state 2"]^2/p["Depth.sd", "state 2"]^2,
                        rate = p["Depth.mean", "state 2"]/p["Depth.sd", "state 2"]^2)
  res <- data.frame(Depth = grid, pdf = c(pdf1, pdf2, pdf1 + pdf2),
                    State = factor(c(paste0("State ", rep(1:2, each = length(grid))), paste0(rep("Total", each = length(grid))))))
  return(res) })
### Turn into dataframe for ggplot
pdf_df <- do.call(rbind, pdf_ls)

gf_TOD <- ggplot(pdf_df, aes(Depth, pdf)) +
  geom_histogram(aes(x = Depth, y=..density..), bins = 30,
                 fill = "white", col = "#901a1E", data = d17) +
  geom_line(aes(col = State), linewidth = 1, alpha = 0.9) +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Gorm - TOD (F)") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.73)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Gorm Bayesian - TOD -----------------------------------------------------------

n_iter <- 2000
### Weights for state-specific density functions
w <- table(bayes_17_TOD$viterbi())/nrow(d17)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_17_TOD$iters()[i,1:8]
  dens1 <- obs_17_TOD$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                          sd = par[3])
  dens2 <- obs_17_TOD$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                          sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines
gbdens_TOD <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d17,
                 color = "#901a1E",fill = "white", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Gorm - TOD (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Ragnar Frequentist - TOD -------------------------------------------------------------

### Get state-dependent parameters for a few values of TOD
par <- freq_6_TOD$predict(what = "obspar")

### Weights for state-dependent distributions
w <- table(freq_6_TOD$viterbi())/nrow(d63)
### Grid of depths (x axis)
grid <- seq(0, 1000, length = 100)

### For each value of TOD, compute state-dependent distributions
pdf_ls <- lapply(1:dim(par)[3], function(i) {
  p <- par[,,i]
  pdf1 <- w[1] * dgamma(grid, shape = p["Depth.mean", "state 1"]^2/p["Depth.sd", "state 1"]^2,
                        rate = p["Depth.mean", "state 1"]/p["Depth.sd", "state 1"]^2)
  pdf2 <- w[2] * dgamma(grid, shape = p["Depth.mean", "state 2"]^2/p["Depth.sd", "state 2"]^2,
                        rate = p["Depth.mean", "state 2"]/p["Depth.sd", "state 2"]^2)
  res <- data.frame(Depth = grid, pdf = c(pdf1, pdf2, pdf1 + pdf2),
                    State = factor(c(paste0("State ", rep(1:2, each = length(grid))), paste0(rep("Total", each = length(grid))))))
  return(res) })
### Turn into dataframe for ggplot
pdf_df <- do.call(rbind, pdf_ls)

rf_TOD <- ggplot(pdf_df, aes(Depth, pdf)) +
  geom_histogram(aes(x = Depth, y=..density..), bins = 30,
                 fill = "white", col = "#901a1E", data = d63) +
  geom_line(aes(col = State), linewidth = 1, alpha = 0.9) +
  theme_bw() +
  labs(x = "Depth", y = "Density") +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  ggtitle("Ragnar - TOD (F)") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.73)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)

#------------------------------------- Ragnar Bayesian - TOD -------------------------------------------------------------

n_iter <- 2000
### Weights for state-specific density functions
w <- table(bayes_6_TOD$viterbi())/nrow(d63)

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))

### For each posterior sample, compute gamma distribution
dens_Depth_df <- data.frame()
Depth <- seq(0.4, 900, length = 100)

for(i in ind_post) {
  par <- bayes_6_TOD$iters()[i,1:8]
  dens1 <- obs_6_TOD$dists()$Depth$pdf()(x = Depth, mean = par[1],
                                         sd = par[3])
  dens2 <- obs_6_TOD$dists()$Depth$pdf()(x = Depth, mean = par[2],
                                         sd = par[4])
  dens_Depth <- data.frame(state = c(paste0("State ", rep(1:2, each = 100)), paste0(rep("Total", each = 100))),
                           dens_Depth = c(w[1] * dens1, w[2] * dens2, w[1] * dens1 + w[2] * dens2))
  
  dens_Depth$group <- paste0("iter ", i, " - ", dens_Depth$state)
  dens_Depth_df <- rbind(dens_Depth_df, dens_Depth)
}
dens_Depth_df$Depth <- Depth

### Plot histogram of Depth and density lines
rbdens_TOD <- ggplot(dens_Depth_df, aes(Depth, dens_Depth)) +
  geom_histogram(aes(y = ..density..), data = d63,
                 col = "#901a1E", fill = "white", breaks = seq(0, 900, length = 30)) +
  geom_line(aes(col = state, group = group), linewidth = 1, alpha = 0.5) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() + ggtitle("Ragnar - TOD (B)") + ylab("Density") +
  theme(legend.text = element_text(size = 11)) +
  theme(legend.position = c(0.85, 0.75)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  theme(axis.title = element_text(size=17)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0,900)) +
  scale_y_continuous(breaks=NULL)




#------------------------------------- Collected plot ---------------------------------------------------------------



grid.arrange(pfdens_TOY, pbdens_TOY, pf_TOD, pbdens_TOD,
             gfdens_TOY, gbdens_TOY, gf_TOD, gbdens_TOD,
             rfdens_TOY, rbdens_TOY, rf_TOD, rbdens_TOD,
             ncol = 4)





























####################################################################################################################
####################################################################################################################
#------------------------------------- TOY model -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################


####################################################################################################################
#------------------------------------- Stationary state probability plot -------------------------------------------
####################################################################################################################





#------------------------------------- Per Frequentist -------------------------------------------------------------

perf_df <- data.frame(Y = 0.5, X = 7)
pf_state <- ggplot(perf_df, aes(X, Y)) +
  geom_point(size = 1, alpha = 1, color = "#901a1E") +
  geom_abline(aes(intercept = -1/13, slope = 1/12.1), linewidth = 2, color = "#901a1E") +
  geom_abline(aes(intercept = 14/13, slope = -1/12.1), linewidth = 2, color = "#901a1E") +
  labs(x = "TOY", y = "Stationary state probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw() +
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                     label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1)) +
  ggtitle("Per - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Per Bayesian -----------------------------------------------------------

n_iter <- 200
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of year
probs_df <- data.frame()
newdata <- data.frame(TOY = seq(1, 13, length = 100))
for(i in ind_post) {
  bayes_30_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_30_TOY$predict(what = "delta",
                                                            newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOY <- newdata$TOY

### Plot stationary state probs against time of year
pb_state <- ggplot(probs_df, aes(TOY, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  labs(x = "Time Of Year", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Per - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                     label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Gorm Frequentist -------------------------------------------------------------

gf_state <- freq_17_TOY$plot(what = "delta", var = "TOY") +
  theme(legend.position="top") +
  ggtitle("Gorm - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(7, 8, 9, 9.9), label = c("July", "Aug", "Sep", "Oct")) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Gorm Bayesian -------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of year
probs_df <- data.frame()
newdata <- data.frame(TOY = seq(6.8, 10, length = 100))
for(i in ind_post) {
  bayes_17_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_17_TOY$predict(what = "delta",
                                                            newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOY <- newdata$TOY

### Plot stationary state probs against time of year
gb_state <- ggplot(probs_df, aes(TOY, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  labs(x = "Time Of Year", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Gorm - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Ragnar Frequentist -----------------------------------------------------------

rf_state <- freq_6_TOY$plot(what = "delta", var = "TOY") +
  theme(legend.position="top") +
  ggtitle("Ragnar - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(5.3, 6, 7, 8, 9, 10, 10.9), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Ragnar Bayesian -------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of year
probs_df <- data.frame()
newdata <- data.frame(TOY = seq(5, 11, length = 100))
for(i in ind_post) {
  bayes_6_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_6_TOY$predict(what = "delta",
                                                           newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOY <- newdata$TOY

### Plot stationary state probs against time of year
rb_state <- ggplot(probs_df, aes(TOY, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  labs(x = "Time Of Year", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Ragnar - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(5.1, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
  theme(legend.text=element_text(size=17))






#------------------------------------- Collected plot --------------------------------------------------------------



grid.arrange(pf_state, gf_state, rf_state,
             pb_state, gb_state, rb_state,
             ncol = 3)





















####################################################################################################################
#------------------------------------- Transition probability matrix plot ------------------------------------------
####################################################################################################################




#------------------------------------- Per Frequentist -------------------------------------------------------------

perf_df <- data.frame(Y = 0.5, X = 7)

pf_toytpm <- grid.arrange(
  ggplot(perf_df, aes(X, Y)) +
    geom_point(size = 1, alpha = 1, color = "#901a1E") +
    geom_abline(aes(intercept = -1/13, slope = 1/12.1), linewidth = 1, color = "#901a1E") +
    geom_abline(aes(intercept = 14/13, slope = -1/12.1), linewidth = 1, color = "#901a1E") +
    labs(x = "Time of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                       label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.5, 1), limits = c(0,1)) +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)),
  ggplot(perf_df, aes(X, Y)) +
    geom_point(size = 1, alpha = 1, color = "#901a1E") +
    geom_abline(aes(intercept = -1/13, slope = 1/12.1), linewidth = 1, color = "#901a1E") +
    geom_abline(aes(intercept = 14/13, slope = -1/12.1), linewidth = 1, color = "#901a1E") +
    labs(x = "Time of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                       label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.5, 1), limits = c(0,1)) +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)),
  ggplot(perf_df, aes(X, Y)) +
    geom_point(size = 1, alpha = 1, color = "#901a1E") +
    geom_abline(aes(intercept = -1/13, slope = 1/12.1), linewidth = 1, color = "#901a1E") +
    geom_abline(aes(intercept = 14/13, slope = -1/12.1), linewidth = 1, color = "#901a1E") +
    labs(x = "Time of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                       label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.5, 1), limits = c(0,1)) +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)),
  ggplot(perf_df, aes(X, Y)) +
    geom_point(size = 1, alpha = 1, color = "#901a1E") +
    geom_abline(aes(intercept = -1/13, slope = 1/12.1), linewidth = 1, color = "#901a1E") +
    geom_abline(aes(intercept = 14/13, slope = -1/12.1), linewidth = 1, color = "#901a1E") +
    labs(x = "Time of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), limits = c(1,13),
                       label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.5, 1), limits = c(0,1)) +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)),
  ncol = 2, top = grid::textGrob("Per - Frequentist", x = 0.33, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Per Bayesian --------------------------------------------------------------

n_iter <- 200
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of year
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOY = seq(1, 13, length = 100))
for(i in ind_post) {
  bayes_30_TOY$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_30_TOY$predict(what = "tpm",
                                                               newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_30_TOY$predict(what = "tpm",
                                                                 newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_30_TOY$predict(what = "tpm",
                                                                 newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_30_TOY$predict(what = "tpm",
                                                               newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOY <- newdata$TOY
probs12_df$TOY <- newdata$TOY
probs21_df$TOY <- newdata$TOY
probs22_df$TOY <- newdata$TOY

### Plot stationary state probs against time of year
pb_toytpm <- grid.arrange(
  ggplot(probs11_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ggplot(probs12_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs21_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0, 0.1)),
  ggplot(probs22_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12), label = c("Feb", "Apr", "June", "Aug", "Oct", "Dec")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ncol = 2, top = grid::textGrob("Per - Bayesian", x = 0.3, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Gorm Frequentist -------------------------------------------------------------

probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
### Compute transition probabilities over
### grid of time of year
newdata <- data.frame(TOY = seq(6.8, 10, length = 100))

probs11 <- data.frame(prob = as.vector((freq_17_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector((freq_17_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector((freq_17_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs11_df <- rbind(probs11_df, probs11)
probs11_df$TOY <- newdata$TOY

probs12 <- data.frame(prob = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs12_df <- rbind(probs12_df, probs12)
probs12_df$TOY <- newdata$TOY

probs21 <- data.frame(prob = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector(1-(freq_17_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs21_df <- rbind(probs21_df, probs21)
probs21_df$TOY <- newdata$TOY

probs22 <- data.frame(prob = as.vector((freq_17_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector((freq_17_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector((freq_17_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs22_df <- rbind(probs22_df, probs22)
probs22_df$TOY <- newdata$TOY

### Plot trans probs against time of year
gf_toytpm <- grid.arrange(
  ggplot(probs11_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(6, 10, by = 4)) +
    labs(x = "Time Of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0.75, 0.875, 1), limits = c(0.75,1)),
  ggplot(probs12_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0, 0.125, 0.25), limits = c(0, 0.25)),
  ggplot(probs21_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0, 0.125, 0.25), limits = c(0,0.25)),
  ggplot(probs22_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0.75, 0.875, 1), limits = c(0.75,1)),
  ncol = 2, top = grid::textGrob("Gorm - Frequentist", x = 0.33, gp = grid::gpar(fontsize = 23)))


#------------------------------------- Gorm Bayesian --------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute transition probabilities over
### grid of time of year
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOY = seq(6.8, 10, length = 100))
for(i in ind_post) {
  bayes_17_TOY$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_17_TOY$predict(what = "tpm",
                                                               newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_17_TOY$predict(what = "tpm",
                                                                 newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_17_TOY$predict(what = "tpm",
                                                                 newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_17_TOY$predict(what = "tpm",
                                                               newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOY <- newdata$TOY
probs12_df$TOY <- newdata$TOY
probs21_df$TOY <- newdata$TOY
probs22_df$TOY <- newdata$TOY

### Plot trans probs against time of year
gb_toytpm <- grid.arrange(
  ggplot(probs11_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0.75, 0.875, 1), limits = c(0.75,1)),
  ggplot(probs12_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0, 0.125, 0.25), limits = c(0,0.25)),
  ggplot(probs21_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0, 0.125, 0.25), limits = c(0, 0.25)),
  ggplot(probs22_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(7, 8, 9, 10), label = c("July", "Aug", "Sep", "Oct")) +
    scale_y_continuous(breaks=c(0.75, 0.875, 1), limits = c(0.75,1)),
  ncol = 2, top = grid::textGrob("Gorm - Bayesian", x = 0.3, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Ragnar Frequentist -------------------------------------------------------------

probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
### Compute transition probabilities over
### grid of time of year
newdata <- data.frame(TOY = seq(5, 11, length = 100))

probs11 <- data.frame(prob = as.vector((freq_6_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector((freq_6_TOY$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector((freq_6_TOY$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs11_df <- rbind(probs11_df, probs11)
probs11_df$TOY <- newdata$TOY

probs12 <- data.frame(prob = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs12_df <- rbind(probs12_df, probs12)
probs12_df$TOY <- newdata$TOY

probs21 <- data.frame(prob = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector(1-(freq_6_TOY$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs21_df <- rbind(probs21_df, probs21)
probs21_df$TOY <- newdata$TOY

probs22 <- data.frame(prob = as.vector((freq_6_TOY$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector((freq_6_TOY$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector((freq_6_TOY$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs22_df <- rbind(probs22_df, probs22)
probs22_df$TOY <- newdata$TOY

### Plot trans probs against time of year
rf_toytpm <- grid.arrange(
  ggplot(probs11_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(5, 11, by = 4)) +
    labs(x = "Time Of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ggplot(probs12_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs21_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(5, 11, by = 4)) +
    labs(x = "Time Of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs22_df, aes(TOY, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(5, 11, by = 4)) +
    labs(x = "Time Of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ncol = 2, top = grid::textGrob("Ragnar - Frequentist", x = 0.33, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Ragnar Bayesian -------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute transition probabilities over
### grid of time of year
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOY = seq(5, 11, length = 100))
for(i in ind_post) {
  bayes_6_TOY$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_6_TOY$predict(what = "tpm",
                                                              newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_6_TOY$predict(what = "tpm",
                                                                newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_6_TOY$predict(what = "tpm",
                                                                newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_6_TOY$predict(what = "tpm",
                                                              newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOY <- newdata$TOY
probs12_df$TOY <- newdata$TOY
probs21_df$TOY <- newdata$TOY
probs22_df$TOY <- newdata$TOY

### Plot trans probs against time of year
rb_toytpm <- grid.arrange(
  ggplot(probs11_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ggplot(probs12_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs21_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs22_df, aes(TOY, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    labs(x = "Time Of Year", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=13, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(5, 6, 7, 8, 9, 10, 11), label = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov")) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ncol = 2, top = grid::textGrob("Ragnar - Bayesian", x = 0.33, gp = grid::gpar(fontsize = 23)))




#------------------------------------- Collected plot --------------------------------------------------------------





grid.arrange(pf_toytpm, pb_toytpm,
             gf_toytpm, gb_toytpm,
             rf_toytpm, rb_toytpm,
             ncol = 2)





























####################################################################################################################
####################################################################################################################
#------------------------------------- TOD model -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################



####################################################################################################################
#------------------------------------- Stationary state probability plot -------------------------------------------
####################################################################################################################


#------------------------------------- Per Frequentist -------------------------------------------------------------

perf <- freq_30_TOD$plot(what = "delta", var = "TOD") +
  theme(legend.position="top") +
  ggtitle("Per - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Per Bayesian --------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of day
probs_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_30_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_30_TOD$predict(what = "delta",
                                                            newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOD <- newdata$TOD

### Plot stationary state probs against time of day
perb <- ggplot(probs_df, aes(TOD, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  labs(x = "Time Of Day", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Per - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Gorm Frequentist -------------------------------------------------------------

gormf <- freq_17_TOD$plot(what = "delta", var = "TOD") +
  theme(legend.position="top") +
  ggtitle("Gorm - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Gorm Bayesian -------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of day
probs_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_17_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_17_TOD$predict(what = "delta",
                                                            newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOD <- newdata$TOD

### Plot stationary state probs against time of day
gormb <- ggplot(probs_df, aes(TOD, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  labs(x = "Time Of Day", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Gorm - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Ragnar Frequentist -------------------------------------------------------------

ragf <- freq_6_TOD$plot(what = "delta", var = "TOD") +
  theme(legend.position="top") +
  ggtitle("Ragnar - Frequentist") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))

#------------------------------------- Ragnar Bayesian ---------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities over
### grid of time of day
probs_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_6_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ", 1:2), each = 100),
                      prob = as.vector(bayes_6_TOD$predict(what = "delta",
                                                           newdata = newdata)))
  probs$group <- paste0("iter ", i, " - ", probs$state)
  probs_df <- rbind(probs_df, probs)
}
probs_df$TOD <- newdata$TOD

### Plot stationary state probs against time of day
ragb <- ggplot(probs_df, aes(TOD, prob, group = group, col = state)) +
  geom_line(size = 0.2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  labs(x = "Time Of Day", y = "Stationary State Probabilities", col = NULL) +
  scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
  guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
  theme_bw()  +
  theme(legend.position="top") +
  ggtitle("Ragnar - Bayesian") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
  theme(legend.text=element_text(size=17))


#------------------------------------- Collected plot ----------------------------------------------------------------

grid.arrange(perf, gormf, ragf,
             perb, gormb, ragb,
             ncol = 3)























####################################################################################################################
#------------------------------------- Transition probability matrix plot ------------------------------------------
####################################################################################################################





#------------------------------------- Per Frequentist -------------------------------------------------------------

probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
### Compute transition probabilities over
### grid of time of day
newdata <- data.frame(TOD = seq(0, 24, length = 100))


probs11 <- data.frame(prob = as.vector((freq_30_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector((freq_30_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector((freq_30_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs11_df <- rbind(probs11_df, probs11)
probs11_df$TOD <- newdata$TOD

probs12 <- data.frame(prob = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs12_df <- rbind(probs12_df, probs12)
probs12_df$TOD <- newdata$TOD

probs21 <- data.frame(prob = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector(1-(freq_30_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs21_df <- rbind(probs21_df, probs21)
probs21_df$TOD <- newdata$TOD

probs22 <- data.frame(prob = as.vector((freq_30_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector((freq_30_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector((freq_30_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs22_df <- rbind(probs22_df, probs22)
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
pf_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.99, 0.995, 1), limits = c(0.99,1)),
  ggplot(probs12_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.005, 0.01), limits = c(0, 0.01)),
  ggplot(probs21_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.005, 0.01), limits = c(0,0.01)),
  ggplot(probs22_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.99, 0.995, 1), limits = c(0.99,1)),
  ncol = 2, top = grid::textGrob("Per - Frequentist", x = 0.3, gp=grid::gpar(fontsize=23)))

#------------------------------------- Per Bayesian ----------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute trnsition probabilities over
### grid of time of day
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_30_TOD$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_30_TOD$predict(what = "tpm",
                                                               newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_30_TOD$predict(what = "tpm",
                                                                 newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_30_TOD$predict(what = "tpm",
                                                                 newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_30_TOD$predict(what = "tpm",
                                                               newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOD <- newdata$TOD
probs12_df$TOD <- newdata$TOD
probs21_df$TOD <- newdata$TOD
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
pb_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.99, 0.995, 1), limits = c(0.99,1)),
  ggplot(probs12_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.005, 0.01), limits = c(0,0.01)),
  ggplot(probs21_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.005, 0.01), limits = c(0,0.01)),
  ggplot(probs22_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.99, 0.995, 1), limits = c(0.99,1)),
  ncol = 2, top = grid::textGrob("Per - Bayesian", x = 0.3, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Gorm Frequentist ---------------------------------------------------------------

probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
### Compute transition probabilities over
### grid of time of day
newdata <- data.frame(TOD = seq(0, 24, length = 100))

probs11 <- data.frame(prob = as.vector((freq_17_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector((freq_17_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector((freq_17_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs11_df <- rbind(probs11_df, probs11)
probs11_df$TOD <- newdata$TOD

probs12 <- data.frame(prob = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs12_df <- rbind(probs12_df, probs12)
probs12_df$TOD <- newdata$TOD

probs21 <- data.frame(prob = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                              newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector(1-(freq_17_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs21_df <- rbind(probs21_df, probs21)
probs21_df$TOD <- newdata$TOD

probs22 <- data.frame(prob = as.vector((freq_17_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector((freq_17_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector((freq_17_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs22_df <- rbind(probs22_df, probs22)
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
gf_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.97, 0.985, 1), limits = c(0.97,1)),
  ggplot(probs12_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.015, 0.03), limits = c(0,0.03)),
  ggplot(probs21_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.015, 0.03), limits = c(0, 0.03)),
  ggplot(probs22_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.97, 0.985, 1), limits = c(0.97,1)),
  ncol = 2, top = grid::textGrob("Gorm - Frequentist", x = 0.3, gp=grid::gpar(fontsize=23)))

#------------------------------------- Gorm Bayesian ---------------------------------------------------------------

n_iter <- 2000
### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute transition probabilities over
### grid of time of day
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_17_TOD$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_17_TOD$predict(what = "tpm",
                                                               newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_17_TOD$predict(what = "tpm",
                                                                 newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_17_TOD$predict(what = "tpm",
                                                                 newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_17_TOD$predict(what = "tpm",
                                                               newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOD <- newdata$TOD
probs12_df$TOD <- newdata$TOD
probs21_df$TOD <- newdata$TOD
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
gb_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.97, 0.985, 1), limits = c(0.97,1)),
  ggplot(probs12_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.015, 0.03), limits = c(0,0.03)),
  ggplot(probs21_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.015, 0.03), limits = c(0, 0.03)),
  ggplot(probs22_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.97, 0.985, 1), limits = c(0.97,1)),
  ncol = 2, top = grid::textGrob("Gorm - Bayesian", x = 0.3, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Ragnar Frequentist -------------------------------------------------------------

probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
### Compute transition probabilities over
### grid of time of day
newdata <- data.frame(TOD = seq(0, 24, length = 100))

probs11 <- data.frame(prob = as.vector((freq_6_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector((freq_6_TOD$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector((freq_6_TOD$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs11_df <- rbind(probs11_df, probs11)
probs11_df$TOD <- newdata$TOD

probs12 <- data.frame(prob = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$mean[1,1,])^(1)),
                      ucl = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$ucl[1,1,])^(1)),
                      lcl = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$lcl[1,1,])^(1)))
probs12_df <- rbind(probs12_df, probs12)
probs12_df$TOD <- newdata$TOD

probs21 <- data.frame(prob = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                             newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector(1-(freq_6_TOD$predict(what = "tpm",
                                                            newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs21_df <- rbind(probs21_df, probs21)
probs21_df$TOD <- newdata$TOD

probs22 <- data.frame(prob = as.vector((freq_6_TOD$predict(what = "tpm",
                                                           newdata = newdata, n_post = 1000)$mean[2,2,])^(1)),
                      ucl = as.vector((freq_6_TOD$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$ucl[2,2,])^(1)),
                      lcl = as.vector((freq_6_TOD$predict(what = "tpm",
                                                          newdata = newdata, n_post = 1000)$lcl[2,2,])^(1)))
probs22_df <- rbind(probs22_df, probs22)
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
rf_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ggplot(probs12_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs21_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs22_df, aes(TOD, prob)) +
    geom_line(size = 1, alpha = 1, color = "#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.1, color = "steelblue", fill = "steelblue") +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ncol = 2, top = grid::textGrob("Ragnar - Frequentist", x = 0.33, gp = grid::gpar(fontsize = 23)))

#------------------------------------- Ragnar Bayesian -------------------------------------------------------------

n_iter <- 2000
# Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute transition probabilities over
### grid of time of day
probs11_df <- data.frame()
probs12_df <- data.frame()
probs21_df <- data.frame()
probs22_df <- data.frame()
newdata <- data.frame(TOD = seq(0, 24, length = 100))
for(i in ind_post) {
  bayes_6_TOD$update_par(iter = i)
  probs11 <- data.frame(prob = as.vector((bayes_6_TOD$predict(what = "tpm",
                                                              newdata = newdata)[1,1,])^(1)))
  probs11$group <- paste0("iter ", i)
  probs11_df <- rbind(probs11_df, probs11)
  
  probs12 <- data.frame(prob = as.vector(1-(bayes_6_TOD$predict(what = "tpm",
                                                                newdata = newdata)[1,1,])^(1)))
  probs12$group <- paste0("iter ", i)
  probs12_df <- rbind(probs12_df, probs12)
  
  probs21 <- data.frame(prob = as.vector(1-(bayes_6_TOD$predict(what = "tpm",
                                                                newdata = newdata)[2,2,])^(1)))
  probs21$group <- paste0("iter ", i)
  probs21_df <- rbind(probs21_df, probs21)
  
  probs22 <- data.frame(prob = as.vector((bayes_6_TOD$predict(what = "tpm",
                                                              newdata = newdata)[2,2,])^(1)))
  probs22$group <- paste0("iter ", i)
  probs22_df <- rbind(probs22_df, probs22)
}
probs11_df$TOD <- newdata$TOD
probs12_df$TOD <- newdata$TOD
probs21_df$TOD <- newdata$TOD
probs22_df$TOD <- newdata$TOD

### Plot trans probs against time of day
rb_1tpm <- grid.arrange(
  ggplot(probs11_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ggplot(probs12_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(1 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs21_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 1)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0, 0.05, 0.1), limits = c(0,0.1)),
  ggplot(probs22_df, aes(TOD, prob, group = group)) +
    geom_line(size = 0.2, alpha = 0.2, color ="#901a1E") +
    scale_x_continuous(breaks = seq(0, 24, by = 4)) +
    labs(x = "Time Of Day", y = "P(2 -> 2)", col = NULL) +
    scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
    guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1))) +
    theme_bw() +
    theme(axis.title.x = element_text(size=17)) +
    theme(axis.title.y = element_text(size=15)) +
    theme(axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=13)) +
    scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) +
    scale_y_continuous(breaks=c(0.9, 0.95, 1), limits = c(0.9,1)),
  ncol = 2, top = grid::textGrob("Ragnar - Bayesian", x = 0.33, gp = grid::gpar(fontsize = 23)))



#------------------------------------- Collected plot --------------------------------------------------------------

grid.arrange(pf_1tpm, pb_1tpm,
             gf_1tpm, gb_1tpm,
             rf_1tpm, rb_1tpm,
             ncol = 2)




