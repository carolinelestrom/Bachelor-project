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
#------------------------------------- Data ------------------------------------------------------------------------
####################################################################################################################
####################################################################################################################




#------------------------------------- Per -------------------------------------------------------------------------

data30 = read.csv("Shark30-100982-Series.csv", header = TRUE)
x = paste(data30$Day, data30$Time)
data30$POSIXct = as.POSIXct(strptime(x, "%d-%b-%Y %H:%M:%S"))

M <- length(data30$Depth)
diff <- rep(NA,M)
for (i in 2:M){
  diff[i] <- data30$Depth[i-1] - data30$Depth[i]
}
data30 <- data30 %>%
  mutate(Speed = diff,
         Speed_abs = abs(diff),
         TOD = as.numeric(format(data30$POSIXct, "%H")) +
           as.numeric(format(data30$POSIXct, "%M"))/60,
         TOY = as.numeric(month(POSIXct)) + as.numeric(day(POSIXct)/31))
d30 <- data30 %>%
  dplyr::select(c("Day", "Time", "Depth", "Temperature", "POSIXct", "Speed", "Speed_abs", "TOD", "TOY"))

d30$Depth[1] <- 1




#------------------------------------- Gorm ------------------------------------------------------------------------

data17 = read.csv("Shark17-tag1-158793-Series.csv", header = TRUE)
x = paste(data17$Day, data17$Time)
data17$POSIXct = as.POSIXct(strptime(x, "%d-%b-%Y %H:%M:%S"))

M <- length(data17$Depth)
diff <- rep(NA,M)
for (i in 2:M){
  diff[i] <- data17$Depth[i-1] - data17$Depth[i]
}
data17 <- data17 %>%
  mutate(Speed = diff,
         Speed_abs = abs(diff),
         TOD = as.numeric(format(data17$POSIXct, "%H")) +
           as.numeric(format(data17$POSIXct, "%M"))/60,
         TOY = month(POSIXct) + day(POSIXct)/31)
d17 <- data17 %>%
  dplyr::select(c("Day", "Time", "Depth", "Temperature", "POSIXct", "Speed", "Speed_abs", "TOD", "TOY"))




#------------------------------------- Ragnar -----------------------------------------------------------------------

data63 = read.csv2("Shark6-tag3-138258.csv")

x = paste(data63$Day, data63$Time)
data63$POSIXct = as.POSIXct(strptime(x, "%d/%b/%Y %H.%M.%S"))

M <- length(data63$Depth)
diff <- rep(NA,M)
for (i in 2:M){
  diff[i] <- data63$Depth[i-1] - data63$Depth[i]
}

data63 <- data63 %>%
  dplyr::mutate(Speed = diff,
                Speed_abs = abs(diff),
                TOD = as.numeric(format(data63$POSIXct, "%H")) +
                  as.numeric(format(data63$POSIXct, "%M"))/60,
                TOY = month(POSIXct) + day(POSIXct)/31)
d63 <- data63 %>%
  dplyr::select(c("Day", "Time", "Depth", "Temperature", "POSIXct", "Speed", "Speed_abs", "TOD", "TOY"))







####################################################################################################################
####################################################################################################################
#------------------------------------- Timeseries plot -------------------------------------------------------------
####################################################################################################################
####################################################################################################################



grid.arrange(
  d30 %>%
    ggplot(aes(x = POSIXct, y = Depth)) +
    geom_line(color = "#901a1E") +
    theme_bw() +
    xlab("Time") +
    ggtitle("Per - Depth") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_x_continuous(label = c("Sep", "Nov", "Jan", "Mar", "May", "June")) +
    scale_y_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 1000)),
  d17 %>%
    ggplot(aes(x = POSIXct, y = Depth)) +
    geom_line(color = "#901a1E") +
    theme_bw() +
    xlab("Time") +
    ggtitle("Gorm - Depth") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7, hjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_y_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 1000)),
  d63 %>%
    ggplot(aes(x = POSIXct, y = Depth)) +
    geom_line(color = "#901a1E") +
    theme_bw() +
    xlab("Time") +
    ggtitle("Ragnar - Depth") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7, hjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_y_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 1000)),
  d30 %>%
    ggplot(aes(x = POSIXct, y = Speed_abs)) +
    geom_point(alpha = 0.5, color = "#901a1E") +
    theme_bw() +
    xlab("Time") + ylab("Speed") +
    ggtitle("Per - Speed") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_x_continuous(label = c("Sep", "Nov", "Jan", "Mar", "May", "June")) +
    scale_y_continuous(breaks=c(0, 50, 100, 150, 200), limits = c(0, 200)),
  d17 %>%
    ggplot(aes(x = POSIXct, y = Speed_abs)) +
    geom_point(alpha = 0.5, color = "#901a1E") +
    theme_bw() +
    xlab("Time") + ylab("Speed") +
    ggtitle("Gorm - Speed") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7, hjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_y_continuous(breaks=c(0, 50, 100, 150, 200), limits = c(0, 200)),
  d63 %>%
    ggplot(aes(x = POSIXct, y = Speed_abs)) +
    geom_point(alpha = 0.5, color = "#901a1E") +
    theme_bw() +
    xlab("Time") + ylab("Speed") +
    ggtitle("Ragnar - Speed") +
    theme(plot.title = element_text(size=37, hjust=0)) +
    theme(axis.title = element_text(size=27)) +
    theme(axis.text.x = element_text(size=17, angle = 45, vjust = 0.7, hjust = 0.7),
          axis.text.y = element_text(size=17)) +
    scale_y_continuous(breaks=c(0, 50, 100, 150, 200), limits = c(0, 200)),
  ncol = 3)










####################################################################################################################
####################################################################################################################
#------------------------------------- Histograms ------------------------------------------------------------------
####################################################################################################################
####################################################################################################################

d30 %>%
  ggplot(aes(x = Depth)) +
  geom_histogram(fill = "white", color = "#901a1E", binwidth = 20) +
  theme_bw() +
  xlab("Depth") + ylab("Density") +
  ggtitle("Per") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 900)) +
  scale_y_continuous(breaks=NULL)


d17 %>%
  ggplot(aes(x = Depth)) +
  geom_histogram(fill = "white", color = "#901a1E", binwidth = 20) +
  theme_bw() +
  xlab("Depth") + ylab("Density") +
  ggtitle("Gorm") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 900)) +
  scale_y_continuous(breaks=NULL) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))




d63 %>%
  ggplot(aes(x = Depth)) +
  geom_histogram(fill = "white", color = "#901a1E", binwidth = 30) +
  theme_bw() +
  xlab("Depth") + ylab("Density") +
  ggtitle("Ragnar") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=20)) +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000), limits = c(0, 900)) +
  scale_y_continuous(breaks=NULL) +
  theme(axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))











####################################################################################################################
####################################################################################################################
#------------------------------------- Mixture models for starting values ------------------------------------------
####################################################################################################################
####################################################################################################################

#------------------------------------- Per -------------------------------------------------------------------------
fit30_Depth <- mixfit(na.omit(d30$Depth), ncomp = 2, family = "gamma", max_iter = 1000)
fit30_Depth


#------------------------------------- Gorm ------------------------------------------------------------------------

fit17_Depth <- mixfit(d17$Depth, ncomp = 2, family = "gamma")
fit17_Depth




#------------------------------------- Ragnar -----------------------------------------------------------------------

fit6_Depth <- mixfit(na.omit(d63$Depth), ncomp = 2, family = "gamma")
fit6_Depth














####################################################################################################################
####################################################################################################################
#------------------------------------- Fitting the models ----------------------------------------------------------
####################################################################################################################
####################################################################################################################




####################################################################################################################
#------------------------------------- Initial settings ------------------------------------------------------------
####################################################################################################################


#------------------------------------- Initial values --------------------------------------------------------------
### Depth
mean30_Depth <- c(220.7267738, 596.1396368)
sd30_Depth <- c(55.2805743, 63.6978207)


mean17_Depth <- c(341.1282154, 433.7704086)
sd17_Depth <- c(92.4953976, 36.6602395)


mean6_Depth <- c(326.9257489, 735.8588410)
sd6_Depth <- c(97.0872998, 116.3449541)













####################################################################################################################
#------------------------------------- Frequentist fit -------------------------------------------------------------
####################################################################################################################


####################################################################################################################
#------------------------------------- Null model ------------------------------------------------------------------
####################################################################################################################



### Distribution
dists <- list(Depth = "gamma2")


#------------------------------------- Per -------------------------------------------------------------------------

### Hidden model
hid_30_freq <- MarkovChain$new(data = d30, n_states = 2, initial_state = "stationary",
                               formula = ~1)

### Parameters
par0_30_freq <- list(Depth = list(mean = c(mean30_Depth[1], mean30_Depth[2]),
                                  sd = c(sd30_Depth[1], sd30_Depth[2])))

### Observation model
obs_30_freq <- Observation$new(data = d30, dists = dists,
                               n_states = 2, par = par0_30_freq)

### Define HMM
freq_30 <- HMM$new(obs = obs_30_freq, hid = hid_30_freq)

freq_30$fit(silent = TRUE)




#------------------------------------- Gorm ------------------------------------------------------------------------

### Hidden model
hid_17_freq <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                               formula = ~1)

### Parameters
par0_17_freq <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                                  sd = c(sd17_Depth[1], sd17_Depth[2])))

### Observation model
obs_17_freq <- Observation$new(data = d17, dists = dists,
                               n_states = 2, par = par0_17_freq)

### Define HMM
freq_17 <- HMM$new(obs = obs_17_freq, hid = hid_17_freq)

freq_17$fit(silent = TRUE)


#------------------------------------- Ragnar -----------------------------------------------------------------------


### Hidden model
hid_6_freq <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                              formula = ~1)

### Parameters
par0_6_freq <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                                 sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6_freq <- Observation$new(data = d63, dists = dists,
                              n_states = 2, par = par0_6_freq)

### Define HMM
freq_6 <- HMM$new(obs = obs_6_freq, hid = hid_6_freq)

freq_6$fit(silent = TRUE)



####################################################################################################################
#------------------------------------- TOY model -------------------------------------------------------------------
####################################################################################################################




### Distribution
dists <- list(Depth = "gamma2")

### formula for transition probabilities
form <- ~mSpline(TOY, df = 3, periodic = TRUE, Boundary.knots = c(1,13))


#------------------------------------- Gorm ------------------------------------------------------------------------

### Hidden model
hid_17_TOY_freq <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                                   formula = form)

### Parameters
par0_17_TOY_freq <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                                      sd = c(sd17_Depth[1], sd17_Depth[2])))


### Observation model
obs_17_TOY_freq <- Observation$new(data = d17, dists = dists,
                                   n_states = 2, par = par0_17_TOY_freq)

### Define HMM
freq_17_TOY <- HMM$new(obs = obs_17_TOY_freq, hid = hid_17_TOY_freq)

freq_17_TOY$fit(silent = TRUE)

#------------------------------------- Ragnar -----------------------------------------------------------------------

### Hidden model
hid_6_TOY_freq <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                                  formula = form)

### Parameters
par0_6_TOY_freq <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                                     sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6_TOY_freq <- Observation$new(data = d63, dists = dists,
                                  n_states = 2, par = par0_6_TOY_freq)

### Define HMM
freq_6_TOY <- HMM$new(obs = obs_6_TOY_freq, hid = hid_6_TOY_freq)

freq_6_TOY$fit(silent = TRUE)




####################################################################################################################
#------------------------------------- TOD model -------------------------------------------------------------------
####################################################################################################################




### Distribution
dists <- list(Depth = "gamma2")

### formula for transition probabilities
form <- ~mSpline(TOD, df = 3, periodic = TRUE, Boundary.knots = c(0,24))


#------------------------------------- Per -------------------------------------------------------------------------

### Hidden model
hid_30_TOD_freq <- MarkovChain$new(data = d30, n_states = 2, initial_state = "stationary",
                                   formula = form)


### Parameters
par0_30_TOD_freq <- list(Depth = list(mean = c(mean30_Depth[1], mean30_Depth[2]),
                                      sd = c(sd30_Depth[1], sd30_Depth[2])))

### Observation model
obs_30_TOD_freq <- Observation$new(data = d30, dists = dists,
                                   n_states = 2, par = par0_30_TOD_freq)


### Define HMM
freq_30_TOD <- HMM$new(obs = obs_30_TOD_freq, hid = hid_30_TOD_freq)

freq_30_TOD$fit(silent = TRUE)




#------------------------------------- Gorm ------------------------------------------------------------------------

### Hidden model
hid_17_TOD_freq <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                                   formula = form)


### Parameters
par0_17_TOD_freq <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                                      sd = c(sd17_Depth[1], sd17_Depth[2])))

### Observation model
obs_17_TOD_freq <- Observation$new(data = d17, dists = dists,
                                   n_states = 2, par = par0_17_TOD_freq)


### Define HMM
freq_17_TOD <- HMM$new(obs = obs_17_TOD_freq, hid = hid_17_TOD_freq)

freq_17_TOD$fit(silent = TRUE)


#------------------------------------- Ragnar -----------------------------------------------------------------------

### Hidden model
hid_6_TOD_freq <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                                   formula = form)


### Parameters
par0_6_TOD_freq <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                                      sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6_TOD_freq <- Observation$new(data = d63, dists = dists,
                                   n_states = 2, par = par0_6_TOD_freq)


### Define HMM
freq_6_TOD <- HMM$new(obs = obs_6_TOD_freq, hid = hid_6_TOD_freq)

freq_6_TOD$fit(silent = TRUE)










####################################################################################################################
#------------------------------------- Bayesian fit ----------------------------------------------------------------
####################################################################################################################



####################################################################################################################
#------------------------------------- Null model ------------------------------------------------------------------
####################################################################################################################



### Distribution
dists <- list(Depth = "gamma2")

#------------------------------------- Per -------------------------------------------------------------------------


### Hidden model
hid_30 <- MarkovChain$new(data = d30, n_states = 2, initial_state = "stationary",
                          formula = ~1)

### Parameters
par0_30 <- list(Depth = list(mean = c(mean30_Depth[1], mean30_Depth[2]),
                             sd = c(sd30_Depth[1], sd30_Depth[2])))

### Observation model
obs_30 <- Observation$new(data = d30, dists = dists,
                          n_states = 2, par = par0_30)


### Define HMM
bayes_30 <- HMM$new(obs = obs_30, hid = hid_30)

### Priors
prior_obs <- matrix(c(NA, NA,
                      NA, NA,
                      NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)
prior_hid <- matrix(c(NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)

### Set priors
bayes_30$set_priors(list(coeff_fe_obs = prior_obs,
                         coeff_fe_hid = prior_hid))

n_chain <- 4
n_iter <- 2000
bayes_30$fit_stan(chains = n_chain, iter = n_iter)




#------------------------------------- Gorm ------------------------------------------------------------------------


### Hidden model
hid_17 <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                          formula = ~1)

### Parameters
par0_17 <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                             sd = c(sd17_Depth[1], sd17_Depth[2])))

### Observation model
obs_17 <- Observation$new(data = d17, dists = dists,
                          n_states = 2, par = par0_17)

### Define HMM
bayes_17 <- HMM$new(obs = obs_17, hid = hid_17)

### Priors
prior_obs <- matrix(c(NA, NA,
                      NA, NA,
                      NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)
prior_hid <- matrix(c(NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)

### Set priors
bayes_17$set_priors(list(coeff_fe_obs = prior_obs,
                         coeff_fe_hid = prior_hid))

n_chain <- 4
n_iter <- 2000
bayes_17$fit_stan(chains = n_chain, iter = n_iter)





#------------------------------------- Ragnar -----------------------------------------------------------------------



### Hidden model
hid_6 <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                         formula = ~1)

### Parameters
par0_6 <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                            sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6 <- Observation$new(data = d63, dists = dists,
                         n_states = 2, par = par0_6)

### Define HMM
bayes_6 <- HMM$new(obs = obs_6, hid = hid_6)

### Priors
prior_obs <- matrix(c(NA, NA,
                      NA, NA,
                      NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)
prior_hid <- matrix(c(NA, NA,
                      NA, NA),
                    ncol = 2, byrow = TRUE)

### Set priors
bayes_6$set_priors(list(coeff_fe_obs = prior_obs,
                        coeff_fe_hid = prior_hid))

n_chain <- 4
n_iter <- 2000
bayes_6$fit_stan(chains = n_chain, iter = n_iter)





####################################################################################################################
#------------------------------------- TOY model -------------------------------------------------------------------
####################################################################################################################




### Distribution
dists <- list(Depth = "gamma2")

### formula for transition probabilities
form <- ~mSpline(TOY, df = 3, periodic = TRUE, Boundary.knots = c(1,13))

#------------------------------------- Per -------------------------------------------------------------------------

### Hidden model
hid_30_TOY <- MarkovChain$new(data = d30, n_states = 2, initial_state = "stationary",
                              formula = form)

### Parameters
par0_30_TOY <- list(Depth = list(mean = c(mean30_Depth[1], mean30_Depth[2]),
                                 sd = c(sd30_Depth[1], sd30_Depth[2])))

### Observation model
obs_30_TOY <- Observation$new(data = d30, dists = dists,
                              n_states = 2, par = par0_TOY)

### Define HMM
bayes_30_TOY <- HMM$new(obs = obs_30_TOY, hid = hid_30_TOY)

### Set priors
bayes_30_TOY$set_priors()

n_chain <- 4
n_iter <- 200
bayes_30_TOY$fit_stan(chains = n_chain, iter = n_iter)



#------------------------------------- Gorm ------------------------------------------------------------------------

### Hidden model
hid_17_TOY <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                              formula = form)

### Parameters
par0_17_TOY <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                                 sd = c(sd17_Depth[1], sd17_Depth[2])))

### Observation model
obs_17_TOY <- Observation$new(data = d17, dists = dists,
                              n_states = 2, par = par0_17_TOY)

### Define HMM
bayes_17_TOY <- HMM$new(obs = obs_17_TOY, hid = hid_17_TOY)

### Set priors
bayes_17_TOY$set_priors()



n_chain <- 4
n_iter <- 2000

bayes_17_TOY$fit_stan(chains = n_chain, iter = n_iter)



#------------------------------------- Ragnar -----------------------------------------------------------------------

### Hidden model
hid_6_TOY <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                             formula = form)

### Parameters
par0_6_TOY <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                                sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6_TOY <- Observation$new(data = d63, dists = dists,
                             n_states = 2, par = par0_6_TOY)

### Define HMM
bayes_6_TOY <- HMM$new(obs = obs_6_TOY, hid = hid_6_TOY)


### Set priors
bayes_6_TOY$set_priors()

n_chain <- 4
n_iter <- 2000

bayes_6_TOY$fit_stan(chains = n_chain, iter = n_iter)






####################################################################################################################
#------------------------------------- TOD model -------------------------------------------------------------------
####################################################################################################################






### Distribution
dists <- list(Depth = "gamma2")


### formula for transition probabilities
form <- ~mSpline(TOD, df = 3, periodic = TRUE, Boundary.knots = c(0,24))


#------------------------------------- Per -------------------------------------------------------------------------



### Hidden model
hid_30_TOD <- MarkovChain$new(data = d30, n_states = 2, initial_state = "stationary",
                              formula = form)


### Parameters
par0_30_TOD <- list(Depth = list(mean = c(mean30_Depth[1], mean30_Depth[2]),
                                 sd = c(sd30_Depth[1], sd30_Depth[2])))

### Observation model
obs_30_TOD <- Observation$new(data = d30, dists = dists,
                              n_states = 2, par = par0_30_TOD)

### Define HMM
bayes_30_TOD <- HMM$new(obs = obs_30_TOD, hid = hid_30_TOD)

### Set priors
bayes_30_TOD$set_priors()

n_chain <- 4
n_iter <- 2000
bayes_30_TOD$fit_stan(chains = n_chain, iter = n_iter)


#------------------------------------- Gorm ------------------------------------------------------------------------


### Hidden model
hid_17_TOD <- MarkovChain$new(data = d17, n_states = 2, initial_state = "stationary",
                              formula = form)

### Parameters
par0_17_TOD <- list(Depth = list(mean = c(mean17_Depth[1], mean17_Depth[2]),
                                 sd = c(sd17_Depth[1], sd17_Depth[2])))

### Observation model
obs_17_TOD <- Observation$new(data = d17, dists = dists,
                              n_states = 2, par = par0_17_TOD)

### Define HMM
bayes_17_TOD <- HMM$new(obs = obs_17_TOD, hid = hid_17_TOD)


### Set priors
bayes_17_TOD$set_priors()

n_chain <- 4
n_iter <- 2000
bayes_17_TOD$fit_stan(chains = n_chain, iter = n_iter)


#------------------------------------- Ragnar -----------------------------------------------------------------------


### Hidden model
hid_6_TOD <- MarkovChain$new(data = d63, n_states = 2, initial_state = "stationary",
                             formula = form)

### Parameters
par0_6_TOD <- list(Depth = list(mean = c(mean6_Depth[1], mean6_Depth[2]),
                                sd = c(sd6_Depth[1], sd6_Depth[2])))

### Observation model
obs_6_TOD <- Observation$new(data = d63, dists = dists,
                             n_states = 2, par = par0_6_TOD)

### Define HMM
bayes_6_TOD <- HMM$new(obs = obs_6_TOD, hid = hid_6_TOD)


### Set priors
bayes_6_TOD$set_priors()

n_chain <- 4
n_iter <- 2000
bayes_6_TOD$fit_stan(chains = n_chain, iter = n_iter)



####################################################################################################################
####################################################################################################################
####################################################################################################################


### Access to fitted models and posterior samples saved under names specified as above
load("sims_02_06.RData")

