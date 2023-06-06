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
####################################################################################################################
#------------------------------------- Extract parameter estimates -------------------------------------------------
####################################################################################################################
####################################################################################################################




####################################################################################################################
#------------------------------------- Frequentist TOY -------------------------------------------------------------
####################################################################################################################




#------------------------------------- Gorm ------------------------------------------------------------------------


### State-dependent parameters with confidence intervals

freq_17_TOY$predict(what = "obspar", n_post = 1000)

### Dataframe for different times of year
newdata <- data.frame(TOY = c(7.5, 8.5, 9.5))

### Stationary state probabilities at different times of day with confidence intervals

freq_17_TOY$predict(what = "delta", newdata = newdata, n_post = 1000)

### Dataframe for different times of year
newdata <- data.frame(TOY = c(8, 9.75))

### Transition probability matrix at different times of day with confidence intervals

freq_17_TOY$predict(what = "tpm", n_post = 1000, newdata = newdata)

#------------------------------------- Ragnar -----------------------------------------------------------------------


### State-dependent parameters with confidence intervals

freq_6_TOY$predict(what = "obspar", n_post = 1000)

### Dataframe for different times of year
newdata <- data.frame(TOY = c(5.5, 6.5, 7.5, 8.5, 9.5, 10.5))

### Stationary state probabilities at different times of day with confidence intervals

freq_6_TOY$predict(what = "delta", newdata = newdata, n_post = 1000)

### Dataframe for different times of year
newdata <- data.frame(TOY = c(6, 8, 9.75))

### Transition probability matrix at different times of day with confidence intervals

freq_6_TOY$predict(what = "tpm", n_post = 1000, newdata = newdata)




####################################################################################################################
#------------------------------------- Frequentist TOD -------------------------------------------------------------
####################################################################################################################



#------------------------------------- Per --------------------------------------------------------------------------


#------------ State-dependent parameters with confidence intervals --------------------------------------------------

freq_30_TOD$predict(what = "obspar", n_post = 1000)

#------------- Dataframe for different times of day -----------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#--------- Stationary state probabilities at different times of day with confidence intervals -----------------------

freq_30_TOD$predict(what = "delta", n_post = 1000, newdata = newdata)


#--------- Transition probability matrix at different times of day with confidence intervals ------------------------

freq_30_TOD$predict(what = "tpm", n_post = 1000, newdata = newdata)


#------------------------------------- Gorm ------------------------------------------------------------------------

#------------- State-dependent parameters with confidence intervals ------------------------------------------------

freq_17_TOD$predict(what = "obspar", n_post = 1000)

#----------------- Dataframe for different times of day ------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#----------- Stationary state probabilities at different times of day with confidence intervals --------------------

freq_17_TOD$predict(what = "delta", n_post = 1000, newdata = newdata)


#----------- Transition probability matrix at different times of day with confidence intervals ---------------------

freq_17_TOD$predict(what = "tpm", n_post = 1000, newdata = newdata)




#------------------------------------- Ragnar -----------------------------------------------------------------------


#----------- State-dependent parameters with confidence intervals ---------------------------------------------------

freq_6_TOD$predict(what = "obspar", n_post = 1000)

#------------- Dataframe for different times of day -----------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#-------- Stationary state probabilities at different times of day with confidence intervals ------------------------

freq_6_TOD$predict(what = "delta", n_post = 1000, newdata = newdata)


#------ Transition probability matrix at different times of day with confidence intervals ---------------------------

freq_6_TOD$predict(what = "tpm", n_post = 1000, newdata = newdata)



####################################################################################################################
#------------------------------------- Bayesian TOY ----------------------------------------------------------------
####################################################################################################################






#------------------------------------- Per --------------------------------------------------------------------------

#----------- State-dependent parameters with confidence intervals ---------------------------------------------------

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_30_TOY$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))

#--------------- Dataframe for different times of year -----------------------------------------------------------------
newdata <- data.frame(TOY = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 10.5, 11.5, 12.5))

#---------- Stationary state probabilities at different times of day with confidence intervals -------------------------

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - Jan", "1 - Feb", "1 - Mar", "1 - Apr", "1 - May", "1 - June", "1 - Sep", "1 - Oct", "1 - Nov", "1 - Dec",
                                           "2 - Jan", "2 - Feb", "2 - Mar", "2 - Apr", "2 - May", "2 - June", "2 - Sep", "2 - Oct", "2 - Nov", "2 - Dec")),
                                  each = 1),
                      prob = as.vector(bayes_30_TOY$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))

#------------- Dataframe for different times of year -------------------------------------------------------------
newdata <- data.frame(TOY = c(2, 4, 6, 9.75, 12))

#--------- Transition probability matrix at different times of day with confidence intervals ---------------------

### Select 100 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 100))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOY = 2", "2->1 - TOY = 2", "1->2 - TOY = 2", "2->2 - TOY = 2",
                                           "1->1 - TOY = 4", "2->1 - TOY = 4", "1->2 - TOY = 4", "2->2 - TOY = 4",
                                           "1->1 - TOY = 6", "2->1 - TOY = 6", "1->2 - TOY = 6", "2->2 - TOY = 6",
                                           "1->1 - TOY = 9.75", "2->1 - TOY = 9.75", "1->2 - TOY = 9.75", "2->2 - TOY = 9.75",
                                           "1->1 - TOY = 12", "2->1 - TOY = 12", "1->2 - TOY = 12", "2->2 - TOY = 12")),
                                  each = 1),
                      prob = as.vector(bayes_30_TOY$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))

#------------------------------------- Gorm ------------------------------------------------------------------------

#------------ State-dependent parameters with confidence intervals -------------------------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_17_TOY$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))

#------------ Dataframe for different times of year ------------------------------------------------------------------
newdata <- data.frame(TOY = c(7.5, 8.5, 9.5))

#--------- Stationary state probabilities at different times of day with confidence intervals ------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - July", "1 - Aug", "1 - Sep",
                                           "2 - July", "2 - Aug", "2 - Sep")),
                                  each = 1),
                      prob = as.vector(bayes_17_TOY$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))

#----------- Dataframe for different times of year ---------------------------------------------------------------------
newdata <- data.frame(TOY = c(8, 9.75))

#------ Transition probability matrix at different times of day with confidence intervals ------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOY = 8", "2->1 - TOY = 8", "1->2 - TOY = 8", "2->2 - TOY = 8",
                                           "1->1 - TOY = 9.75", "2->1 - TOY = 9.75", "1->2 - TOY = 9.75", "2->2 - TOY = 9.75")),
                                  each = 1),
                      prob = as.vector(bayes_17_TOY$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))


#------------------------------------- Ragnar -----------------------------------------------------------------------

#----------- State-dependent parameters with confidence intervals ---------------------------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_6_TOY$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))

#--------------- Dataframe for different times of year ---------------------------------------------------------------
newdata <- data.frame(TOY = c(5.5, 6.5, 7.5, 8.5, 9.5, 10.5))

#------- Stationary state probabilities at different times of day with confidence intervals --------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - May", "1 - June", "1 - July", "1 - Aug", "1 - Sep", "1 - Oct",
                                           "2 - May", "2 - June", "2 - July", "2 - Aug", "2 - Sep", "2 - Oct")),
                                  each = 1),
                      prob = as.vector(bayes_6_TOY$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))

#------------- Dataframe for different times of year ------------------------------------------------------------------
newdata <- data.frame(TOY = c(6, 8, 9.75))

#------- Transition probability matrix at different times of day with confidence intervals ----------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOY$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOY = 6", "2->1 - TOY = 6", "1->2 - TOY = 6", "2->2 - TOY = 6",
                                           "1->1 - TOY = 8", "2->1 - TOY = 8", "1->2 - TOY = 8", "2->2 - TOY = 8",
                                           "1->1 - TOY = 9.75", "2->1 - TOY = 9.75", "1->2 - TOY = 9.75", "2->2 - TOY = 9.75")),
                                  each = 1),
                      prob = as.vector(bayes_6_TOY$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))









####################################################################################################################
#------------------------------------- Bayesian TOD ----------------------------------------------------------------
####################################################################################################################




#------------------------------------- Per -------------------------------------------------------------------------

#-------------- State-dependent parameters with confidence intervals -----------------------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_30_TOD$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))


#--------------- Dataframe for different times of day -----------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#-------- Stationary state probabilities at different times of day with confidence intervals --------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - TOD = 0", "1 - TOD = 6", "1 - TOD = 12", "1 - TOD = 18", "1 - TOD = 24",
                                           "2 - TOD = 0", "2 - TOD = 6", "2 - TOD = 12", "2 - TOD = 18", "2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_30_TOD$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))


#--------- Transition probability matrix at different times of day with confidence intervals -------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_30_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOD = 0", "1->2 - TOD = 0", "2->1 - TOD = 0", "2->2 - TOD = 0",
                                           "1->1 - TOD = 6", "1->2 - TOD = 6", "2->1 - TOD = 6", "2->2 - TOD = 6",
                                           "1->1 - TOD = 12", "1->2 - TOD = 12", "2->1 - TOD = 12", "2->2 - TOD = 12",
                                           "1->1 - TOD = 18", "1->2 - TOD = 18", "2->1 - TOD = 18", "2->2 - TOD = 18",
                                           "1->1 - TOD = 24", "1->2 - TOD = 24", "2->1 - TOD = 24", "2->2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_30_TOD$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))


#------------------------------------- Gorm ------------------------------------------------------------------------

#------------------- State-dependent parameters with confidence intervals ------------------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_17_TOD$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))

#-------------------- Dataframe for different times of day ---------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#--------- Stationary state probabilities at different times of day with confidence intervals ----------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - TOD = 0", "1 - TOD = 6", "1 - TOD = 12", "1 - TOD = 18", "1 - TOD = 24",
                                           "2 - TOD = 0", "2 - TOD = 6", "2 - TOD = 12", "2 - TOD = 18", "2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_17_TOD$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))




#-------------- Transition probability matrix at different times of day with confidence intervals ------------------------


### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_17_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOD = 0", "1->2 - TOD = 0", "2->1 - TOD = 0", "2->2 - TOD = 0",
                                           "1->1 - TOD = 6", "1->2 - TOD = 6", "2->1 - TOD = 6", "2->2 - TOD = 6",
                                           "1->1 - TOD = 12", "1->2 - TOD = 12", "2->1 - TOD = 12", "2->2 - TOD = 12",
                                           "1->1 - TOD = 18", "1->2 - TOD = 18", "2->1 - TOD = 18", "2->2 - TOD = 18",
                                           "1->1 - TOD = 24", "1->2 - TOD = 24", "2->1 - TOD = 24", "2->2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_17_TOD$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))


#------------------------------------- Ragnar -----------------------------------------------------------------------


#-------------------- State-dependent parameters with confidence intervals ------------------------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute obspar
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0(c("mean1", "sd1", "mean2", "sd2")), each = 1),
                      prob = as.vector(bayes_6_TOD$predict(what = "obspar")))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_obs = mean(prob),
            ucl_obs = quantile(prob, probs = 0.975),
            lcl_obs = quantile(prob, probs = 0.025))

#------------------- Dataframe for different times of day ---------------------------------------------------------------
newdata <- data.frame(TOD = seq(0, 24, length = 5))

#----------- Stationary state probabilities at different times of day with confidence intervals -------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute stationary state probabilities
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1 - TOD = 0", "1 - TOD = 6", "1 - TOD = 12", "1 - TOD = 18", "1 - TOD = 24",
                                           "2 - TOD = 0", "2 - TOD = 6", "2 - TOD = 12", "2 - TOD = 18", "2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_6_TOD$predict(what = "delta", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_delta = mean(prob),
            ucl_delta = quantile(prob, probs = 0.975),
            lcl_delta = quantile(prob, probs = 0.025))


#-------------- Transition probability matrix at different times of day with confidence intervals --------------------------

### Select 1000 random posterior samples
ind_post <- sort(sample(1:(n_iter*n_chain/2), size = 1000))
### For each posterior sample, compute tpm
probs_df <- data.frame()
for(i in ind_post) {
  bayes_6_TOD$update_par(iter = i)
  probs <- data.frame(state = rep(paste0("state ",
                                         c("1->1 - TOD = 0", "1->2 - TOD = 0", "2->1 - TOD = 0", "2->2 - TOD = 0",
                                           "1->1 - TOD = 6", "1->2 - TOD = 6", "2->1 - TOD = 6", "2->2 - TOD = 6",
                                           "1->1 - TOD = 12", "1->2 - TOD = 12", "2->1 - TOD = 12", "2->2 - TOD = 12",
                                           "1->1 - TOD = 18", "1->2 - TOD = 18", "2->1 - TOD = 18", "2->2 - TOD = 18",
                                           "1->1 - TOD = 24", "1->2 - TOD = 24", "2->1 - TOD = 24", "2->2 - TOD = 24")),
                                  each = 1),
                      prob = as.vector(bayes_6_TOD$predict(what = "tpm", newdata = newdata)))
  probs_df <- rbind(probs_df, probs)
}

group_by(probs_df, state) %>%
  summarize(mid_tpm = mean(prob),
            ucl_tpm = quantile(prob, probs = 0.975),
            lcl_tpm = quantile(prob, probs = 0.025))






