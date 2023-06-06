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
load("sims_02_06.RData")

####################################################################################################################
####################################################################################################################
####################################################################################################################













####################################################################################################################
####################################################################################################################
#------------------------------------- Model checking -------------------------------------------------------------
####################################################################################################################
####################################################################################################################



####################################################################################################################
#------------------------------------- Bayesian Convergence --------------------------------------------------------
####################################################################################################################

### Rhat and ESS values

bayes_30$out_stan()

bayes_17$out_stan()

bayes_6$out_stan()

bayes_30_TOY$out_stan()

bayes_17_TOY$out_stan()

bayes_6_TOY$out_stan()

bayes_30_TOD$out_stan()

bayes_17_TOD$out_stan()

bayes_6_TOD$out_stan()




####################################################################################################################
#------------------------------------- Traceplots and density plots ------------------------------------------------
####################################################################################################################


####################################################################################################################
#------------------------------------- Null model ------------------------------------------------------------------
####################################################################################################################




#------------------------------------- Per -------------------------------------------------------------------------

stan_trace(bayes_30$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__")) + ggtitle("Per - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_30$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Per - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))



#------------------------------------- Gorm ------------------------------------------------------------------------


stan_trace(bayes_17$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__")) + ggtitle("Gorm - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_17$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Gorm - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))


#------------------------------------- Ragnar ----------------------------------------------------------------------


stan_trace(bayes_6$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__")) + ggtitle("Ragnar - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_6$out_stan(), pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Ragnar - Null") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))




####################################################################################################################
#------------------------------------- TOY model -------------------------------------------------------------------
####################################################################################################################




#------------------------------------- Per -------------------------------------------------------------------------


stan_trace(bayes_30_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Per - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_30_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Per - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 13),
                                             title.theme = element_text(size = 17)))







#------------------------------------- Gorm ------------------------------------------------------------------------



stan_trace(bayes_17_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Gorm - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_17_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Gorm - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 13),
                                             title.theme = element_text(size = 17)))


#------------------------------------- Ragnar ----------------------------------------------------------------------


stan_trace(bayes_6_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Ragnar - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_6_TOY$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Ragnar - TOY") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 13),
                                             title.theme = element_text(size = 17)))


####################################################################################################################
#------------------------------------- TOD model -------------------------------------------------------------------
####################################################################################################################


#------------------------------------- Per -------------------------------------------------------------------------




stan_trace(bayes_30_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Per - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_30_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Per - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 13),
                                             title.theme = element_text(size = 17)))


#------------------------------------- Gorm ------------------------------------------------------------------------


stan_trace(bayes_17_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Gorm - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_17_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Gorm - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))



#------------------------------------- Ragnar ----------------------------------------------------------------------


stan_trace(bayes_6_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__")) + ggtitle("Ragnar - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 17),
                                             title.theme = element_text(size = 17)))

stan_dens(bayes_6_TOD$out_stan(),pars = c("coeff_fe_obs[1]", "coeff_fe_obs[2]", "coeff_fe_obs[3]", "coeff_fe_obs[4]", "coeff_fe_hid[1]", "coeff_fe_hid[2]", "coeff_fe_hid[3]", "coeff_fe_hid[4]", "coeff_fe_hid[5]", "coeff_fe_hid[6]", "coeff_fe_hid[7]", "coeff_fe_hid[8]", "lp__"), fill = "#901a1E", color = "steelblue") + ggtitle("Ragnar - TOD") +
  theme(legend.text = element_text(size = 17, angle = 0)) +
  theme(plot.title = element_text(size=23, hjust=0)) +
  scale_x_continuous(breaks=NULL) +
  scale_y_continuous(breaks=NULL) +
  theme(legend.position = c(0.9, 0.11)) + 
  scale_colour_discrete(guide = guide_legend(direction = "horizontal", title.position = "top",
                                             label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                             label.theme = element_text(angle = 0, size = 13),
                                             title.theme = element_text(size = 17)))






####################################################################################################################
#------------------------------------- AIC and WAIC values ---------------------------------------------------------
####################################################################################################################




AIC(freq_30,freq_17,freq_6)

AIC(freq_17_TOY,freq_6_TOY)

AIC(freq_30_TOD,freq_17_TOD,freq_6_TOD)


waic(extract_log_lik(bayes_30$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_17$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_6$out_stan(), parameter_name = "lp__"))


waic(extract_log_lik(bayes_30_TOY$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_17_TOY$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_6_TOY$out_stan(), parameter_name = "lp__"))


waic(extract_log_lik(bayes_30_TOD$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_17_TOD$out_stan(), parameter_name = "lp__"))

waic(extract_log_lik(bayes_6_TOD$out_stan(), parameter_name = "lp__"))





####################################################################################################################
#------------------------------------- Pseudo-residuals ------------------------------------------------------------
####################################################################################################################





pb_TOY_res <- bayes_30_TOY$pseudores()
gb_TOY_res <- bayes_17_TOY$pseudores()
rb_TOY_res <- bayes_6_TOY$pseudores()

pb_TOD_res <- bayes_30_TOD$pseudores()
gb_TOD_res <- bayes_17_TOD$pseudores()
rb_TOD_res <- bayes_6_TOD$pseudores()

gf_TOY_res <- freq_17_TOY$pseudores()
rf_TOY_res <- freq_6_TOY$pseudores()

pf_TOD_res <- freq_30_TOD$pseudores()
gf_TOD_res <- freq_17_TOD$pseudores()
rf_TOD_res <- freq_6_TOD$pseudores()


par(oma=c(3,3,0,0),mar=c(3,3,2,2),mfrow=c(3,3))


plot(x = 0, y = 0, col = "#901a1E", main = "Per - TOY - Frequentist", xlab = "", ylab = "", cex.main=2)
abline(a = 0, b = 1, col = "#901a1E", lwd = 5)
abline(a = 0, b = -1, col = "#901a1E", lwd = 5)


qqnorm(gf_TOY_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Gorm - TOY - Frequentist", cex.main=2)
qqline(gf_TOY_res, col = "steelblue", lwd = 3)

qqnorm(rf_TOY_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Ragnar - TOY - Frequentist", cex.main=2)
qqline(rf_TOY_res, col = "steelblue", lwd = 3)

qqnorm(pb_TOY_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Per - TOY - Bayesian", cex.main=2)
qqline(pb_TOY_res, col = "steelblue", lwd = 3)

qqnorm(gb_TOY_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Gorm - TOY - Bayesian", cex.main=2)
qqline(gb_TOY_res, col = "steelblue", lwd = 3)

qqnorm(rb_TOY_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Ragnar - TOY - Bayesian", cex.main=2)
qqline(rb_TOY_res, col = "steelblue", lwd = 3)

qqnorm(pf_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Per - TOD - Frequentist", cex.main=2)
qqline(pf_TOD_res, col = "steelblue", lwd = 3)

qqnorm(gf_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Gorm - TOD - Frequentist", cex.main=2)
qqline(gf_TOD_res, col = "steelblue", lwd = 3)

qqnorm(rf_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Ragnar - TOD - Frequentist", cex.main=2, xlab = "", ylab = "")
qqline(rf_TOD_res, col = "steelblue", lwd = 3)

qqnorm(pb_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Per - TOD - Bayesian", cex.main=2)
qqline(pb_TOD_res, col = "steelblue", lwd = 3)

qqnorm(gb_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Gorm - TOD - Bayesian", cex.main=2)
qqline(gb_TOD_res, col = "steelblue", lwd = 3)

qqnorm(rb_TOD_res, pch = 1, frame = FALSE, col = "#901a1E", main = "Ragnar - TOD - Bayesian", cex.main=2)
qqline(rb_TOD_res, col = "steelblue", lwd = 3)

mtext(text="Theoretical Quantiles",side=1,line=0,outer=TRUE, cex = 2)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2, adj = 0.7)









