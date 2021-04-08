require(ggplot2)
require(ggpubr)
setwd("")

area_dmat_italy <- read.csv("Italy_Distances.csv")[,-1] # load distances between study sites for Italy
area_dmat_ukr <- read.csv("Ukraine_Distances.csv")[,-1] # load distances between study sites for Ukraine
time_dmat_italy <- time_dmat_italy <- matrix(c( # create a matrix of time between primary periods for Italy
  # P=1,
  0, 3, 6, 12, 15,
  # P=2,
  3, 0, 3, 9, 12,
  # P=3,
  6, 3, 0, 6, 9,
  # P=4,
  12, 9, 6, 0, 3,
  #P=5
  15, 12, 9, 3, 0
),
nrow = 5, ncol=5, byrow = TRUE
)

time_dmat_ukr <- matrix(c( # create a matrix of time between primary periods for Ukraine
  # P=1,
  0, 3, 6, 12, 15,
  # P=2,
  3, 0, 3, 9, 12,
  # P=3,
  6, 3, 0, 6, 9,
  # P=4,
  12, 9, 6, 0, 3,
  #P=5
  15, 12, 9, 3, 0
),
nrow = 5, ncol=5, byrow = TRUE
)

logit <- function(x) log(x / (1-x)) # function for determining results on logit scale

inv_logit <- function(x) exp(x) / (1 + exp(x)) # function for taking the inverse logit (e.g. to determine results on probability scale)
HDI <- function(x) coda::HPDinterval(coda::as.mcmc(x), prob=0.95)  # function for extracting the highest density intervals

ps_simple <- read.csv("mcmc-samples-simple.csv") # read in model results

###Extracting parameters####

# mean and highest density intervals (HDI) of detection (mu_lodd_delta), recruitment (mu_lodd_psi = nuisance parameter)
# and removal probabilities (mu_lodd_phi) for _italy = italy, and _ukr = Ukraine)
#Log odds scale
mean(ps_simple$mu_lodd_delta_italy)
HDI(ps_simple$mu_lodd_delta_italy)
#probability scale
mean(inv_logit(ps_simple$mu_lodd_delta_italy))
HDI(inv_logit(ps_simple$mu_lodd_delta_italy))

#Superpopulation
mean(ps_simple$superpop_italy.1) # estimate the mean of the superpopulation for _italy or _ukr for study sites 1-4 (.1 to .4)
HDI(ps_simple$superpop_italy.1)

#N_hat
mean(ps_simple$N_hat_italy.1.1) # estimate the mean of the population size (N_hat) for primary periods 1 to 5 and study sites 1 to 4 (.primary_period.study_site)
HDI(ps_simple$N_hat_italy.1.1)

# Between country comparisons for detection probability (delta), survival (phi) and recruitment (psi)
mean( ps_simple$mu_lodd_phi_italy - ps_simple$mu_lodd_phi_italy ) # on log odds scale
HDI( ps_simple$mu_lodd_phi_italy - ps_simple$mu_lodd_phi_ukr )
mean( (inv_logit(ps_simple$mu_lodd_phi_italy)) - (inv_logit(ps_simple$mu_lodd_phi_ukr))) # on probability scale
HDI( (inv_logit(ps_simple$mu_lodd_phi_italy)) - (inv_logit(ps_simple$mu_lodd_phi_ukr)))
mean(exp ( (ps_simple$mu_lodd_delta_italy) - ((ps_simple$mu_lodd_delta_ukr)))) # as odds
HDI(exp ( (ps_simple$mu_lodd_delta_italy) - ((ps_simple$mu_lodd_delta_ukr))))

#Sigma dogs [study site, 1recruitment/2survival/3detection]
mean( ps_simple$sigma_dog_ukr.1.3 ) # random effect of dog per study site on recruitment/survival/detection parameters
HDI( ps_simple$sigma_dog_ukr.1.3 )

###Explanatory variables ####

###log odds of weekend effect (beta_weekend) on detection probability (delta) (anchored on weekend) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_weekend_ukr)
HDI(ps_simple$beta_weekend_ukr)
#Probability of detecting a dog on the weekday
mean(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
#Weekend probability
mean(inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr))

###log odds of beta_market effect on detection (delta) (anchored on non-market day) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_market_ukr)
HDI(ps_simple$beta_market_ukr)
#Probability of detecting a dog on a market day
mean(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))

###log odds of rain (beta_rain) effect on detection (delta) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_rain_italy)
HDI(ps_simple$beta_rain_italy)

###log odds of temperature (beta_temp) effect on detection (delta) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_temp_italy)
HDI(ps_simple$beta_temp_italy)

###log odds of sex effect on survival (.1 phi) (female -> male) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_sex_ukr.1)
HDI(ps_simple$beta_sex_ukr.1)
# total probability females' survival 
mean(inv_logit(ps_simple$mu_lodd_phi_ukr))
HDI(inv_logit(ps_simple$mu_lodd_phi_ukr))
# total probability of males' survival
mean(inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))
HDI(inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))

###log odds of sex effect on detection (.2 delta) (female -> male)####
mean(ps_simple$beta_sex_italy.2)
HDI(ps_simple$beta_sex_italy.2)

###Estimates for primary periods####
#Log odds of detection (delta) or survival (phi) for Italy/Ukraine for primary periods x.2 = phi, x.3 = delta
UAPP2phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.2.2
UAPP3phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.3.2
UAPP4phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.4.2
UAPP5phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.5.2
mean(UAPP2phi)
HDI(UAPP2phi)
#probability scale
mean(inv_logit(UAPP2phi))
HDI(inv_logit(UAPP2phi))


###Estimates for areas points####
#Log odds of detection (delta) or survival (phi) for Italy/Ukraine for study areas x.2 = phi, x.3 = delta
UASA1phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.1.2
UASA2phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.2.2
UASA3phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.3.2
UASA4phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.4.2
mean(UASA4phi)
HDI(UASA4phi)
#probability
mean(inv_logit(UASA4phi))
HDI(inv_logit(UASA4phi))

### Entry probability primary period.study area ####
mean(ps_simple$entry_prob_italy.2.1)
HDI(ps_simple$entry_prob_italy.2.1)
mean(inv_logit(ps_simple$entry_prob_italy.2.1))
HDI(inv_logit(ps_simple$entry_prob_italy.2.1))
