require(ggplot2)
require(ggpubr)
setwd("")

area_dmat_italy <- read.csv("italy_Distances.csv")[,-1] # load distances between study sites for study sites in Italy
area_dmat_ukr <- read.csv("Ukraine_Distances.csv")[,-1] # load distances between study sites for study sites in Ukraine
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

ps_simple <- read.csv("lauren-mcmc-samples-simple-19-12-2019.csv") # read in model results

###Extracting parameters####

#Log odds of xxx
mean(ps_simple$mu_lodd_delta_italy)
HDI(ps_simple$mu_lodd_delta_italy)
#probability mu_psi
mean(inv_logit(ps_simple$mu_lodd_delta_italy))
HDI(inv_logit(ps_simple$mu_lodd_delta_italy))

inv_logit(-1.10)

#Superpopulation
mean(ps_simple$superpop_italy.1) # estimate the mean of the superpopulation for _italy or _ukraine for study sites 1-4 (.1 to .4)
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

hist((ps_simple$mu_lodd_delta_italy) - (ps_simple$mu_lodd_delta_ukr))

#Sigma dogs [study site, 1recruitment/2survival/3detection]
mean( ps_simple$sigma_dog_ukr.4.3 ) # random effect of dog per study site on recruitment/survival/detection parameters
HDI( ps_simple$sigma_dog_ukr.4.3 )

###Explanatory variables ####

###log odds of weekend effect (beta_weekend) on detection probability (delta) (anchored on weekend) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_weekend_ukr)
HDI(ps_simple$beta_weekend_ukr)

#Converting the weekend effect on detection to the odds scale
mean(exp(ps_simple$beta_weekend_italy))
HDI(exp(ps_simple$beta_weekend_italy))
#Converting the weekend effect on detection to the probability scale
mean(inv_logit(ps_simple$beta_weekend_italy))
HDI(inv_logit(ps_simple$beta_weekend_italy))
#Calculation probability of detecting a dog on the weekday
mean(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
#Weekend probability
mean(inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr))
#Calcuating the odds of detecting a dog on the weekday
mean(exp(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
HDI(exp(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr))
#Odds of weekend detection
mean(exp(ps_simple$mu_lodd_delta_ukr))
HDI(exp(ps_simple$mu_lodd_delta_ukr))
#Difference in means (probability scale)
mean( (inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr)) - inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI( (inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr)) - inv_logit(ps_simple$mu_lodd_delta_ukr))
# Difference in odds
mean(exp ((ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr) - (ps_simple$mu_lodd_delta_ukr)))
HDI(exp ((ps_simple$mu_lodd_delta_ukr + ps_simple$beta_weekend_ukr) - (ps_simple$mu_lodd_delta_ukr)))


###log odds of beta_market effect on detection (delta) (anchored on non-market day) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_market_ukr)
HDI(ps_simple$beta_market_ukr)
#Converting to odds ratio
mean(exp(ps_simple$beta_market_italy))
HDI(exp(ps_simple$beta_market_italy))
#Marketday probability
mean(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))
# Market day odds of detection
mean(exp(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))
HDI(exp(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr))
#Non_marketday probability
mean(inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr))
# Non market-day odds of detection
mean(exp(ps_simple$mu_lodd_delta_ukr))
HDI(exp(ps_simple$mu_lodd_delta_ukr))
#Difference in means
mean( (inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr)) - inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI( (inv_logit(ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr)) - inv_logit(ps_simple$mu_lodd_delta_ukr))
# Difference in odds
mean(exp ((ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr) - (ps_simple$mu_lodd_delta_ukr)))
HDI(exp ((ps_simple$mu_lodd_delta_ukr + ps_simple$beta_market_ukr) - (ps_simple$mu_lodd_delta_ukr)))

###log odds of rain (beta_rain) effect on detection (delta) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_rain_italy)
HDI(ps_simple$beta_rain_italy)
#converting to odds ratios
mean(exp(ps_simple$beta_rain_ukr))
HDI(exp(ps_simple$beta_rain_ukr))
#Probability beta_rain effect
mean(inv_logit(ps_simple$beta_rain_italy))
HDI(inv_logit(ps_simple$beta_rain_italy))
#Difference in means
mean( (ps_simple$mu_lodd_delta_italy + ps_simple$beta_rain_italy) - ps_simple$mu_lodd_delta_italy)
mean( (inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_rain_italy)) - inv_logit(ps_simple$mu_lodd_delta_italy))
HDI( (inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_rain_italy)) - inv_logit(ps_simple$mu_lodd_delta_italy))


###log odds of temperature (beta_temp) effect on detection (delta) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_temp_italy)
HDI(ps_simple$beta_temp_italy)
#Converting to odds ratios
mean(exp(ps_simple$beta_temp_ukr))
HDI(exp(ps_simple$beta_temp_ukr))
#Probability beta_temp effect
mean(inv_logit(ps_simple$beta_temp_italy))
HDI(inv_logit(ps_simple$beta_temp_italy))
#Difference in means
mean( (ps_simple$mu_lodd_delta_italy + ps_simple$beta_temp_italy) - ps_simple$mu_lodd_delta_italy)
mean( (inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_temp_italy)) - inv_logit(ps_simple$mu_lodd_delta_italy))
HDI( (inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_temp_italy)) - inv_logit(ps_simple$mu_lodd_delta_italy))


###log odds of sex effect on survival (.1 phi) (female -> male) for Italy (_italy) and Ukraine (_ukr) ####
mean(ps_simple$beta_sex_ukr.1)
HDI(ps_simple$beta_sex_ukr.1)
# effect on probability scale
mean(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$sigma_dog_italy.1.2)) # sigma_dog[study site, 1recruitment/2survival/3detection]
sd(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$sigma_dog_italy.1.2))

#Converting to odds ratios
mean(exp(ps_simple$beta_sex_ukr.1))
HDI(exp(ps_simple$beta_sex_ukr.1))
# 
# #probabilty of sex effect on survival
# mean(inv_logit(ps_simple$beta_sex_italy.1))
# HDI(inv_logit(ps_simple$beta_sex_italy.1))

# total log odds of females' survival 
mean(ps_simple$mu_lodd_phi_ukr)
HDI(ps_simple$mu_lodd_phi_ukr)
#odds
mean(exp(ps_simple$mu_lodd_phi_ukr))
HDI(exp(ps_simple$mu_lodd_phi_ukr))

# total probability females' survival 
mean(inv_logit(ps_simple$mu_lodd_phi_ukr))
HDI(inv_logit(ps_simple$mu_lodd_phi_ukr))

# total log odds of males' survival
mean(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1)
HDI(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1)
#Odds
mean(exp(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))
HDI(exp(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))

# total probability of males' survival
mean(inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))
HDI(inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1))

# difference between males and females
mean(ps_simple$mu_lodd_phi_italy) - mean(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1)
HDI(ps_simple$mu_lodd_phi_italy) - mean(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1)
mean(ps_simple$mu_lodd_phi_ukr) - mean(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)
HDI(ps_simple$mu_lodd_phi_ukr) - mean(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)

mean(
  inv_logit(ps_simple$mu_lodd_phi_ukr) - inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)
) 
HDI(
  inv_logit(ps_simple$mu_lodd_phi_ukr) - inv_logit(ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)
)

mean(
  exp(ps_simple$mu_lodd_phi_ukr) - (ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)
) 
HDI(
  exp(ps_simple$mu_lodd_phi_ukr) - (ps_simple$mu_lodd_phi_ukr + ps_simple$beta_sex_ukr.1)
)



###log odds of sex effect on detection (.2 delta) (female -> male)####
mean(ps_simple$beta_sex_italy.2)
HDI(ps_simple$beta_sex_italy.2)
#Converting to odds ratios
mean(exp(ps_simple$beta_sex_ukr.2))
HDI(exp(ps_simple$beta_sex_ukr.2))

#probabilty of sex effect on detection (delta)
mean(inv_logit(ps_simple$beta_sex_italy.2))
HDI(inv_logit(ps_simple$beta_sex_italy.2))

# total log odds of females' detection (delta)
mean(ps_simple$mu_lodd_delta_italy)
HDI(ps_simple$mu_lodd_delta_italy)

#Odds of female detection (delta)
mean(exp(ps_simple$mu_lodd_delta_italy))
HDI(exp(ps_simple$mu_lodd_delta_italy))

# probability of female delta
mean(inv_logit(ps_simple$mu_lodd_delta_italy))
HDI(inv_logit(ps_simple$mu_lodd_delta_italy))

# total log odds of males' delta
mean(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2)
HDI(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2)
# Odds
mean(exp(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2))
HDI(exp(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2))
# probability of male delta
mean(inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2))
HDI(inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2))

# difference between males and females delta
mean(ps_simple$mu_lodd_delta_italy) - mean(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2)
mean(
  inv_logit(ps_simple$mu_lodd_delta_italy) - inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2)
) 
HDI(
  inv_logit(ps_simple$mu_lodd_delta_italy) - inv_logit(ps_simple$mu_lodd_delta_italy + ps_simple$beta_sex_italy.2)
) 

###Estimates for primary periods####
#Log odds of psi for italy/ukraine primary periods x(of five). 1= psi, or recruitment
ITPP1psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_time_italy.1.1
ITPP2psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_time_italy.2.1
ITPP3psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_time_italy.3.1
ITPP4psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_time_italy.4.1
ITPP5psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_time_italy.5.1
mean(ITPP5)
HDI(ITPP5)
#probability
mean(inv_logit(ITPP5))
HDI(inv_logit(ITPP5))


#Log odds of phi for Italy/Ukraine for primary periods x.2 = phi, or survival
# UAPP1phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.1.2 # This value doesn't really mean anything
UAPP2phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.2.2
UAPP3phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.3.2
UAPP4phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.4.2
UAPP5phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_time_ukr.5.2
mean(UAPP2phi)
HDI(UAPP2phi)
#probability
mean(inv_logit(UAPP2phi))
HDI(inv_logit(UAPP2phi))

# To work out the monthly survival probabilities - 1/3 for a 3-month interval, 1/6 for a 6-month interval
(0.9988881)^{1/3}
(0.9998149)^{1/6}

#Log odds of delta for italy/Ukraine primary periods x.3 = delta, or detection
UAPP1delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_time_ukr.1.3
UAPP2delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_time_ukr.2.3
UAPP3delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_time_ukr.3.3
UAPP4delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_time_ukr.4.3
UAPP5delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_time_ukr.5.3
mean(ITPP5delta)
HDI(ITPP5delta)
#probability
mean(inv_logit(ITPP5delta))
HDI(inv_logit(ITPP5delta))

#Comparing different primary periods
mean(ITPP4delta - ITPP5delta)
HDI(ITPP4delta - ITPP5delta)
# Odds ratios
mean(exp(UAPP4phi - UAPP5phi))
HDI(exp(UAPP4phi - UAPP5phi))

###Estimates for areas points####
# #Log odds of psi for italy areas x.1
# ITSA1psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.1.1
# ITSA2psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.2.1
# ITSA3psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.3.1
# ITSA4psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.4.1
# mean(ITSA4)
# HDI(ITSA4)
# #probability
# mean(inv_logit(ITSA4))
# HDI(inv_logit(ITSA4))
# 
# UASA1psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.1.1
# UASA2psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.2.1
# UASA3psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.3.1
# UASA4psi = ps_simple$mu_lodd_psi_italy + ps_simple$r_area_italy.4.1
# mean(UASA4psi)
# HDI(UASA4psi)
# #probability
# mean(inv_logit(UASA4psi))
# HDI(inv_logit(UASA4psi))

#Log odds of phi for ukr areas x.2
UASA1phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.1.2
UASA2phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.2.2
UASA3phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.3.2
UASA4phi = ps_simple$mu_lodd_phi_ukr + ps_simple$r_area_ukr.4.2
mean(UASA4phi)
HDI(UASA4phi)
#probabilUAy
mean(inv_logUA(UASA4phi))
HDI(inv_logUA(UASA4phi))

#Log odds of delta for ukr areas x.3
UASA1delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_area_ukr.1.3
UASA2delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_area_ukr.2.3
UASA3delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_area_ukr.3.3
UASA4delta = ps_simple$mu_lodd_delta_ukr + ps_simple$r_area_ukr.4.3
mean(UASA4delta)
HDI(UASA4delta)
#probabilUAy
mean(inv_logUA(UASA4delta))
HDI(inv_logUA(UASA4delta))

#Comparing different study areas
mean(ITSA1delta - ITSA3delta)
HDI(ITSA1delta - ITSA3delta)
#Probability
mean(inv_logit(ITSA3delta - ITSA4delta))
HDI(inv_logit(ITSA3delta - ITSA4delta))
#Odds ratio
mean(exp(UASA3delta - UASA4delta))
HDI(exp(UASA3delta - UASA4delta))

###Entry probability italy.primary period.study area####
mean(ps_simple$entry_prob_italy.5.4)
HDI(ps_simple$entry_prob_italy.5.4)
mean(inv_logit(ps_simple$entry_prob_italy.5.4))
HDI(inv_logit(ps_simple$entry_prob_italy.5.4))


### correlations between time points ####

### primary periods Italy

# phi parameter
decay_phi <- ps_simple[ ,"decay_italy.2"]
etasq_time_phi_italy <- ps_simple[ , "etasq_time_italy.2"]
sigma_time_phi_italy <- ps_simple$sigma_k_time_italy.2

italy_time_phi_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  italy_time_phi_corrs[[i]] <- k_cov_corr(
                                  periodic_decay_kernal(
                                    dmat = time_dmat_italy, decaytime = decay_phi[i], 
                                    etasq = etasq_time_phi_italy[i], 
                                    period = 12.0, 
                                    sigma = sigma_time_phi_italy[i]
                                  )
                                )
}

italy_time_phi_median <- matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy))#Reduce("+", italy_time_phi_corrs)/length(italy_time_phi_corrs)
italy_time_phi_hdi_low <- matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy))
italy_time_phi_hdi_high <- matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy))

for(k in 1:(ncol(time_dmat_italy)-1)){
  for(j in (k+1):ncol(time_dmat_italy)){
    italy_time_phi_median[k,j] <- median(unlist(lapply(italy_time_phi_corrs, function(x) x[k,j])))
    italy_time_phi_median[j,k] <- italy_time_phi_median[k,j]
    hdis <- HDI(unlist(lapply(italy_time_phi_corrs, function(x) x[k,j])))
    italy_time_phi_hdi_low[k,j] <- hdis[1]
    italy_time_phi_hdi_low[j,k] <- italy_time_phi_hdi_low[k,j]
    italy_time_phi_hdi_high[k,j] <- hdis[2]
    italy_time_phi_hdi_high[j,k] <- italy_time_phi_hdi_high[k,j]
  }
}

for(k in 1:ncol(time_dmat_italy)){
  italy_time_phi_hdi_low[k,k] <- italy_time_phi_hdi_high[k,k] <- italy_time_phi_median[k,k] <- 1.0
}

# delta parameter
decay_delta <- ps_simple[ ,"decay_italy.3"]
etasq_time_delta_italy <- ps_simple[ , "etasq_time_italy.3"]
sigma_time_delta_italy <- ps_simple$sigma_k_time_italy.3

italy_time_delta_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  italy_time_delta_corrs[[i]] <- k_cov_corr(
    periodic_decay_kernal(
      dmat = time_dmat_italy, decaytime = decay_delta[i], 
      etasq = etasq_time_delta_italy[i], 
      period = 12.0, 
      sigma = sigma_time_delta_italy[i]
    )
  )
}

italy_time_delta_median <-matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy)) #Reduce("+", italy_time_delta_corrs)/length(italy_time_delta_corrs)
italy_time_delta_hdi_low <- matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy))
italy_time_delta_hdi_high <- matrix(NA, ncol = ncol(time_dmat_italy), nrow = ncol(time_dmat_italy))

for(k in 1:(ncol(time_dmat_italy)-1)){
  for(j in (k+1):ncol(time_dmat_italy)){
    italy_time_delta_median[k,j] <- median(unlist(lapply(italy_time_delta_corrs, function(x) x[k,j])))
    italy_time_delta_median[j,k] <- italy_time_delta_median[k,j]
    hdis <- HDI(unlist(lapply(italy_time_delta_corrs, function(x) x[k,j])))
    italy_time_delta_hdi_low[k,j] <- hdis[1]
    italy_time_delta_hdi_low[j,k] <- italy_time_delta_hdi_low[k,j]
    italy_time_delta_hdi_high[k,j] <- hdis[2]
    italy_time_delta_hdi_high[j,k] <- italy_time_delta_hdi_high[k,j]
  }
}

for(k in 1:ncol(time_dmat_italy)){
  italy_time_delta_hdi_low[k,k] <- italy_time_delta_hdi_high[k,k] <- italy_time_delta_median[k,k] <- 1.0
}

### primary periods Ukraine

# phi parameter
decay_phi <- ps_simple[ ,"decay_ukr.2"]
etasq_time_phi_ukr <- ps_simple[ , "etasq_time_ukr.2"]
sigma_time_phi_ukr <- ps_simple$sigma_k_time_ukr.2

ukr_time_phi_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  ukr_time_phi_corrs[[i]] <- k_cov_corr(
    periodic_decay_kernal(
      dmat = time_dmat_ukr, decaytime = decay_phi[i], 
      etasq = etasq_time_phi_ukr[i], 
      period = 12.0, 
      sigma = sigma_time_phi_ukr[i]
    )
  )
}

ukr_time_phi_median <- matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr))#Reduce("+", ukr_time_phi_corrs)/length(ukr_time_phi_corrs)
ukr_time_phi_hdi_low <- matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr))
ukr_time_phi_hdi_high <- matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr))

for(k in 1:(ncol(time_dmat_ukr)-1)){
  for(j in (k+1):ncol(time_dmat_ukr)){
    ukr_time_phi_median[k,j] <- median(unlist(lapply(ukr_time_phi_corrs, function(x) x[k,j])))
    ukr_time_phi_median[j,k] <- ukr_time_phi_median[k,j]
    hdis <- HDI(unlist(lapply(ukr_time_phi_corrs, function(x) x[k,j])))
    ukr_time_phi_hdi_low[k,j] <- hdis[1]
    ukr_time_phi_hdi_low[j,k] <- ukr_time_phi_hdi_low[k,j]
    ukr_time_phi_hdi_high[k,j] <- hdis[2]
    ukr_time_phi_hdi_high[j,k] <- ukr_time_phi_hdi_high[k,j]
  }
}

for(k in 1:ncol(time_dmat_ukr)){
  ukr_time_phi_hdi_low[k,k] <- ukr_time_phi_hdi_high[k,k] <- ukr_time_phi_median[k,k] <- 1.0
}

# delta parameter
decay_delta <- ps_simple[ ,"decay_ukr.3"]
etasq_time_delta_ukr <- ps_simple[ , "etasq_time_ukr.3"]
sigma_time_delta_ukr <- ps_simple$sigma_k_time_ukr.3

ukr_time_delta_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  ukr_time_delta_corrs[[i]] <- k_cov_corr(
    periodic_decay_kernal(
      dmat = time_dmat_ukr, decaytime = decay_delta[i], 
      etasq = etasq_time_delta_ukr[i], 
      period = 12.0, 
      sigma = sigma_time_delta_ukr[i]
    )
  )
}

ukr_time_delta_median <-matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr)) #Reduce("+", ukr_time_delta_corrs)/length(ukr_time_delta_corrs)
ukr_time_delta_hdi_low <- matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr))
ukr_time_delta_hdi_high <- matrix(NA, ncol = ncol(time_dmat_ukr), nrow = ncol(time_dmat_ukr))

for(k in 1:(ncol(time_dmat_ukr)-1)){
  for(j in (k+1):ncol(time_dmat_ukr)){
    ukr_time_delta_median[k,j] <- median(unlist(lapply(ukr_time_delta_corrs, function(x) x[k,j])))
    ukr_time_delta_median[j,k] <- ukr_time_delta_median[k,j]
    hdis <- HDI(unlist(lapply(ukr_time_delta_corrs, function(x) x[k,j])))
    ukr_time_delta_hdi_low[k,j] <- hdis[1]
    ukr_time_delta_hdi_low[j,k] <- ukr_time_delta_hdi_low[k,j]
    ukr_time_delta_hdi_high[k,j] <- hdis[2]
    ukr_time_delta_hdi_high[j,k] <- ukr_time_delta_hdi_high[k,j]
  }
}

for(k in 1:ncol(time_dmat_ukr)){
  ukr_time_delta_hdi_low[k,k] <- ukr_time_delta_hdi_high[k,k] <- ukr_time_delta_median[k,k] <- 1.0
}

### Spatial Correlations ####

### areas italy

# phi parameter
etasq_area_phi_italy <- ps_simple[ , "etasq_area_italy.2"]
rhosq_area_phi_italy <- ps_simple$rhosq_area_italy.2
sigma_area_phi_italy <- ps_simple$sigma_k_area_italy.2

italy_area_phi_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  italy_area_phi_corrs[[i]] <- k_cov_corr(
    sqexp_cov_kernel(
      dmat = area_dmat_italy,
      etasq = etasq_area_phi_italy[i],
      rhosq = rhosq_area_phi_italy[i],
      sigma = sigma_area_phi_italy[i]
    )
  )
}

italy_area_phi_median <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))#Reduce("+", italy_area_phi_corrs)/length(italy_area_phi_corrs)
italy_area_phi_hdi_low <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))
italy_area_phi_hdi_high <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))


for(k in 1:(ncol(area_dmat_italy)-1)){
  for(j in (k+1):ncol(area_dmat_italy)){
    italy_area_phi_median[k,j] <- median(unlist(lapply(italy_area_phi_corrs, function(x) x[k,j])))
    italy_area_phi_median[j,k] <- italy_area_phi_median[k,j]
    hdis <- HDI(unlist(lapply(italy_area_phi_corrs, function(x) x[k,j])))
    italy_area_phi_hdi_low[k,j] <- hdis[1]
    italy_area_phi_hdi_low[j,k] <- italy_area_phi_hdi_low[k,j]
    italy_area_phi_hdi_high[k,j] <- hdis[2]
    italy_area_phi_hdi_high[j,k] <- italy_area_phi_hdi_high[k,j]
  }
}

for(k in 1:ncol(area_dmat_italy)){
  italy_area_phi_hdi_low[k,k] <- italy_area_phi_hdi_high[k,k] <- italy_area_phi_median[k,k] <- 1.0
}

# delta parameter
etasq_area_delta_italy <- ps_simple[ , "etasq_area_italy.3"]
rhosq_area_delta_italy <- ps_simple$rhosq_area_italy.3
sigma_area_delta_italy <- ps_simple$sigma_k_area_italy.3

italy_area_delta_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  italy_area_delta_corrs[[i]] <- k_cov_corr(
    sqexp_cov_kernel(
      dmat = area_dmat_italy,  
      etasq = etasq_area_delta_italy[i], 
      rhosq = rhosq_area_delta_italy[i],
      sigma = sigma_area_delta_italy[i]
    )
  )
}

italy_area_delta_median <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))#Reduce("+", italy_area_delta_corrs)/length(italy_area_delta_corrs)
italy_area_delta_hdi_low <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))
italy_area_delta_hdi_high <- matrix(NA, ncol = ncol(area_dmat_italy), nrow = ncol(area_dmat_italy))


for(k in 1:(ncol(area_dmat_italy)-1)){
  for(j in (k+1):ncol(area_dmat_italy)){
    italy_area_delta_median[k,j] <- median(unlist(lapply(italy_area_delta_corrs, function(x) x[k,j])))
    italy_area_delta_median[j,k] <- italy_area_delta_median[k,j]
    hdis <- HDI(unlist(lapply(italy_area_delta_corrs, function(x) x[k,j])))
    italy_area_delta_hdi_low[k,j] <- hdis[1]
    italy_area_delta_hdi_low[j,k] <- italy_area_delta_hdi_low[k,j]
    italy_area_delta_hdi_high[k,j] <- hdis[2]
    italy_area_delta_hdi_high[j,k] <- italy_area_delta_hdi_high[k,j]
  }
}

for(k in 1:ncol(area_dmat_italy)){
  italy_area_delta_hdi_low[k,k] <- italy_area_delta_hdi_high[k,k] <- italy_area_delta_median[k,k] <- 1.0
}

### areas Ukraine

# phi parameter
etasq_area_phi_ukr <- ps_simple[ , "etasq_area_ukr.2"]
rhosq_area_phi_ukr <- ps_simple$rhosq_area_ukr.2
sigma_area_phi_ukr <- ps_simple$sigma_k_area_ukr.2

ukr_area_phi_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  ukr_area_phi_corrs[[i]] <- k_cov_corr(
    sqexp_cov_kernel(
      dmat = area_dmat_ukr,
      etasq = etasq_area_phi_ukr[i],
      rhosq = rhosq_area_phi_ukr[i],
      sigma = sigma_area_phi_ukr[i]
    )
  )
}

ukr_area_phi_median <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))#Reduce("+", ukr_area_phi_corrs)/length(ukr_area_phi_corrs)
ukr_area_phi_hdi_low <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))
ukr_area_phi_hdi_high <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))


for(k in 1:(ncol(area_dmat_ukr)-1)){
  for(j in (k+1):ncol(area_dmat_ukr)){
    ukr_area_phi_median[k,j] <- median(unlist(lapply(ukr_area_phi_corrs, function(x) x[k,j])))
    ukr_area_phi_median[j,k] <- ukr_area_phi_median[k,j]
    hdis <- HDI(unlist(lapply(ukr_area_phi_corrs, function(x) x[k,j])))
    ukr_area_phi_hdi_low[k,j] <- hdis[1]
    ukr_area_phi_hdi_low[j,k] <- ukr_area_phi_hdi_low[k,j]
    ukr_area_phi_hdi_high[k,j] <- hdis[2]
    ukr_area_phi_hdi_high[j,k] <- ukr_area_phi_hdi_high[k,j]
  }
}

for(k in 1:ncol(area_dmat_ukr)){
  ukr_area_phi_hdi_low[k,k] <- ukr_area_phi_hdi_high[k,k] <- ukr_area_phi_median[k,k] <- 1.0
}

# delta parameter
etasq_area_delta_ukr <- ps_simple[ , "etasq_area_ukr.3"]
rhosq_area_delta_ukr <- ps_simple$rhosq_area_ukr.3
sigma_area_delta_ukr <- ps_simple$sigma_k_area_ukr.3

ukr_area_delta_corrs <- rep(list(list()), nrow(ps_simple))

for(i in 1:nrow(ps_simple)){
  ukr_area_delta_corrs[[i]] <- k_cov_corr(
    sqexp_cov_kernel(
      dmat = area_dmat_ukr,  
      etasq = etasq_area_delta_ukr[i], 
      rhosq = rhosq_area_delta_ukr[i],
      sigma = sigma_area_delta_ukr[i]
    )
  )
}

ukr_area_delta_median <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))#Reduce("+", ukr_area_delta_corrs)/length(ukr_area_delta_corrs)
ukr_area_delta_hdi_low <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))
ukr_area_delta_hdi_high <- matrix(NA, ncol = ncol(area_dmat_ukr), nrow = ncol(area_dmat_ukr))


for(k in 1:(ncol(area_dmat_ukr)-1)){
  for(j in (k+1):ncol(area_dmat_ukr)){
    ukr_area_delta_median[k,j] <- median(unlist(lapply(ukr_area_delta_corrs, function(x) x[k,j])))
    ukr_area_delta_median[j,k] <- ukr_area_delta_median[k,j]
    hdis <- HDI(unlist(lapply(ukr_area_delta_corrs, function(x) x[k,j])))
    ukr_area_delta_hdi_low[k,j] <- hdis[1]
    ukr_area_delta_hdi_low[j,k] <- ukr_area_delta_hdi_low[k,j]
    ukr_area_delta_hdi_high[k,j] <- hdis[2]
    ukr_area_delta_hdi_high[j,k] <- ukr_area_delta_hdi_high[k,j]
  }
}

for(k in 1:ncol(area_dmat_ukr)){
  ukr_area_delta_hdi_low[k,k] <- ukr_area_delta_hdi_high[k,k] <- ukr_area_delta_median[k,k] <- 1.0
}



### Plots of posterior N PP1 ####

Plot_Italy_1.1. <- ggplot(data=ps_simple, aes(N_hat_italy.1.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 1")                         #Primary period 1. Study site 1.

Plot_Italy_1.2. <- ggplot(data=ps_simple, aes(N_hat_italy.1.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 2")                         #Primary period 1. Study site 2.

Plot_Italy_1.3. <- ggplot(data=ps_simple, aes(N_hat_italy.1.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 3")                         #Primary period 1. Study site 3.

Plot_Italy_1.4. <- ggplot(data=ps_simple, aes(N_hat_italy.1.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 4")                         #Primary period 1. Study site 4.

Alt_Plot_Italy_1.1. <- Plot_Italy_1.1. + ylim(0, 2000)
Alt_Plot_Italy_1.2. <- Plot_Italy_1.2. + ylim(0, 2000)
Alt_Plot_Italy_1.3. <- Plot_Italy_1.3. + ylim(0, 2000)
Alt_Plot_Italy_1.4. <- Plot_Italy_1.4. + ylim(0, 2000)

Plot_ukr_1.1. <- ggplot(data=ps_simple, aes(N_hat_ukr.1.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 1")                         #Primary period 1. Study site 1.

Plot_ukr_1.2. <- ggplot(data=ps_simple, aes(N_hat_ukr.1.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 2")                         #Primary period 1. Study site 2.

Plot_ukr_1.3. <- ggplot(data=ps_simple, aes(N_hat_ukr.1.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 3")                         #Primary period 1. Study site 3.

Plot_ukr_1.4. <- ggplot(data=ps_simple, aes(N_hat_ukr.1.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 4")                         #Primary period 1. Study site 4.

Alt_Plot_ukr_1.1. <- Plot_ukr_1.1. + ylim(0, 2000)
Alt_Plot_ukr_1.2. <- Plot_ukr_1.2. + ylim(0, 2000)
Alt_Plot_ukr_1.3. <- Plot_ukr_1.3. + ylim(0, 2000)
Alt_Plot_ukr_1.4. <- Plot_ukr_1.4. + ylim(0, 2000)

Nhat_posterior_PP1 <- ggarrange(Alt_Plot_Italy_1.1., Alt_Plot_Italy_1.2., Alt_Plot_Italy_1.3., Alt_Plot_Italy_1.4.,
                       Alt_Plot_ukr_1.1., Alt_Plot_ukr_1.2., Alt_Plot_ukr_1.3., Alt_Plot_ukr_1.4.,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 4)

png("Nhat_posterior_PP1.png", res=800, width = 25, height =30, units="cm")
Nhat_posterior_PP1
dev.off()

### PP2 ####

Plot_Italy_2.1. <- ggplot(data=ps_simple, aes(N_hat_italy.2.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 1")                         #Primary period 2. Study site 1.

Plot_Italy_2.2. <- ggplot(data=ps_simple, aes(N_hat_italy.2.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 2")                         #Primary period 2. Study site 2.

Plot_Italy_2.3. <- ggplot(data=ps_simple, aes(N_hat_italy.2.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 3")                         #Primary period 2. Study site 3.

Plot_Italy_2.4. <- ggplot(data=ps_simple, aes(N_hat_italy.2.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 4")                         #Primary period 2. Study site 4.

Alt_Plot_Italy_2.1. <- Plot_Italy_2.1. + ylim(0, 2000)
Alt_Plot_Italy_2.2. <- Plot_Italy_2.2. + ylim(0, 2000)
Alt_Plot_Italy_2.3. <- Plot_Italy_2.3. + ylim(0, 2000)
Alt_Plot_Italy_2.4. <- Plot_Italy_2.4. + ylim(0, 2000)

Plot_ukr_2.1. <- ggplot(data=ps_simple, aes(N_hat_ukr.2.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 1")                         #Primary period 2. Study site 1.

Plot_ukr_2.2. <- ggplot(data=ps_simple, aes(N_hat_ukr.2.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 2")                         #Primary period 2. Study site 2.

Plot_ukr_2.3. <- ggplot(data=ps_simple, aes(N_hat_ukr.2.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 3")                         #Primary period 2. Study site 3.

Plot_ukr_2.4. <- ggplot(data=ps_simple, aes(N_hat_ukr.2.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 4")                         #Primary period 2. Study site 4.

Alt_Plot_ukr_2.1. <- Plot_ukr_2.1. + ylim(0, 2000)
Alt_Plot_ukr_2.2. <- Plot_ukr_2.2. + ylim(0, 2000)
Alt_Plot_ukr_2.3. <- Plot_ukr_2.3. + ylim(0, 2000)
Alt_Plot_ukr_2.4. <- Plot_ukr_2.4. + ylim(0, 2000)

Nhat_posterior_PP2 <- ggarrange(Alt_Plot_Italy_2.1., Alt_Plot_Italy_2.2., Alt_Plot_Italy_2.3., Alt_Plot_Italy_2.4.,
                                Alt_Plot_ukr_2.1., Alt_Plot_ukr_2.2., Alt_Plot_ukr_2.3., Alt_Plot_ukr_2.4.,
                                labels = c("A", "B", "C", "D"),
                                ncol = 2, nrow = 4)

png("Nhat_posterior_PP2.png", res=800, width = 25, height =30, units="cm")
Nhat_posterior_PP1
dev.off()

### PP3 ####

Plot_Italy_3.1. <- ggplot(data=ps_simple, aes(N_hat_italy.3.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 1")                         #Primary period 3. Study site 1.

Plot_Italy_3.2. <- ggplot(data=ps_simple, aes(N_hat_italy.3.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 2")                         #Primary period 3. Study site 2.

Plot_Italy_3.3. <- ggplot(data=ps_simple, aes(N_hat_italy.3.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 3")                         #Primary period 3. Study site 3.

Plot_Italy_3.4. <- ggplot(data=ps_simple, aes(N_hat_italy.3.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 4")                         #Primary period 3. Study site 4.

Alt_Plot_Italy_3.1. <- Plot_Italy_3.1. + ylim(0, 2000)
Alt_Plot_Italy_3.2. <- Plot_Italy_3.2. + ylim(0, 2000)
Alt_Plot_Italy_3.3. <- Plot_Italy_3.3. + ylim(0, 2000)
Alt_Plot_Italy_3.4. <- Plot_Italy_3.4. + ylim(0, 2000)

Plot_ukr_3.1. <- ggplot(data=ps_simple, aes(N_hat_ukr.3.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 1")                         #Primary period 3. Study site 1.

Plot_ukr_3.2. <- ggplot(data=ps_simple, aes(N_hat_ukr.3.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 2")                         #Primary period 3. Study site 2.

Plot_ukr_3.3. <- ggplot(data=ps_simple, aes(N_hat_ukr.3.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 3")                         #Primary period 3. Study site 3.

Plot_ukr_3.4. <- ggplot(data=ps_simple, aes(N_hat_ukr.3.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 4")                         #Primary period 3. Study site 4.

Alt_Plot_ukr_3.1. <- Plot_ukr_3.1. + ylim(0, 2000)
Alt_Plot_ukr_3.2. <- Plot_ukr_3.2. + ylim(0, 2000)
Alt_Plot_ukr_3.3. <- Plot_ukr_3.3. + ylim(0, 2000)
Alt_Plot_ukr_3.4. <- Plot_ukr_3.4. + ylim(0, 2000)

Nhat_posterior_PP3 <- ggarrange(Alt_Plot_Italy_3.1., Alt_Plot_Italy_3.2., Alt_Plot_Italy_3.3., Alt_Plot_Italy_3.4.,
                                Alt_Plot_ukr_3.1., Alt_Plot_ukr_3.2., Alt_Plot_ukr_3.3., Alt_Plot_ukr_3.4.,
                                labels = c("A", "B", "C", "D"),
                                ncol = 2, nrow = 4)

png("Nhat_posterior_PP3.png", res=800, width = 25, height =30, units="cm")
Nhat_posterior_PP1
dev.off()


### PP4 ####

Plot_Italy_4.1. <- ggplot(data=ps_simple, aes(N_hat_italy.4.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 1")                         #Primary period 4. Study site 1.

Plot_Italy_4.2. <- ggplot(data=ps_simple, aes(N_hat_italy.4.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 2")                         #Primary period 4. Study site 2.

Plot_Italy_4.3. <- ggplot(data=ps_simple, aes(N_hat_italy.4.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 3")                         #Primary period 4. Study site 3.

Plot_Italy_4.4. <- ggplot(data=ps_simple, aes(N_hat_italy.4.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 4")                         #Primary period 4. Study site 4.

Alt_Plot_Italy_4.1. <- Plot_Italy_4.1. + ylim(0, 2000)
Alt_Plot_Italy_4.2. <- Plot_Italy_4.2. + ylim(0, 2000)
Alt_Plot_Italy_4.3. <- Plot_Italy_4.3. + ylim(0, 2000)
Alt_Plot_Italy_4.4. <- Plot_Italy_4.4. + ylim(0, 2000)

Plot_ukr_4.1. <- ggplot(data=ps_simple, aes(N_hat_ukr.4.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 1")                         #Primary period 4. Study site 1.

Plot_ukr_4.2. <- ggplot(data=ps_simple, aes(N_hat_ukr.4.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 2")                         #Primary period 4. Study site 2.

Plot_ukr_4.3. <- ggplot(data=ps_simple, aes(N_hat_ukr.4.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 3")                         #Primary period 4. Study site 3.

Plot_ukr_4.4. <- ggplot(data=ps_simple, aes(N_hat_ukr.4.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 4")                         #Primary period 4. Study site 4.

Alt_Plot_ukr_4.1. <- Plot_ukr_4.1. + ylim(0, 2000)
Alt_Plot_ukr_4.2. <- Plot_ukr_4.2. + ylim(0, 2000)
Alt_Plot_ukr_4.3. <- Plot_ukr_4.3. + ylim(0, 2000)
Alt_Plot_ukr_4.4. <- Plot_ukr_4.4. + ylim(0, 2000)

Nhat_posterior_PP4 <- ggarrange(Alt_Plot_Italy_4.1., Alt_Plot_Italy_4.2., Alt_Plot_Italy_4.3., Alt_Plot_Italy_4.4.,
                                Alt_Plot_ukr_4.1., Alt_Plot_ukr_4.2., Alt_Plot_ukr_4.3., Alt_Plot_ukr_4.4.,
                                labels = c("A", "B", "C", "D"),
                                ncol = 2, nrow = 4)

png("Nhat_posterior_PP4.png", res=800, width = 25, height =30, units="cm")
Nhat_posterior_PP1
dev.off()


### PP5 ####

Plot_Italy_5.1. <- ggplot(data=ps_simple, aes(N_hat_italy.5.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 1")                         #Primary period 5. Study site 1.

Plot_Italy_5.2. <- ggplot(data=ps_simple, aes(N_hat_italy.5.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 2")                         #Primary period 5. Study site 2.

Plot_Italy_5.3. <- ggplot(data=ps_simple, aes(N_hat_italy.5.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 3")                         #Primary period 5. Study site 3.

Plot_Italy_5.4. <- ggplot(data=ps_simple, aes(N_hat_italy.5.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Pescara Site 4")                         #Primary period 5. Study site 4.

Alt_Plot_Italy_5.1. <- Plot_Italy_5.1. + ylim(0, 2000)
Alt_Plot_Italy_5.2. <- Plot_Italy_5.2. + ylim(0, 2000)
Alt_Plot_Italy_5.3. <- Plot_Italy_5.3. + ylim(0, 2000)
Alt_Plot_Italy_5.4. <- Plot_Italy_5.4. + ylim(0, 2000)

Plot_ukr_5.1. <- ggplot(data=ps_simple, aes(N_hat_ukr.5.1)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 1")                         #Primary period 5. Study site 1.

Plot_ukr_5.2. <- ggplot(data=ps_simple, aes(N_hat_ukr.5.2)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 2")                         #Primary period 5. Study site 2.

Plot_ukr_5.3. <- ggplot(data=ps_simple, aes(N_hat_ukr.5.3)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 3")                         #Primary period 5. Study site 3.

Plot_ukr_5.4. <- ggplot(data=ps_simple, aes(N_hat_ukr.5.4)) + 
  geom_histogram() +
  ylab("Frequency") +
  xlab("Population Size") +
  ggtitle("Lviv Site 4")                         #Primary period 5. Study site 4.

Alt_Plot_ukr_5.1. <- Plot_ukr_5.1. + ylim(0, 2000)
Alt_Plot_ukr_5.2. <- Plot_ukr_5.2. + ylim(0, 2000)
Alt_Plot_ukr_5.3. <- Plot_ukr_5.3. + ylim(0, 2000)
Alt_Plot_ukr_5.4. <- Plot_ukr_5.4. + ylim(0, 2000)

Nhat_posterior_PP5 <- ggarrange(Alt_Plot_Italy_5.1., Alt_Plot_Italy_5.2., Alt_Plot_Italy_5.3., Alt_Plot_Italy_5.4.,
                                Alt_Plot_ukr_5.1., Alt_Plot_ukr_5.2., Alt_Plot_ukr_5.3., Alt_Plot_ukr_5.4.,
                                labels = c("A", "B", "C", "D"),
                                ncol = 2, nrow = 4)

png("Nhat_posterior_PP5.png", res=800, width = 25, height =30, units="cm")
Nhat_posterior_PP1
dev.off()

# Plot detection and survival probabilities

fig1 = ggplot(ps_simple, aes(y=f_prob))+ 
  geom_dotplot(method = "histodot", binwidth = 0.75, 
               stackdir = "centerwhole")



# total probability females' survival 
f_prob <- mean(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1*-0.5))
f_HDI <- HDI(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1*-0.5))

# total log odds of males' survival
m_prob <- mean(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1*0.5))
m_HDI <- HDI(inv_logit(ps_simple$mu_lodd_phi_italy + ps_simple$beta_sex_italy.1*0.5))

sex_effect <- c(f_prob, m_prob)
