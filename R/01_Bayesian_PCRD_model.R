require(rstan)
rstan_options(auto_write = TRUE)
require(lubridate)
require(data.table)

logit <- function(x) log(x / (1-x))
  
inv_logit <- function(x) exp(x) / (1 + exp(x))

#------------------------------------------------------------------------------------------------------
# calculate population abundance statistics
# z_mat is a matrix of latent state estimates for each dog in an area across primary periods
abundance <- function(z_mat){
  n_i <- nrow(z_mat)
  n_t <- ncol(z_mat)
  ever_alive <- numeric(n_i)
  offsite <- onsite <- array(0, dim=c(n_i, n_t))
  N_hat_offsite <- N_hat_onsite <- numeric(n_t)
  superpop <- numeric()
  for(n in 1:n_i){
    if(any(z_mat[n,] == 3)  | any(z_mat[n,] == 4))
      ever_alive[n] = 1;
    for(t in 1:n_t){
      offsite[n,t] <- ifelse(z_mat[n,t] == 3, 1, 0)
      onsite[n,t]  <- ifelse(z_mat[n,t] == 4, 1, 0)
    }
  }
  for(t in 1:n_t){
    N_hat_offsite[t] <- sum(offsite[,t])
    N_hat_onsite[t] <- sum(onsite[,t])
  }
  superpop = sum(ever_alive)
  return(list(
    N_hat_offsite = N_hat_offsite, N_hat_onsite = N_hat_onsite, 
    superpop = superpop,
    offsite = offsite, onsite = onsite,
    ever_alive = ever_alive
  ))
}

#------------------------------------------------------------------------------------------------------
# Viterbi algorithm to get the most likely path for an individual

viterbi <- function(TM_array, EM_array, ch)
{
  ch_l <- length(ch)
  z_star <- numeric(ch_l)
  backptr <- bestprob <- array(NA, dim=c(ch_l, K))
  
  for(k in 1:K){
    bestprob[1,k] <- TM_array[1, k, 1] * EM_array[ch[1]+1, k, 1]
  }
  
  for(t in 2:ch_l){
    for(k in 1:K){
      bestprob[t,k] = 0;
      for(j in 1:K){
        p <- bestprob[t-1, j] * TM_array[j,k,t] * EM_array[ch[t]+1, k,t]
        if(p > bestprob[t,k]){
          backptr[t,k] = j;
          bestprob[t,k] = p;
        }
      }
    }
  }
  
  p_z_star <- max(bestprob[ch_l,])
  for(k in 1:K){
    if(bestprob[ch_l, k] == p_z_star){
      z_star[ch_l] = k
    }
  }
  
  for(t in 1:(ch_l-1)){
    z_star[ch_l - t] <- backptr[ch_l - t + 1, z_star[ch_l - t + 1]]
  }
  return(list(
    z_star = z_star,
    ch = ch,
    bestprob = bestprob
  ))
}
#-------------------------------------------------------------


first_cap <- function(y){
  n_p = length(y[1,])
  n_s = length(y[,1])
  cap = matrix(0, nrow = n_p, ncol=2)
  for(p in 1:n_p){
    cap[p,1] = p;
    for(s in 1:n_s){
      if(y[s,p] == 1){
        cap[p,2] = s;
        break;
      }
    }
  }
  if(sum(cap[,2]) > 0){
    for(p in 1:n_p){
      if(cap[p,2] > 0){
        return(cap[p,])
      }
    }
  }
  else{
    return(0)
  }
}

HDI <- function(x) coda::HPDinterval(coda::as.mcmc(x), prob=0.95) 

# ARRANGE ITALY DATA------------------------------------------------------------------------------------------------

d_italy = read.csv("PRDMR.csv")

###Italy dataset####


# Setting up arrays for Italy, by study area and sampling period.
# The dimensions are set based on the total number of individuals observed in a study area.

# Splitting up the columns for the different primary periods.
ch_cols_italy <- c(5:7,9:11,13:15,17:19,21:23)

# number of dogs in each area
n_area_italy <- c(
  nrow(d_italy[ d_italy$Study.Area %in% "ITSA1" , ch_cols_italy ]),
  nrow(d_italy[ d_italy$Study.Area %in% "ITSA2" , ch_cols_italy ]),
  nrow(d_italy[ d_italy$Study.Area %in% "ITSA3" , ch_cols_italy ]),
  nrow(d_italy[ d_italy$Study.Area %in% "ITSA4" , ch_cols_italy ])
)

# Splitting up for the different study areas.
ditaly_ch_s1 <- array(
  unlist(d_italy[ d_italy$Study.Area %in% "ITSA1" , ch_cols_italy ]),
  dim = c(n_area_italy[1], 3, 5)
)

ditaly_ch_s2 <- array(
  unlist(d_italy[ d_italy$Study.Area %in% "ITSA2" , ch_cols_italy ]),
  dim = c(n_area_italy[2], 3, 5)
)

ditaly_ch_s3 <- array(
  unlist(d_italy[ d_italy$Study.Area %in% "ITSA3" , ch_cols_italy ]),
  dim = c(n_area_italy[3], 3, 5)
)

ditaly_ch_s4 <- array(
  unlist(d_italy[ d_italy$Study.Area %in% "ITSA4" , ch_cols_italy ]),
  dim = c(n_area_italy[4], 3, 5)
)

# add in augmented individuals
# number of total individuals (seen and not seen)
n_augment_italy <- 50

area_list_italy <- rep(list(list()), 5) # list to hold the area data in the for loop
ditaly_ch_zi <- array(NA, dim=c(n_augment_italy, 3, 5, 4)) # new data array for the analysis

for(p in 1:5){ # loop through primary periods
  # fill in the area list
  area_list_italy[[1]][[p]] <- rbind(ditaly_ch_s1[,,p], array(0, dim=c(n_augment_italy - nrow(ditaly_ch_s1), 3)))
  area_list_italy[[2]][[p]] <- rbind(ditaly_ch_s2[,,p], array(0, dim=c(n_augment_italy - nrow(ditaly_ch_s2), 3)))
  area_list_italy[[3]][[p]] <- rbind(ditaly_ch_s3[,,p], array(0, dim=c(n_augment_italy - nrow(ditaly_ch_s3), 3)))
  area_list_italy[[4]][[p]] <- rbind(ditaly_ch_s4[,,p], array(0, dim=c(n_augment_italy - nrow(ditaly_ch_s4), 3)))
  for(a in 1:4){ # loop through areas
    # fill in zero inflated individuals into one array
    ditaly_ch_zi[,,,a] <- array(unlist(area_list_italy[[a]]), dim=c(n_augment_italy,3,5) )
  }
}


# now create matrix of predictors for each individual

### Weekend or weekday Italy ####
WeatherItaly = read.csv("Italy_weather.csv")
WeatherItaly[,"Date"] = as_date(dmy(WeatherItaly$Date))
days_italy = wday(WeatherItaly$Date)
WeatherItaly[,"Day"] = days_italy
SatSun_italy <- as.numeric(c(1, 7))
WeatherItaly[, "Weekend"] <- ifelse(WeatherItaly$Day %in% SatSun_italy ,
                                    "1",
                                    "0")

# put day-level predictors in the right format
mean_day_temp_italy <- array(WeatherItaly$Mean.temperature, dim=c(3,4,5))
temp_mean_italy <- mean(mean_day_temp_italy)
rain_italy <- array(ifelse(WeatherItaly$Rain == "No", 0, 1), dim=c(3,4,5))
market_italy <- array(ifelse(WeatherItaly$Marketday == "No", 0, 1), dim=c(3,4,5))
weekend_italy <- array(as.numeric(WeatherItaly$Weekend), dim=c(3,4,5))

## INDIVIDUAL-LEVEL PREDICTORS
#Sex

sex_a1_italy <- ifelse(d_italy[ d_italy$Study.Area %in% "ITSA1" , "Sex" ] == "Female", 0,
                 ifelse(d_italy[ d_italy$Study.Area %in% "ITSA1" , "Sex" ] == "Male", 1,
                        -100
                 )
)


sex_a2_italy <- ifelse(d_italy[ d_italy$Study.Area %in% "ITSA2" , "Sex" ] == "Female", 0,
                 ifelse(d_italy[ d_italy$Study.Area %in% "ITSA2" , "Sex" ] == "Male", 1,
                        -100
                 )
)

sex_a3_italy <- ifelse(d_italy[ d_italy$Study.Area %in% "ITSA3" , "Sex" ] == "Female", 0,
                 ifelse(d_italy[ d_italy$Study.Area %in% "ITSA3" , "Sex" ] == "Male", 1,
                        -100
                 )
)

sex_a4_italy <- ifelse(d_italy[ d_italy$Study.Area %in% "ITSA4" , "Sex" ] == "Female", 0,
                 ifelse(d_italy[ d_italy$Study.Area %in% "ITSA4" , "Sex" ] == "Male", 1,
                        -100
                 )
)

# add sex (unknown) for zi individuals
sex_zi_A1_italy <- c(sex_a1_italy, rep(-100, n_augment_italy - nrow(ditaly_ch_s1)))
sex_zi_A2_italy <- c(sex_a2_italy, rep(-100, n_augment_italy - nrow(ditaly_ch_s2)))
sex_zi_A3_italy <- c(sex_a3_italy, rep(-100, n_augment_italy - nrow(ditaly_ch_s3)))
sex_zi_A4_italy <- c(sex_a4_italy, rep(-100, n_augment_italy - nrow(ditaly_ch_s4)))

# create combined sex predictor
sex_italy <- matrix(c(sex_zi_A1_italy, sex_zi_A2_italy, sex_zi_A3_italy, sex_zi_A4_italy),
              nrow = n_augment_italy, ncol = 4)


## Area distances
area_dmat_italy <- read.csv("Italy_Distances.csv")[,-1]

# Time distances

time_dmat_italy <- matrix(c(
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

# ARRANGE UKRAINE DATA------------------------------------------------------------------------------------------------

d_ukr = read.csv("PRDMR.csv")

###ukr dataset####


# Setting up arrays for ukr, by study area and sampling period.
# The dimensions are set based on the total number of individuals observed in a study area.

# Splitting up the columns for the different primary periods.
ch_cols_ukr <- c(5:7,9:11,13:15,17:19,21:23)

# number of dogs in each area
n_area_ukr <- c(
  nrow(d_ukr[ d_ukr$Study.Area %in% "UASA1" , ch_cols_ukr ]),
  nrow(d_ukr[ d_ukr$Study.Area %in% "UASA2" , ch_cols_ukr ]),
  nrow(d_ukr[ d_ukr$Study.Area %in% "UASA3" , ch_cols_ukr ]),
  nrow(d_ukr[ d_ukr$Study.Area %in% "UASA4" , ch_cols_ukr ])
)

# Splitting up for the different study areas.
dukr_ch_s1 <- array(
  unlist(d_ukr[ d_ukr$Study.Area %in% "UASA1" , ch_cols_ukr ]),
  dim = c(n_area_ukr[1], 3, 5)
)

# change NA to -1
dukr_ch_s1[,2,3] <- -1

dukr_ch_s2 <- array(
  unlist(d_ukr[ d_ukr$Study.Area %in% "UASA2" , ch_cols_ukr ]),
  dim = c(n_area_ukr[2], 3, 5)
)

dukr_ch_s3 <- array(
  unlist(d_ukr[ d_ukr$Study.Area %in% "UASA3" , ch_cols_ukr ]),
  dim = c(n_area_ukr[3], 3, 5)
)

dukr_ch_s4 <- array(
  unlist(d_ukr[ d_ukr$Study.Area %in% "UASA4" , ch_cols_ukr ]),
  dim = c(n_area_ukr[4], 3, 5)
)

# add in augmented individuals
# number of total individuals (seen and not seen)
n_augment_ukr <- 300

area_list_ukr <- rep(list(list()), 5) # list to hold the area data in the for loop
dukr_ch_zi <- array(NA, dim=c(n_augment_ukr, 3, 5, 4)) # new data array for the analysis

for(p in 1:5){ # loop through primary periods
  # fill in the area list
  area_list_ukr[[1]][[p]] <- rbind(dukr_ch_s1[,,p], array(0, dim=c(n_augment_ukr - nrow(dukr_ch_s1), 3)))
  area_list_ukr[[2]][[p]] <- rbind(dukr_ch_s2[,,p], array(0, dim=c(n_augment_ukr - nrow(dukr_ch_s2), 3)))
  area_list_ukr[[3]][[p]] <- rbind(dukr_ch_s3[,,p], array(0, dim=c(n_augment_ukr - nrow(dukr_ch_s3), 3)))
  area_list_ukr[[4]][[p]] <- rbind(dukr_ch_s4[,,p], array(0, dim=c(n_augment_ukr - nrow(dukr_ch_s4), 3)))
  for(a in 1:4){ # loop through areas
    # fill in zero inflated individuals into one array
    dukr_ch_zi[,,,a] <- array(unlist(area_list_ukr[[a]]), dim=c(n_augment_ukr,3,5) )
  }
}

dukr_ch_zi[,2,3,1] <- -1


# now create matrix of predictors for each individual

### Weekend or weekday ukr ####
Weatherukr = read.csv("Ukraine_weather.csv")
Weatherukr[,"Date"] = as_date(dmy(Weatherukr$Date))
days_ukr = wday(Weatherukr$Date)
Weatherukr[,"Day"] = days_ukr
SatSun_ukr <- as.numeric(c(1, 7))
Weatherukr[, "Weekend"] <- ifelse(Weatherukr$Day %in% SatSun_ukr ,
                                    "1",
                                    "0")

# put day-level predictors in the right format
mean_day_temp_ukr <- array(Weatherukr$Mean.Temperature, dim=c(3,4,5))
temp_mean_ukr <- mean(mean_day_temp_ukr)
rain_ukr <- array(ifelse(Weatherukr$Rain == "No", 0, 1), dim=c(3,4,5))
market_ukr <- array(ifelse(Weatherukr$Marketday == "No", 0, 1), dim=c(3,4,5))
market_ukr[2,1,3] <- 0
weekend_ukr <- array(as.numeric(Weatherukr$Weekend), dim=c(3,4,5))

## INDIVIDUAL-LEVEL PREDICTORS
#Sex

sex_a1_ukr <- ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA1" , "Sex" ] == "Female", 0,
                       ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA1" , "Sex" ] == "Male", 1,
                              -100
                       )
)


sex_a2_ukr <- ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA2" , "Sex" ] == "Female", 0,
                       ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA2" , "Sex" ] == "Male", 1,
                              -100
                       )
)

sex_a3_ukr <- ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA3" , "Sex" ] == "Female", 0,
                       ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA3" , "Sex" ] == "Male", 1,
                              -100
                       )
)

sex_a4_ukr <- ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA4" , "Sex" ] == "Female", 0,
                       ifelse(d_ukr[ d_ukr$Study.Area %in% "UASA4" , "Sex" ] == "Male", 1,
                              -100
                       )
)

# add sex (unknown) for zi individuals
sex_zi_A1_ukr <- c(sex_a1_ukr, rep(-100, n_augment_ukr - nrow(dukr_ch_s1)))
sex_zi_A2_ukr <- c(sex_a2_ukr, rep(-100, n_augment_ukr - nrow(dukr_ch_s2)))
sex_zi_A3_ukr <- c(sex_a3_ukr, rep(-100, n_augment_ukr - nrow(dukr_ch_s3)))
sex_zi_A4_ukr <- c(sex_a4_ukr, rep(-100, n_augment_ukr - nrow(dukr_ch_s4)))

# create combined sex predictor
sex_ukr <- matrix(c(sex_zi_A1_ukr, sex_zi_A2_ukr, sex_zi_A3_ukr, sex_zi_A4_ukr),
                    nrow = n_augment_ukr, ncol = 4)

## Area distances

area_dmat_ukr <- read.csv("Ukraine_Distances.csv")[,-1]
#area_dmat_ukr_scaled_10 <- area_dmat_ukr/10

# Time distances

time_dmat_ukr <- matrix(c(
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


#-FIT STAN three state, both countries-------------------------------------------------------------------------------------------------------------
# prepare data for Stan

stan_data <- list(
  n_pop = 2,
  n_m = c(n_augment_italy, n_augment_ukr),
  n_s = 3,
  n_p = 5,
  n_k = 3,
  n_a = 4,
  # italy
  ch_italy = ditaly_ch_zi,
  sex_italy = sex_italy,
  dmat_area_italy = area_dmat_italy,
  dmat_time_italy = time_dmat_italy,
  temp_italy = mean_day_temp_italy,
  mean_temp_italy = temp_mean_italy,
  rain_italy = rain_italy,
  market_italy = market_italy,
  weekend_italy = weekend_italy,
  # ukraine
  ch_ukr = dukr_ch_zi,
  sex_ukr = sex_ukr,
  dmat_area_ukr = area_dmat_ukr,
  dmat_time_ukr = time_dmat_ukr,
  temp_ukr = mean_day_temp_ukr,
  mean_temp_ukr = temp_mean_ukr,
  rain_ukr = rain_ukr,
  market_ukr = market_ukr,
  weekend_ukr = weekend_ukr
)

fit_stan <- stan(
  file = "lauren-hmm-three-state-dapx-two-pop.stan",
  data = stan_data, 
  include = FALSE, pars = c("TM_italy","EM_italy","TM_ukr","EM_ukr"),
  warmup = 2500, iter = 5000, 
  chains = 4, cores = 4,
  control = list(adapt_delta = 0.95), 
  seed = 2019, 
  sample_file = "lauren-mcmc-samples.csv"
)

capture.output(print(fit_stan, 
                     pars = c(
                       "mu_lodd_psi_italy", "mu_lodd_phi_italy", "mu_lodd_delta_italy", "sigma_dog_italy",
                       "beta_sex_italy", "beta_weekend_italy", "beta_market_italy", "beta_temp_italy", "beta_rain_italy",
                       "etasq_area_italy", "etasq_time_italy", "rhosq_area_italy", "decay_italy", "sigma_k_area_italy", 
                       "sigma_k_time_italy", "kappa_sex_italy", "N_hat_italy", "superpop_italy", "entry_prob_italy",
                       "mu_lodd_psi_ukr", "mu_lodd_phi_ukr", "mu_lodd_delta_ukr", "sigma_dog_ukr",
                       "beta_sex_ukr", "beta_weekend_ukr", "beta_market_ukr", "beta_temp_ukr", "beta_rain_ukr",
                       "etasq_area_ukr", "etasq_time_ukr", "rhosq_area_ukr", "decay_ukr", "sigma_k_area_ukr", 
                       "sigma_k_time_ukr", "kappa_sex_ukr",  "N_hat_ukr", "superpop_ukr", "entry_prob_ukr",
                       "psi_italy", "phi_italy", "delta_italy",
                       "psi_ukr", "phi_ukr", "delta_ukr"
                     )),
               file = "two-pop-output-summary.csv")

fwrite(x = as.matrix(fit_stan), file = "two-pop-posterior-dist.csv", row.names = FALSE)

ps <- as.matrix(fread("lauren-mcmc-samples_1.csv"))

simple_pars <- c("mu_lodd_psi_italy", "mu_lodd_phi_italy", "mu_lodd_delta_italy", "sigma_dog_italy",
                 "beta_sex_italy", "beta_weekend_italy", "beta_market_italy", "beta_temp_italy", "beta_rain_italy",
                 "N_hat_italy", "superpop_italy", "entry_prob_italy",
                 "mu_lodd_psi_ukr", "mu_lodd_phi_ukr", "mu_lodd_delta_ukr", "sigma_dog_ukr",
                 "beta_sex_ukr", "beta_weekend_ukr", "beta_market_ukr", "beta_temp_ukr", "beta_rain_ukr",
                 "N_hat_ukr", "superpop_ukr", "entry_prob_ukr")

ps_simple <- data.frame(rows = 1:nrow(ps))

for(i in 1:length(simple_pars)){
  get_cols <- as.matrix(ps[,grep(simple_pars[i], colnames(ps))])
  if(ncol(get_cols)==1){
    colname_ <- simple_pars[i]
    ps_simple <- data.frame(ps_simple, ps[,grep(simple_pars[i], colnames(ps))])
    names(ps_simple)[ncol(ps_simple)] <- colname_
  }
  else{
    ps_simple <- data.frame(ps_simple, ps[,grep(simple_pars[i], colnames(ps))])
  }
  
}

write.csv(ps_simple, "lauren-mcmc-samples-simple.csv", row.names = FALSE)

## read in other files

ps_simple <- read.csv("lauren-mcmc-samples-simple.csv")

ps <- as.matrix(fread("lauren-mcmc-samples_2.csv"))

simple_pars <- c("mu_lodd_psi_italy", "mu_lodd_phi_italy", "mu_lodd_delta_italy", "sigma_dog_italy",
                 "beta_sex_italy", "beta_weekend_italy", "beta_market_italy", "beta_temp_italy", "beta_rain_italy",
                 "N_hat_italy", "superpop_italy", "entry_prob_italy",
                 "mu_lodd_psi_ukr", "mu_lodd_phi_ukr", "mu_lodd_delta_ukr", "sigma_dog_ukr",
                 "beta_sex_ukr", "beta_weekend_ukr", "beta_market_ukr", "beta_temp_ukr", "beta_rain_ukr",
                 "N_hat_ukr", "superpop_ukr", "entry_prob_ukr")

temp_ps <- data.frame(rows = 1:nrow(ps))

for(i in 1:length(simple_pars)){
  get_cols <- as.matrix(ps[,grep(simple_pars[i], colnames(ps))])
  if(ncol(get_cols)==1){
    colname_ <- simple_pars[i]
    temp_ps <- data.frame(temp_ps, ps[,grep(simple_pars[i], colnames(ps))])
    names(temp_ps)[ncol(temp_ps)] <- colname_
  }
  else{
    temp_ps <- data.frame(temp_ps, ps[,grep(simple_pars[i], colnames(ps))])
  }
}

rbind(ps_simple, temp_ps)



#-FIT STAN-------------------------------------------------------------------------------------------------------------
# prepare data for Stan

stan_italy_data <- list(
  n_m = n_augment,
  n_s = 3,
  n_p = 5,
  n_k = 4,
  n_a = 4,
  y = ditaly_ch_zi,
  sex = sex,
  distinct = distinct,
  distinct_mat = distinct_model_matrix,
  dmat_area = area_dmat_scaled_10,
  dmat_time = time_dmat,
  temp = mean_temp,
  rain = rain,
  market = market,
  weekend = weekend,
  first_caps = first_caps
)

fit_italy <- stan(
  file = "Lauren-PCRD-ITALY-HMM-DAPX.stan",
  data = stan_italy_data,
  warmup = 250, iter = 500, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

capture.output(print(fit_italy), file = "Italy-output-summary.csv")

write.csv(as.matrix(fit_italy), file = "Italy-posterior-dist.csv")

### without misclassification
# prepare data for Stan

stan_italy_data <- list(
  n_m = n_augment,
  n_s = 3,
  n_p = 5,
  n_k = 4,
  n_a = 4,
  y = ditaly_ch_zi,
  sex = sex,
  dmat_area = area_dmat_scaled_10,
  dmat_time = time_dmat,
  temp = mean_temp,
  rain = rain,
  market = market,
  weekend = weekend
)

fit_italy <- stan(
  file = "Lauren-PCRD-ITALY-withoutmisclass-HMM-DAPX.stan",
  data = stan_italy_data,
  warmup = 2500, iter = 5000, chains = 4, cores = 4,
  control = list(adapt_delta=0.95)
)

capture.output(print(fit_italy, 
                     pars = c("mu_psi","mu_phi","mu_lambda","mu_gamma_prime","mu_gamma_double_prime",
                              "mu_delta", "beta_weekend", "beta_market", "beta_sex", "superpop", "N_hat")
                     ), 
               file = "Italy-output-summary.csv")

fwrite(as.matrix(fit_italy), file = "D:/Lauren_Model_Output/Italy-posterior-dist.csv")

italy = as.data.frame(fread("E:/Lauren_Model_Output/Italy-posterior-dist.csv"))

inv_logit(italy$mu_lodd_delta)

###Extracting parameters####

#Log odds of xxx
mean(ps_simple$mu_lodd_delta_ukr)
HDI(ps_simple$mu_lodd_delta_ukr)
#probability mu_psi
mean(inv_logit(ps_simple$mu_lodd_delta_ukr))
HDI(inv_logit(ps_simple$mu_lodd_delta_ukr))


#Superpopulation
mean(ps_simple$superpop_ukr.4)
HDI(ps_simple$superpop_ukr.4)

#N_hat
mean(ps_simple$N_hat_ukr.5.4)
HDI(ps_simple$N_hat_ukr.5.4)


### Explanatory variables ####

samples <- sample(1:nrow(italy), 500, replace=TRUE) #nrow(italy)
n_samples <- length(samples)
TMs <- italy[,grep("TM", colnames(italy))]
EMs <- italy[,grep("EM", colnames(italy))]
TM_all_array <- array(NA, dim=c(4,4,n_augment,n_samples,5,4))
EM_all_array <- array(NA, dim=c(2,4,n_augment,n_samples,3,5,4))
z_star <- array(NA, dim=c(n_augment, 5, n_samples, 4))
N_hats_onsite <- array(NA, dim=c(5, 4, n_samples))
N_hats_offsite <- array(NA, dim=c(5, 4, n_samples))
superpops <- array(NA, dim=c(4,n_samples))

K = 4

for(i in seq_along(samples)){# MCMC sample i
  for(a in 1:4){# area a
    for(n in 1:n_augment){# individual n
      for(t in 1:5){# primary period t
        
        # GET TRANSMISSION MATRIX
        for(k in 1:(K-1)){# rows
          for(j in (k+1):K){ # columns
            colname_1 <-  paste("TM[",k,",",j,",",n,",",t,",",a,"]",sep="")
            TM_all_array[k,j,n,i,t,a] <- TMs[samples[i],colname_1]
            colname_2 <-  paste("TM[",j,",",k,",",n,",",t,",",a,"]",sep="")
            TM_all_array[j,k,n,i,t,a] <-  TMs[samples[i],colname_2]
          }
        }
        for(k in 1:K){# row == column
          colname_ <-  paste("TM[",k,",",k,",",n,",",t,",",a,"]",sep="")
          TM_all_array[k,k,n,i,t,a] <- TMs[samples[i],colname_]
        }
        
        # GET EMISSIONS MATRIX
        for(s in 1:3){# for each secondary period
          for(k in 1:2){# rows
            for(j in 1:K){ # columns
              colname_1 <-  paste("EM[",k,",",j,",",n,",",s,",",t,",",a,"]",sep="")
              EM_all_array[k,j,n,i,s,t,a] <- EMs[samples[i],colname_1]
              #colname_2 <-  paste("EM[",j,",",k,",",n,",",s,",",t,",",a,"]",sep="")
              #EM_all_array[j,k,n,i,s,t,a] <-  EMs[i,colname_2]
            }#j
          }#k
        }#s
        # calculate average emissions matrix
        average_EM <- array(NA, dim=c(2,K,5)) 
        for(t in 1:5){
          average_EM[,,t] <- apply(EM_all_array[,,n,i,,t,a], c(1,2), mean)
        }
      }#t
      # get the capture history for primary periods
      ch_ <- numeric(5)
      for(t in 1:5){
        # if seen during primary period, ch = 1
        ch_[t] <- ifelse(sum(ditaly_ch_zi[n,,t,a]) > 0, 1, 0)
      }
      # get the latent states
      z_star[n,,i,a] <- viterbi(TM_array = TM_all_array[,,n,i,,a], 
                                EM_array = average_EM,
                                ch = ch_
                                )$z_star
    }#n
    
    # CALCULATE ABUNDANCE STATISTICS
    N_hats_onsite[,a,i] <- abundance(z_mat = z_star[,,i,a])$N_hat_onsite
    N_hats_offsite[,a,i] <- abundance(z_mat = z_star[,,i,a])$N_hat_offsite
    superpops[a,i] <- abundance(z_mat = z_star[,,i,a])$superpop
  }#a
}#i

superpops_mean <- rowMeans(superpops)
superpops_hdi <- array(NA, dim=c(4,2))
for(a in 1:4){
  superpops_hdi[a,] <- HDI(superpops[a,])
} 


  



