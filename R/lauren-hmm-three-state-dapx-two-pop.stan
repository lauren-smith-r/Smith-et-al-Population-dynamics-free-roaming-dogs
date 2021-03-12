functions{
  int[] row_sums(int[,] y){
    int n_p = size(y[1,]);
    int n_s = size(y[,1]);
    int rs[n_p];
    for(p in 1:n_p){
      int c = 0;
      for(s in 1:n_s){
        if(y[s,p] > 0)
          c += 1;
      }
      rs[p] = c;
    }
    return rs;
  }

  int last(int[,] y){
    int n_p = size(y[1,]);
    int lc = 0;
    int rs[n_p] = row_sums(y);
    for(p in 1:n_p){
      if(rs[p] > 0){
        lc = p;
      }
    }
    return lc;
  }
}


data{
  int<lower=1> n_pop;
  int<lower=0> n_m[n_pop];                      // number individuals
  int<lower=0> n_s;                             // number of secondary periods
  int<lower=0> n_p;                             // number of primary periods
  int<lower=1> n_k;                             // number of states
  int<lower=1> n_a;                             // number of areas
  // italy data
  int<upper=1> ch_italy[n_m[1],n_s,n_p,n_a];    // capture histories
  int sex_italy[n_m[1], n_a];
  matrix[n_a, n_a] dmat_area_italy;
  matrix[n_p, n_p] dmat_time_italy;
  real temp_italy[n_s, n_a, n_p];
  real mean_temp_italy;
  int rain_italy[n_s, n_a, n_p];
  int market_italy[n_s, n_a, n_p];
  int weekend_italy[n_s, n_a, n_p];
  // ukraine data
  int<upper=1> ch_ukr[n_m[2],n_s,n_p,n_a];
  int sex_ukr[n_m[2], n_a];
  matrix[n_a, n_a] dmat_area_ukr;
  matrix[n_p, n_p] dmat_time_ukr;
  real temp_ukr[n_s, n_a, n_p];
  real mean_temp_ukr;
  int rain_ukr[n_s, n_a, n_p];
  int market_ukr[n_s, n_a, n_p];
  int weekend_ukr[n_s, n_a, n_p];
}

transformed data{
  int<upper=2> ch_italy_add1[n_m[1], n_s, n_p, n_a];
  int<upper=2> ch_ukr_add1[n_m[2], n_s, n_p, n_a];
  matrix[n_m[1], n_a] sex_italy_c;
  matrix[n_m[2], n_a] sex_ukr_c;
  real temp_italy_c[n_s, n_a, n_p];
  real temp_ukr_c[n_s, n_a, n_p];
  real<lower=0> pi = 3.14;

  for(m in 1:n_m[1]){
    for(s in 1:n_s){
      for(p in 1:n_p){
        for(a in 1:n_a){

          if(ch_italy[m,s,p,a] > -1)
            ch_italy_add1[m,s,p,a] = ch_italy[m,s,p,a] + 1;

          if(sex_italy[m,a] != -100){
            sex_italy_c[m,a] = sex_italy[m,a] - 0.5;
          }
          temp_italy_c[s,a,p] = temp_italy[s,a,p] - mean_temp_italy;
        }
      }
    }
  }
   for(m in 1:n_m[2]){
    for(s in 1:n_s){
      for(p in 1:n_p){
        for(a in 1:n_a){
          if(ch_ukr[m,s,p,a] > -1)
            ch_ukr_add1[m,s,p,a] = ch_ukr[m,s,p,a] + 1;

          if(sex_ukr[m,a] != -100){
            sex_ukr_c[m,a] = sex_ukr[m,a] - 0.5;
          }
          temp_ukr_c[s,a,p] = temp_ukr[s,a,p] - mean_temp_ukr;
        }
      }
    }
  }
}

parameters{
  // italy
  real mu_lodd_psi_italy;                       // recruitment prob
  real mu_lodd_phi_italy;                       // survival probs
  real mu_lodd_delta_italy;                     // detection probs
  matrix[n_m[1],n_a] z_dog_italy[3];
  matrix<lower=0>[n_a,3] sigma_dog_italy;
  matrix[n_p, 3] r_time_italy;
  matrix[n_a, 3] r_area_italy;
  vector[2] beta_sex_italy;
  real beta_weekend_italy;
  real beta_market_italy;
  real beta_temp_italy;
  real beta_rain_italy;
  real<lower=0> etasq_area_italy[3];
  real<lower=0> etasq_time_italy[3];
  real<lower=0> rhosq_area_italy[3];
  real<lower=0> decay_italy[3];
  real<lower=0> sigma_k_area_italy[3];
  real<lower=0> sigma_k_time_italy[3];
  real<lower=0,upper=1> kappa_sex_italy;
  // ukraine
  real mu_lodd_psi_ukr;                       // recruitment prob
  real mu_lodd_phi_ukr;                       // survival probs
  real mu_lodd_delta_ukr;                     // detection probs
  matrix[n_m[2],n_a] z_dog_ukr[3];
  matrix<lower=0>[n_a,3] sigma_dog_ukr;
  matrix[n_p, 3] r_time_ukr;
  matrix[n_a, 3] r_area_ukr;
  vector[2] beta_sex_ukr;
  real beta_weekend_ukr;
  real beta_market_ukr;
  real beta_temp_ukr;
  real beta_rain_ukr;
  real<lower=0,upper=1> kappa_sex_ukr;
  real<lower=0> etasq_area_ukr[3];
  real<lower=0> etasq_time_ukr[3];
  real<lower=0> rhosq_area_ukr[3];
  real<lower=0> decay_ukr[3];
  real<lower=0> sigma_k_area_ukr[3];
  real<lower=0> sigma_k_time_ukr[3];
}

transformed parameters{
  // italy
  real<lower=0,upper=1> EM_italy[2,  n_k,n_m[1],n_s,n_p,n_a];                   // emission matrix
  real<lower=0,upper=1> TM_italy[n_k,n_k,n_m[1],n_p,n_a];                       // transmission matrix
  real<lower=0,upper=1> psi_italy[n_m[1],n_p,n_a];
  real<lower=0,upper=1> phi_italy[n_m[1],n_p,n_a];
  real<lower=0,upper=1> delta_italy[n_m[1],n_s,n_p,n_a];
  matrix[n_a,n_a] K_area_italy[3];
  matrix[n_p,n_p] K_time_italy[3];
  // ukraine
  real<lower=0,upper=1> EM_ukr[2,  n_k,n_m[2],n_s,n_p,n_a];                   // emission matrix
  real<lower=0,upper=1> TM_ukr[n_k,n_k,n_m[2],n_p,n_a];                       // transmission matrix
  real<lower=0,upper=1> psi_ukr[n_m[2],n_p,n_a];
  real<lower=0,upper=1> phi_ukr[n_m[2],n_p,n_a];
  real<lower=0,upper=1> delta_ukr[n_m[2],n_s,n_p,n_a];
  matrix[n_a,n_a] K_area_ukr[3];
  matrix[n_p,n_p] K_time_ukr[3];

  /* Squared exponentional kernel for area random effect */
    for(k in 1:3){
      for(i in 1:(n_a-1)){
        for(j in (i+1):n_a){
          // italy
          K_area_italy[k,i,j] = etasq_area_italy[k]*exp(-rhosq_area_italy[k]*pow(dmat_area_italy[i,j],2));
          K_area_italy[k,j,i] = K_area_italy[k,i,j];
          // ukraine
          K_area_ukr[k,i,j] = etasq_area_ukr[k]*exp(-rhosq_area_ukr[k]*pow(dmat_area_ukr[i,j],2));
          K_area_ukr[k,j,i] = K_area_ukr[k,i,j];
        }
      }
    }

  for(k in 1:3){
    for(a in 1:n_a){
      K_area_italy[k,a,a] = etasq_area_italy[k] + sigma_k_area_italy[k];
      K_area_ukr[k,a,a] = etasq_area_ukr[k] + sigma_k_area_ukr[k];
    }
  }

  /* Periodic-delay kernel for primary period random effect */
    for(k in 1:3){
      for(i in 1:(n_p-1)){
        for(j in (i+1):n_p){
          K_time_italy[k,i,j] = etasq_time_italy[k]*exp(-0.5 * dmat_time_italy[i,j]/decay_italy[k] - (2.0*pow(sin(pi*(dmat_time_italy[i,j])),2))/12.0);
          K_time_italy[k,j,i] = K_time_italy[k,i,j];
          K_time_ukr[k,i,j] = etasq_time_ukr[k]*exp(-0.5 * dmat_time_ukr[i,j]/decay_ukr[k] - (2.0*pow(sin(pi*(dmat_time_ukr[i,j])),2))/12.0);
          K_time_ukr[k,j,i] = K_time_ukr[k,i,j];
        }
      }
    }

  for(k in 1:3){
    for(p in 1:n_p){
      K_time_italy[k,p,p] = etasq_time_italy[k] + sigma_k_time_italy[k];
      K_time_ukr[k,p,p] = etasq_time_ukr[k] + sigma_k_time_ukr[k];
    }
  }

  /* ITALY */
  /* fill the transition & emission matrixes */
    for(m in 1:n_m[1]){
      for(p in 1:n_p){
        for(a in 1:n_a){

          psi_italy[m,p,a] = inv_logit(mu_lodd_psi_italy + r_time_italy[p,1] + r_area_italy[a,1] +  z_dog_italy[1,m,a] * sigma_dog_italy[a,1]);

          /* if sex is observed  */
          if(sex_italy[m,a] != -100)
            phi_italy[m,p,a] = inv_logit(mu_lodd_phi_italy + r_time_italy[p,2] + r_area_italy[a,2] + z_dog_italy[2,m,a] * sigma_dog_italy[a,2] + beta_sex_italy[1] * (sex_italy_c[m,a]));
          /* if sex is missing  */
          if(sex_italy[m,a] == -100)
            phi_italy[m,p,a] = (1 - kappa_sex_italy) * inv_logit(mu_lodd_phi_italy + r_time_italy[p,2] + r_area_italy[a,2] + z_dog_italy[2,m,a] * sigma_dog_italy[a,2] + beta_sex_italy[1]*-0.5) +
                                    kappa_sex_italy  * inv_logit(mu_lodd_phi_italy + r_time_italy[p,2] + r_area_italy[a,2] + z_dog_italy[2,m,a] * sigma_dog_italy[a,2] + beta_sex_italy[1]* 0.5);

          /* TRANSMISSION MATRIX */
          TM_italy[1,1,m,p,a] = 1 - psi_italy[m,p,a];                                     // nye->nye
          TM_italy[1,2,m,p,a] = 0.0;                                                // nye->dead
          TM_italy[1,3,m,p,a] = psi_italy[m,p,a];                                         // nye->alive
          TM_italy[2,1,m,p,a] = 0.0;                                                // dead->nye
          TM_italy[2,2,m,p,a] = 1.0;                                                // dead->dead
          TM_italy[2,3,m,p,a] = 0.0;                                                // dead->alive
          TM_italy[3,1,m,p,a] = 0.0;                                                // alive->nye
          TM_italy[3,2,m,p,a] = 1 - phi_italy[m,p,a];                                     // alive->dead
          TM_italy[3,3,m,p,a] = phi_italy[m,p,a];                                         // alive->alive

          for(s in 1:n_s){

            /* if sex is observed  */
            if(sex_italy[m,a] != -100){
              delta_italy[m,s,p,a] = inv_logit(mu_lodd_delta_italy + // mean
                                              // random effects
                                              r_time_italy[p, 3] + r_area_italy[a,3] + z_dog_italy[3,m,a] * sigma_dog_italy[a,3] +
                                              // predictors
                                              beta_weekend_italy * weekend_italy[s,a,p] + beta_market_italy * market_italy[s,a,p] +
                                              beta_temp_italy * temp_italy_c[s,a,p] + beta_rain_italy * rain_italy[s,a,p] +
                                              beta_sex_italy[2] * sex_italy_c[m,a]
                                          );
             }
            /* if sex is missing */
            if(sex_italy[m,a] == -100){
              delta_italy[m,s,p,a] = (1-kappa_sex_italy) * inv_logit(mu_lodd_delta_italy + // mean
                                                                      // random effects
                                                                      r_time_italy[p, 3] + r_area_italy[a,3] + z_dog_italy[3,m,a] * sigma_dog_italy[a,3] +
                                                                      // predictors
                                                                      beta_weekend_italy * weekend_italy[s,a,p] + beta_market_italy * market_italy[s,a,p] +
                                                                      beta_temp_italy * temp_italy_c[s,a,p] + beta_rain_italy * rain_italy[s,a,p] +
                                                                      // missing
                                                                      beta_sex_italy[2] * -0.5
                                                                    )
                                                          +
                                        kappa_sex_italy  * inv_logit(mu_lodd_delta_italy + // mean
                                                                      // random effects
                                                                      r_time_italy[p, 3] + r_area_italy[a,3] + z_dog_italy[3,m,a] * sigma_dog_italy[a,3] +
                                                                      // predictors
                                                                      beta_weekend_italy * weekend_italy[s,a,p] + beta_market_italy * market_italy[s,a,p] +
                                                                      beta_temp_italy * temp_italy_c[s,a,p] + beta_rain_italy * rain_italy[s,a,p] +
                                                                      // missing
                                                                      beta_sex_italy[2] * 0.5
                                                                    );
             }

              EM_italy[1,1,m,s,p,a] = 1.0;                                    // nye->unobserved
              EM_italy[1,2,m,s,p,a] = 1.0;                                    // dead->unobserved
              EM_italy[1,3,m,s,p,a] = 1.0 - delta_italy[m,s,p,a];             // alive->unobserved
              EM_italy[2,1,m,s,p,a] = 0.0;                                    // nye->observed
              EM_italy[2,2,m,s,p,a] = 0.0;                                    // dead->observed
              EM_italy[2,3,m,s,p,a] = delta_italy[m,s,p,a];                   // alive->observed
          }//s
        }//p
      }//a
    }//m

    /* UKRAINE */
    /* fill the transition & emission matrixes */
    for(m in 1:n_m[2]){
      for(p in 1:n_p){
        for(a in 1:n_a){

          psi_ukr[m,p,a] = inv_logit(mu_lodd_psi_ukr + r_time_ukr[p,1] + r_area_ukr[a,1] +  z_dog_ukr[1,m,a] * sigma_dog_ukr[a,1]);

          /* if sex is observed  */
          if(sex_ukr[m,a] != -100)
            phi_ukr[m,p,a] = inv_logit(mu_lodd_phi_ukr + r_time_ukr[p,2] + r_area_ukr[a,2] + z_dog_ukr[2,m,a] * sigma_dog_ukr[a,2] + beta_sex_ukr[1] * (sex_ukr_c[m,a]));
          /* if sex is missing  */
          if(sex_ukr[m,a] == -100)
            phi_ukr[m,p,a] = (1 - kappa_sex_ukr) * inv_logit(mu_lodd_phi_ukr + r_time_ukr[p,2] + r_area_ukr[a,2] + z_dog_ukr[2,m,a] * sigma_dog_ukr[a,2] + beta_sex_ukr[1]*-0.5) +
                                    kappa_sex_ukr  * inv_logit(mu_lodd_phi_ukr + r_time_ukr[p,2] + r_area_ukr[a,2] + z_dog_ukr[2,m,a] * sigma_dog_ukr[a,2] + beta_sex_ukr[1]* 0.5);


          /* TRANSMISSION MATRIX */
          TM_ukr[1,1,m,p,a] = 1 - psi_ukr[m,p,a];                                     // nye->nye
          TM_ukr[1,2,m,p,a] = 0.0;                                                // nye->dead
          TM_ukr[1,3,m,p,a] = psi_ukr[m,p,a];                                         // nye->alive
          TM_ukr[2,1,m,p,a] = 0.0;                                                // dead->nye
          TM_ukr[2,2,m,p,a] = 1.0;                                                // dead->dead
          TM_ukr[2,3,m,p,a] = 0.0;                                                // dead->alive
          TM_ukr[3,1,m,p,a] = 0.0;                                                // alive->nye
          TM_ukr[3,2,m,p,a] = 1 - phi_ukr[m,p,a];                                     // alive->dead
          TM_ukr[3,3,m,p,a] = phi_ukr[m,p,a];                                         // alive->alive

          for(s in 1:n_s){

            /* if sex is observed  */
            if(sex_ukr[m,a] != -100){
              delta_ukr[m,s,p,a] = inv_logit(mu_lodd_delta_ukr + // mean
                                              // random effects
                                              r_time_ukr[p, 3] + r_area_ukr[a,3] + z_dog_ukr[3,m,a] * sigma_dog_ukr[a,3] +
                                              // predictors
                                              beta_weekend_ukr * weekend_ukr[s,a,p] + beta_market_ukr * market_ukr[s,a,p] +
                                              beta_temp_ukr * temp_ukr_c[s,a,p] + beta_rain_ukr * rain_ukr[s,a,p] +
                                              beta_sex_ukr[2] * sex_ukr_c[m,a]
                                          );
             }
            /* if sex is missing */
            if(sex_ukr[m,a] == -100){
              delta_ukr[m,s,p,a] = (1-kappa_sex_ukr) * inv_logit(mu_lodd_delta_ukr + // mean
                                                                      // random effects
                                                                      r_time_ukr[p, 3] + r_area_ukr[a,3] + z_dog_ukr[3,m,a] * sigma_dog_ukr[a,3] +
                                                                      // predictors
                                                                      beta_weekend_ukr * weekend_ukr[s,a,p] + beta_market_ukr * market_ukr[s,a,p] +
                                                                      beta_temp_ukr * temp_ukr_c[s,a,p] + beta_rain_ukr * rain_ukr[s,a,p] +
                                                                      // missing
                                                                      beta_sex_ukr[2] * -0.5
                                                                    )
                                                          +
                                        kappa_sex_ukr  * inv_logit(mu_lodd_delta_ukr + // mean
                                                                      // random effects
                                                                      r_time_ukr[p, 3] + r_area_ukr[a,3] + z_dog_ukr[3,m,a] * sigma_dog_ukr[a,3] +
                                                                      // predictors
                                                                      beta_weekend_ukr * weekend_ukr[s,a,p] + beta_market_ukr * market_ukr[s,a,p] +
                                                                      beta_temp_ukr * temp_ukr_c[s,a,p] + beta_rain_ukr * rain_ukr[s,a,p] +
                                                                      // missing
                                                                      beta_sex_ukr[2] * 0.5
                                                                    );
              }

              EM_ukr[1,1,m,s,p,a] = 1.0;                                    // nye->unobserved
              EM_ukr[1,2,m,s,p,a] = 1.0;                                    // dead->unobserved
              EM_ukr[1,3,m,s,p,a] = 1.0 - delta_ukr[m,s,p,a];             // alive->unobserved
              EM_ukr[2,1,m,s,p,a] = 0.0;                                    // nye->observed
              EM_ukr[2,2,m,s,p,a] = 0.0;                                    // dead->observed
              EM_ukr[2,3,m,s,p,a] = delta_ukr[m,s,p,a];                   // alive->observed
          }//s
        }//p
      }//a
    }//m
}//transformed parameters

model{
  vector[n_k] alpha_italy[n_p];
  real accumulator_italy[n_k];
  vector[n_k] alpha_ukr[n_p];
  real accumulator_ukr[n_k];

  // priors for ITALY
  mu_lodd_phi_italy ~ normal(0, 1);
  mu_lodd_delta_italy  ~ normal(0, 1);
  mu_lodd_psi_italy    ~ normal(0, 1);
  kappa_sex_italy ~ beta(2, 2);
  for(k in 1:3){
    to_vector(sigma_dog_italy[,k]) ~ normal(0, 1);
    to_vector(z_dog_italy[k]) ~ normal(0, 1);
  }
  etasq_area_italy ~ normal(0, 1);
  etasq_time_italy ~ normal(0, 1);
  rhosq_area_italy ~ normal(0, 1);
  decay_italy ~ normal(0, 5);
  sigma_k_area_italy ~ normal(0, 1);
  sigma_k_time_italy ~ normal(0, 1);
  beta_weekend_italy ~ normal(0, 1);
  beta_market_italy ~ normal(0, 1);
  beta_temp_italy ~ normal(0, 1);
  beta_rain_italy ~ normal(0, 1);
  beta_sex_italy ~ normal(0, 1);
  // ukraine
  mu_lodd_phi_ukr ~ normal(0, 1);
  mu_lodd_delta_ukr ~ normal(0, 1);
  mu_lodd_psi_ukr   ~ normal(0, 1);
  kappa_sex_ukr ~ beta(2, 2);
  for(k in 1:3){
    to_vector(sigma_dog_ukr[,k]) ~ normal(0, 1);
    to_vector(z_dog_ukr[k]) ~ normal(0, 1);
  }
  etasq_area_ukr ~ normal(0, 1);
  etasq_time_ukr ~ normal(0, 1);
  rhosq_area_ukr ~ normal(0, 1);
  decay_ukr ~ normal(0, 5);
  sigma_k_area_ukr ~ normal(0, 1);
  sigma_k_time_ukr ~ normal(0, 1);
  beta_weekend_ukr ~ normal(0, 1);
  beta_market_ukr ~ normal(0, 1);
  beta_temp_ukr ~ normal(0, 1);
  beta_rain_ukr ~ normal(0, 1);
  beta_sex_ukr ~ normal(0, 1);

  // gaussian process MVN
  for(k in 1:3){
    r_time_italy[,k] ~ multi_normal(rep_vector(0, n_p), K_time_italy[k]);
    r_area_italy[,k] ~ multi_normal(rep_vector(0, n_a), K_area_italy[k]);
    r_time_ukr[,k] ~ multi_normal(rep_vector(0, n_p), K_time_ukr[k]);
    r_area_ukr[,k] ~ multi_normal(rep_vector(0, n_a), K_area_ukr[k]);
  }

  /* ITALY LIKELIHOOD */
  for(m in 1:n_m[1]){ // each dog
    for(a in 1:n_a){  // each area

      if(sex_italy[m,a] != -100)
        sex_italy[m,a] ~ bernoulli(kappa_sex_italy);

      // STARTING STATES: assumed to be NYE (state 1) at p = 0
      for(j in 1:n_k){
        alpha_italy[1,j] = TM_italy[1,j,m,1,a];   // nye->nye
        for(s in 1:n_s)
          alpha_italy[1,j] *= EM_italy[ch_italy_add1[m, s, 1, a], j, m, s, 1, a];
        }


      /* FORWARD ALGORITHM */
      // remaining time points
      for(p in 2:n_p){
        //current state
        for(j in 1:n_k){
          // state at p - 1
          for(i in 1:n_k){
            // accumulate possibilites
            accumulator_italy[i] = alpha_italy[p-1,i] * TM_italy[i,j,m,p,a];
            for(s in 1:n_s){
              // compute capture history at p = 1 across secondary periods conditional on state at j
              // if not missing
              if(ch_italy[m,s,p,a] > -1){
                accumulator_italy[i] *= EM_italy[ch_italy_add1[m, s, p, a], j, m, s , p, a];
              }
            }
          }
          // sum possibilites
          alpha_italy[p, j] = sum(accumulator_italy);
        }
      }
      // increment the log probability of possibilites to target
      target += log(sum(alpha_italy[n_p]));
    }//a
  }//m

  /* UKRAINE LIKELIHOOD */
  for(m in 1:n_m[2]){ // each dog
    for(a in 1:n_a){  // each area

      if(sex_ukr[m,a] != -100)
        sex_ukr[m,a] ~ bernoulli(kappa_sex_ukr);

      // STARTING STATES: assumed to be NYE (state 1) at p = 0
      for(j in 1:n_k){
        alpha_ukr[1,j] = TM_ukr[1,j,m,1,a];   // nye->nye
        for(s in 1:n_s)
          alpha_ukr[1,j] *= EM_ukr[ch_ukr_add1[m, s, 1, a], j, m, s, 1, a];
        }


      /* FORWARD ALGORITHM */
      // remaining time points
      for(p in 2:n_p){
        //current state
        for(j in 1:n_k){
          // state at p - 1
          for(i in 1:n_k){
            // accumulate possibilites
            accumulator_ukr[i] = alpha_ukr[p-1,i] * TM_ukr[i,j,m,p,a];
            for(s in 1:n_s){
              // compute capture history at p = 1 across secondary periods conditional on state at j
              // if not missing
              if(ch_ukr[m,s,p,a] > -1){
                accumulator_ukr[i] *= EM_ukr[ch_ukr_add1[m, s, p, a], j, m, s , p, a];
              }
            }
          }
          // sum possibilites
          alpha_ukr[p, j] = sum(accumulator_ukr);
        }
      }
      // increment the log probability of possibilites to target
      target += log(sum(alpha_ukr[n_p]));
    }//a
  }//m
}

generated quantities{
  // italy
  int<lower=1,upper=3> z_italy[n_m[1], n_p, n_a];       // matrix of most possible states
  int<lower=0,upper=1> alive_italy[n_m[1], n_p, n_a];
  int<lower=0,upper=1> ever_alive_italy[n_m[1], n_a];
  int<lower=0> N_hat_italy[n_p,n_a]; //matrix<lower=0>[n_p, n_a] N_hat_italy;
  int<lower=0> superpop_italy[n_a];
  real<lower=0,upper=1> entry_prob_italy[n_p,n_a];
  real<lower=0,upper=1> inclusion_prob_italy[n_a];
  // ukraine
  int<lower=1,upper=3> z_ukr[n_m[2], n_p, n_a];       // matrix of most possible states
  int<lower=0,upper=1> alive_ukr[n_m[2], n_p, n_a];
  int<lower=0,upper=1> ever_alive_ukr[n_m[2], n_a];
  int<lower=0> N_hat_ukr[n_p,n_a]; //matrix<lower=0>[n_p, n_a] N_hat_ukr;
  int<lower=0> superpop_ukr[n_a];
  real<lower=0,upper=1> entry_prob_ukr[n_p,n_a];
  real<lower=0,upper=1> inclusion_prob_ukr[n_a];

  //for(p in 1:n_p){
  //  N_hat_italy[p,1:n_a] = rep_vector(0, n_a)';
  //  N_hat_ukr[p,1:n_a] = rep_vector(0, n_a)';
  //}

  /* italy */
  for(m in 1:n_m[1]){
    for(a in 1:n_a){
      vector[n_k] probs;
      probs = to_vector(TM_italy[1,,m,1,a]);
      z_italy[m,1,a] = categorical_rng(probs);

      for(p in 2:n_p){
        probs = to_vector(TM_italy[z_italy[m,p-1,a],,m,p,a]);
        z_italy[m,p,a] = categorical_rng(probs);
      }

      for(p in 1:n_p){
        if(z_italy[m,p,a] == 3){
          alive_italy[m,p,a] = 1;
        }
        else
          alive_italy[m,p,a] = 0;
        }

        if(sum(alive_italy[m,,a]) > 0){
          ever_alive_italy[m,a] = 1;
        }//if
        else{
          ever_alive_italy[m,a] = 0;
        }
      }//a
    }//m

  for(a in 1:n_a){
    for(p in 1:n_p){
      N_hat_italy[p,a] = sum(alive_italy[1:n_m[1],p,a]);
    }//p
  }//a


  for(a in 1:n_a){
    superpop_italy[a] = sum(ever_alive_italy[1:n_m[1],a]);
  }

  for(a in 1:n_a){
    vector[n_p] q_psi;
    vector[n_p] cprob;
    for(p in 1:n_p){
      q_psi[p] = 1 - mean(to_vector(psi_italy[,p,a]));
    }//p
    cprob[1] = mean(to_vector(psi_italy[,1,a]));
    for(p in 2:n_p){
      cprob[p] = mean(to_vector(psi_italy[,p,a])) * prod(q_psi[1:(p-1)]);
    }//p
    inclusion_prob_italy[a] = sum(cprob[]);
    for(p in 1:n_p)
      entry_prob_italy[p,a] = cprob[p]/inclusion_prob_italy[a];
  }//a

  /* ukr */
  for(m in 1:n_m[2]){
    for(a in 1:n_a){
      vector[n_k] probs;
      probs = to_vector(TM_ukr[1,,m,1,a]);
      z_ukr[m,1,a] = categorical_rng(probs);

      for(p in 2:n_p){
        probs = to_vector(TM_ukr[z_ukr[m,p-1,a],,m,p,a]);
        z_ukr[m,p,a] = categorical_rng(probs);
      }

      for(p in 1:n_p){
        if(z_ukr[m,p,a] == 3){
          alive_ukr[m,p,a] = 1;
        }
        else
          alive_ukr[m,p,a] = 0;
        }//p

        if(sum(alive_ukr[m,,a]) > 0){
          ever_alive_ukr[m,a] = 1;
        }
        else
          ever_alive_ukr[m,a] = 0;
      }//a
    }//m

  for(a in 1:n_a){
    for(p in 1:n_p){
      N_hat_ukr[p,a] = sum(alive_ukr[,p,a]);
    }
  }


  for(a in 1:n_a){
    superpop_ukr[a] = sum(ever_alive_ukr[,a]);
  }

  for(a in 1:n_a){
    vector[n_p] q_psi;
    vector[n_p] cprob;
    for(p in 1:n_p){
      q_psi[p] = 1 - mean(to_vector(psi_ukr[,p,a]));
    }//p
    cprob[1] = mean(to_vector(psi_ukr[,1,a]));
    for(p in 2:n_p){
      cprob[p] = mean(to_vector(psi_ukr[,p,a])) * prod(q_psi[1:(p-1)]);
    }//p
    inclusion_prob_ukr[a] = sum(cprob[]);
    for(p in 1:n_p)
      entry_prob_ukr[p,a] = cprob[p]/inclusion_prob_ukr[a];
  }
}//gq
