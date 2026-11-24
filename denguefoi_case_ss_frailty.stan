//Age-stratified dengue serosurvey and age-stratified case data used to infer 
//annual force of infection with spatial variation

data {
  int <lower=0> N; //the number of individuals in the serosurvey
  int <lower=0,upper=1> Y[N]; //serostatus of each individual
  int firstyear_exposure[N]; // calendar year of first dengue season experienced
  int serosurv_year; // calendar year in which the serosurvey was conducted
  int<lower=2,upper=4>n_serotypes; //how many "endemic" serotypes circulate
  
  int <lower=0> AG; //the number of age classes
  int <lower=0> max_A; //the max age
  int casedata_lastyear; //the latest year for which case data is available
  int casedata_firstyear; // the earliest year for which case data is available
  int <lower=0>	reported_cases[AG, casedata_lastyear-casedata_firstyear+1]; //the number of cases reported at each time point for each age group
  int <lower=0> pop[AG, casedata_lastyear-casedata_firstyear+1]; // the total population at each time point for each age group
  int <lower=0> lr_bound[AG]; // lower age groups bounds
  int <lower=0> ur_bound[AG]; // upper age groups bounds
  real foi_mean;// hyperparameter lambda0: prior mean
  real foi_sd_prior; // hyperparameter lambda0: prior sd
  real reporting_priora;// hyperparameter reporting rate: prior a (beta distribution)
  real reporting_priorb; // hyperparameter reporting rate: prior b (beta distribution)
  int<lower=0,upper=1> vary_reprateratio;// whether to fix the reporting rate ratio of primary vs. secondary cases to 1
  real<lower=0> fixed_reprateratio; //value to fix reprateratio to if fixed
  // if 0, reported cases are only secondary; if 1, reporting rate for primary and secondary is equal
  
  int clust[N]; //cluster of each individual (numeric from 1 to num_clust, must be ordered)
  int <lower=0> num_clust; //number of clusters
  real<lower=0> alpha_priormean;
  real<lower=0> alpha_priorvar;
  
}

transformed data {
  int T_case = casedata_lastyear-casedata_firstyear+1; // number of years for which we have case data
  
  int serosurv_firstyear = min(firstyear_exposure); // earliest birth year of individuals enrolled in the serosurvey
  int serosurv_lastyear = max(firstyear_exposure); // latest birth year of individuals enrolled in the serosurvey
  int T_ss = serosurv_lastyear - serosurv_firstyear + 1; // number of years from which we are estimating lambda from the serosurvey
  
  int lastyear_annualfoi = max(serosurv_lastyear,casedata_lastyear);
  
  int firstyear_foi = casedata_firstyear - max_A;
  
  int T_lambda = lastyear_annualfoi - firstyear_foi + 1; // the number of years for which we will estimate FOI
  
}

parameters {
  
  real <lower=0> invroot_phi; // overdispersion of the case reporting
  real <lower=0,upper=1> reporting_rate; // reporting rate
  real lambda0_logit; // average FOI over time
  real lambdaRE_logit[T_lambda]; // yearly FOI random effect
  real <lower=0, upper=1> alpha; // proportion of first age group experiencing lambda[t-1]
  
  real <lower=0> reprate_pvss_param[vary_reprateratio ? 1:0]; // relative reporting rate for primary vs secondary infections
  
  real<lower=0,upper=5> frailty_alpha; //gamma distribution of hazard by cluster, alpha
  real<lower=0> frailty[num_clust];
  
  real<lower=0> foi_sd; // annual variation in FOI (log-odds scale)
  
}

transformed parameters {
  real <lower=0, upper=1> lambda[max_A + T_case]; // yearly FOI
  real <lower=0> phi;
  real reprate_pvss;
  
  phi = 1/(invroot_phi)^2;

  for (t in 1:(max_A + T_case)) {
    lambda[t] = inv_logit(lambdaRE_logit[t]);
  }
  
  if (vary_reprateratio==1) {
    reprate_pvss = reprate_pvss_param[1];
  } else {
    reprate_pvss = fixed_reprateratio;
  }
  
  // construct individual-level cumulative FOI based on each individual's age and cluster
  real lambda_i[N]; // Lambda for each individual included in the serosurvey
  real exp_sp[N]; // Probability of seropositivity for each individual included in the serosurvey
  real exp_sn[N]; // Probability of being seronegative for each individual included in the serosurvey
  
  real cum_lambda[max_A, T_case]; // Cumulative lambda experienced by residents of BH
  real mono[max_A, T_case]; // Proportion monotypic
  real susc[max_A, T_case]; // Proportion susceptible
  real exp_inc[max_A, T_case]; // Expected incidence
  real exp_reported_cases[AG, T_case]; // Expected reported cases
  real exp_inc_grouped[AG, T_case]; // Expected incidence in age groups

  //Lambda by individual
  for (i in 1:N) {
    // For all children, we can estimate annual FOI from firstyear_exposure to serosurv_year
    lambda_i[i] = 0;
    
    for (y in firstyear_exposure[i]:serosurv_year) {
      
      // a = age is y-serosurv_firstyear+1
      lambda_i[i] = lambda_i[i] + lambda[T_lambda  - (casedata_lastyear - y)];
      
    }

    // Expected probability of seropositivity given individual-level FOI (multiplied by frailty terms)
    exp_sp[i] = 1-exp(-n_serotypes * lambda_i[i] * frailty[clust[i]]) ; 
    exp_sn[i] = exp(-n_serotypes * lambda_i[i] * frailty[clust[i]]) ; 
  }

  // Now get hazard experienced by age cohorts for case data
  
  // Cumulative FOI by age in years (a) and year considered (t)
  cum_lambda[1, 1] = alpha * lambda[T_lambda - (lastyear_annualfoi - casedata_firstyear)];
  for (a in 2:max_A) {
    cum_lambda[a, 1] = sum(lambda[(T_lambda - (lastyear_annualfoi - casedata_firstyear) - a + 1):(T_lambda - (lastyear_annualfoi - casedata_firstyear))]);
  }
  for (t in 2:T_case) {
    cum_lambda[1, t] = alpha * lambda[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1];//first age group born with no past exposure
    for (a in 2:max_A) {
      cum_lambda[a, t] = cum_lambda[a-1, t-1] + lambda[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1];
    }
  }

  // Expected case numbers (representing primary and secondary cases)
  for (t in 1:T_case) {
    for (a in 1:max_A) {
      susc[a, t] = exp(-n_serotypes * cum_lambda[a, t]) ;
      mono[a, t] = n_serotypes * exp(-(n_serotypes-1) * cum_lambda[a, t]) * (1 - exp(-cum_lambda[a, t])) ;
      exp_inc[a, t] = reprate_pvss * n_serotypes * lambda[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1] * susc[a,t] + (n_serotypes-1) * lambda[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1] * mono[a, t];
    }
  }

  // Expected reported case numbers by age groups in data (as expectation of negative binomial for the likelihood of observed case numbers)
  for (t in 1:T_case) {
    for (a_gr in 1:AG) {
      exp_inc_grouped[a_gr, t] = mean(exp_inc[lr_bound[a_gr]:ur_bound[a_gr], t]);
      exp_reported_cases[a_gr, t] = exp_inc_grouped[a_gr, t] * reporting_rate * pop[a_gr, t];
    }
  }
  
}

model {
  
  // Priors
  alpha ~ beta(2, 1);
  
  // Frailty by cluster in serosurvey
  frailty_alpha ~ lognormal(alpha_priormean,alpha_priorvar);
  for (i in 1:num_clust) {
    frailty[i] ~ gamma(frailty_alpha,frailty_alpha);
  }

  reporting_rate ~ beta(reporting_priora,reporting_priorb);
  
  lambda0_logit ~ normal(foi_mean, sqrt(foi_sd_prior ^2 / 2));
  for (i in 1:(max_A+T_case)) {
    lambdaRE_logit[i] ~ normal(lambda0_logit, sqrt(foi_sd ^2 / 2));
  }
  
  invroot_phi ~ normal(0, 100);
  
  for (i in 1:N) {
    //likelihood for seroprevalence data
    target += (Y[i] * log(exp_sp[i]) + (1 - Y[i]) * log(exp_sn[i]));
  }
  
  // likelihood for case data
  for (t in 1:T_case) {
    for (a_gr in 1:AG) {
      reported_cases[a_gr, t] ~ neg_binomial_2(exp_reported_cases[a_gr, t], phi);
    }
  }
  
}

