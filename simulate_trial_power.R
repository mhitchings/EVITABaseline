#############################
## Analysis of dengue serosurvey in Belo Horizonte
## Matt Hitchings 09/08/22
# Parameters
years = 2016
studylength = 3
heterogeneous_immunity = 1

require(dplyr)
require(rstan)
require(tidyverse)
require(data.table)
require(ggplot2); theme_set(theme_classic())
require(iccCounts)

max_A = 49
foi_mean_prior=-1.97
foi_sd_prior=1.22
reporting_priora=1
reporting_priorb=1
vary_reprateratio=0
fixed_reprateratio=1
n_serotypes=4
suffix='_confcases'

sensitivity = 0.7
effect_size=0.5

all_df = read.csv('./Serosurvey_Clean.csv')

##### Function to sample number of infection over trial duration
sample_num_infections = function(cumfoi,annualfoi,intervention,effect_size,n_serotypes,frailty,heterogeneous_immunity) {
  
  p_seroneg = exp(-n_serotypes * cumfoi * (frailty * heterogeneous_immunity + (1 - heterogeneous_immunity) * 1) )
  p_mono = n_serotypes * exp(-(n_serotypes-1) * cumfoi * (frailty * heterogeneous_immunity + (1 - heterogeneous_immunity) * 1)) * 
    (1 - exp(-cumfoi * (frailty * heterogeneous_immunity + (1 - heterogeneous_immunity) * 1)))
  p_multi = 1 - p_seroneg - p_mono
  
  serostatus = sample(0:2,1,prob=c(p_seroneg,p_mono,p_multi))
  
  num_inf = 0
  
  for (y in 1:length(annualfoi)) {
    
    if (serostatus > 1) {break}
    
    foi = (intervention * effect_size + (1 - intervention) * 1) * annualfoi[y] * (n_serotypes - serostatus) * frailty
    
    newinf = rbinom(1,1,1 - exp(-foi))
    
    num_inf = num_inf + newinf
    serostatus = serostatus + newinf
    
  }
  
  return(num_inf)
  
}

estimate_power = function(serosurv_data,
                         studylength,
                         sensitivity,
                         effect_size,
                         years,
                         heterogeneous_immunity,
                         covariates,covariate_names,
                         maternal_immunity,
                         foi_mean_prior,
                         foi_sd_prior,
                         foi_sd_annual,
                         reporting_priora,
                         reporting_priorb,
                         vary_reprateratio,
                         fixed_reprateratio,
                         n_serotypes,
                         suffix) {
  
  filesuffix = paste0('_ns_',n_serotypes,
                      '_fpm_',round(foi_mean_prior,2),
                      '_fps_',round(foi_sd_prior,2),
                      '_rpm_',round(reporting_priora,2),
                      '_rps_',round(reporting_priorb,2),
                      '_vrr_',vary_reprateratio,
                      '_frr_',fixed_reprateratio,
                      suffix)
  
  serosurv_data = serosurv_data %>% arrange(clusterid)
  
  stanmod = readRDS(paste0('./dengue_stanmod',filesuffix,'.rds'))
  
  casedata_firstyear = 2000
  casedata_lastyear = 2020
  
  serosurv_year = 2020
  lastyear_annualfoi = casedata_lastyear
  
  # First, add five children per cluster to get to planned sample size after drop-out
  new_children = data.frame('clusterid'=rep(unique(serosurv_data$clusterid),each=5),
                            'age'=rep(6:10,times=length(unique(serosurv_data$clusterid))))
  
  serosurv_data_aug = serosurv_data %>% dplyr::select(clusterid,age) %>% 
    bind_rows(new_children) %>% arrange(clusterid) ## 2,769 children
  
  lambda = rstan::extract(stanmod,pars="lambda",inc_warmup=F,permute=T)$lambda
  dimnames(lambda)[[2]]=(casedata_firstyear - max_A):lastyear_annualfoi
  
  frailty_terms = rstan::extract(stanmod,pars='frailty',inc_warmup=F,permute=T)$frailty
  dimnames(frailty_terms)[[2]] = unique(serosurv_data$clusterid)
  
  ### Go through all years from 2000 to 2016 and estimate cumulative incidence
  ### over the trial period for a cohort of 6-11-year-olds
  
  ### First, need to get serological profile at baseline
  ### Then, get cumulative FOI over trial period
  ### Then, get CI
  
  ### Need to do this for each posterior draw to get uncertainty
  posterior_length = dim(lambda)[1]
  
  numyears = length(years)
  
  trial_cis_control = matrix(NA,nrow=numyears,ncol=posterior_length)
  trial_irrs = matrix(NA,nrow=numyears,ncol=posterior_length)
  trial_powers = matrix(NA,nrow=numyears,ncol=posterior_length)
  trial_iccs = matrix(NA,nrow=numyears,ncol=posterior_length)
  
  frailty_df = data.frame('clusterid'=unique(serosurv_data$clusterid))
  numclusters = nrow(frailty_df)
  clusternames = unique(serosurv_data$clusterid)
  
  cumfoi_byage = data.frame('age'=6:11)
  
  S = 100
  
  for (year_index in 1:numyears) {
    
    year_start = years[year_index]
    
    tic = Sys.time()
    
    for (draw in 1:posterior_length) {
      
      annualfoi = lambda[draw,as.numeric(names(lambda[draw,])) %in% (year_start):(year_start+studylength-1)]
      
      frailty_df = frailty_df %>% mutate(frailty=frailty_terms[draw,])
      
      serosurv_data_temp = serosurv_data_aug %>% left_join(frailty_df,by='clusterid')
      
      cumfoi_byage$cumfoi = sapply(cumfoi_byage$age,function(x) sum(lambda[draw,as.numeric(names(lambda[draw,])) %in% (year_start-x+1):(year_start-1)]))
      
      serosurv_data_temp = serosurv_data_temp %>% left_join(cumfoi_byage,by='age')
      
      ## Now simulate each trial
      pvals = rep(NA,S)
      irrs = rep(NA,S)
      cis = rep(NA,S)
      iccs = rep(NA,S)
      
      for (sim in 1:S) {
        
        ## Assign half of clusters to get intervention
        serosurv_data_temp$intervention = 0
        serosurv_data_temp$intervention[serosurv_data_temp$clusterid %in% sample(clusternames,numclusters/2,replace=F)] = 1
      
        ## Simulate number of infections for each individual
        serosurv_data_temp = serosurv_data_temp %>% rowwise() %>% 
          mutate(numinf = sample_num_infections(cumfoi,annualfoi,intervention,effect_size,n_serotypes,frailty,heterogeneous_immunity),
                 numdetectedinf = rbinom(1,numinf,prob = sensitivity)) 
        
        iccs[sim] = icc_counts(serosurv_data_temp,y='numinf',id='clusterid')$ICC[1]
        
        serosurv_data_bycluster = serosurv_data_temp %>% 
          group_by(intervention,clusterid) %>% 
          summarise(numevents = sum(numdetectedinf),
                    logpersontime = log(n() * studylength),
                    .groups = 'drop_last')
        
        cis[sim] = sum(serosurv_data_bycluster$numevents[serosurv_data_bycluster$intervention==0]) / 
          sum(serosurv_data_temp$intervention==0)
        
        
        glm.qpois <- glm(numevents ~ offset(logpersontime) + intervention, family = quasipoisson(),data = serosurv_data_bycluster)
        
        irrs[sim] = unname(exp(coef(glm.qpois)[2]))
        pvals[sim] = summary(glm.qpois)$coef[2,4]
        
      }
      
      trial_powers[year_index,draw] = mean(pvals<=0.05)
      trial_irrs[year_index,draw] = mean(irrs)
      trial_cis_control[year_index,draw] = mean(cis)
      trial_iccs[year_index,draw] = mean(iccs,na.rm=T)
      
      cat(draw,'\n')
      print(Sys.time() - tic)
      tic = Sys.time()
      
    }
    
    
    
  }
  
  write.csv(trial_powers,paste0('./trialpower_year_',years,'_studylength_',studylength,'_hetimm_',heterogeneous_immunity,filesuffix,'.csv'),row.names=F)
  write.csv(trial_irrs,paste0('./trialirr_year_',years,'_studylength_',studylength,'_hetimm_',heterogeneous_immunity,filesuffix,'.csv'),row.names=F)
  write.csv(trial_cis_control,paste0('./trialcicontrol_year_',years,'_studylength_',studylength,'_hetimm_',heterogeneous_immunity,filesuffix,'.csv'),row.names=F)
  write.csv(trial_iccs,paste0('./trialICCs_year_',years,'_studylength_',studylength,'_hetimm_',heterogeneous_immunity,filesuffix,'.csv'),row.names=F)
  
}

estimate_power(all_df,
               studylength,
               sensitivity,
               effect_size,
               years,
               heterogeneous_immunity,
               foi_mean_prior,
               foi_sd_prior,
               reporting_priora,
               reporting_priorb,
               vary_reprateratio,
               fixed_reprateratio,
               n_serotypes,
               suffix)
