### Analysis of baseline seroprevalence for EVITA trial
require(dplyr)
require(tidyr)
require(lubridate)
require(ggplot2); theme_set(theme_bw())
require(grid)
require(gridExtra)
require(binom)
require(rstan)
options(mc.cores = parallel::detectCores()) # checks number of cores without having later to specify the cores argument
rstan_options(auto_write = TRUE) # extended packages to use stan


# Serosurvey data (contains de-identified cluster ID, age, and dengue seropositivity)
all_df = read.csv('./Serosurvey_Clean.csv')

pop_data = read.csv('./popdata.csv')
case_data = read.csv('./case_data.csv')
case_data_conf = read.csv('./case_data_labconf.csv')

colnames(case_data) = colnames(case_data_conf) = colnames(pop_data) = 2000:2020
rownames(case_data) = rownames(case_data_conf) = rownames(pop_data) = c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")

source('./calc_parameters.R')

## Input data

## Serosurvey data
## One line per individual
## Analysis code assumes the following variables:
## age - age of child in years
## denpos - dengue serostatus
## clusterid - numerical de-identified cluster ID

## Case data and population data

## This case data is derived from publicly available SINAN data restricted to the municipality of Belo Horizonte
## The input to the model is a matrix with number of columns equal to the time span of the case data
## and number of rows equal to the number of age groups considered
## Each entry of case_data is the number of cases in that year and age group

## The population data (pop_data) is a matrix of the same dimension, with the total population of as in each age group and year as each entry
## We assume population is constant over time. Population data is derived from 2010 Brazilian census

expit = function(x) {1/(1+exp(-x))}
# Function to run the modek
run_regression_model_case_ss = function(compiled_model, # compiled stan model
                                        serosurv_data, # serosurvey data
                                        case_data, # case data (by age group and year)
                                        pop_data, # population data (by age group and year)
                                        lr_bound, # age group lower bounds
                                        ur_bound, # age group upper bounds
                                        foi_mean_prior, # prior mean for annual mean FOI (log-odds scale)
                                        foi_sd_prior, # prior sd for annual mean FOI (log-odds scale)
                                        reporting_priora, # reporting rate: beta distribution prior a
                                        reporting_priorb, # reporting rate: beta distribution prior b
                                        vary_reprateratio, # 0/1: vary reporting rate between primary and secondary infections?
                                        fixed_reprateratio, # if primary/secondary reporting rate ratio is fixed - what is the fixed value?
                                        n_serotypes, # number of circulating serotypes
                                        runmod, # 0/1: whether to fit model
                                        numiter, # number of MCMC iterations
                                        suffix) { # extra file suffix
  
  serosurv_data = serosurv_data %>% arrange(clusterid)
  
  ## Assume no risk in the first year due to maternal immunity
  firstyear_exposure = serosurv_year - serosurv_data$age + 1
  
  # Unique file suffix based on input parameters
  filesuffix = paste0('_ns_',n_serotypes,
                      '_fpm_',round(foi_mean_prior,2),
                      '_fps_',round(foi_sd_prior,2),
                      '_rpm_',round(reporting_priora,2),
                      '_rps_',round(reporting_priorb,2),
                      '_vrr_',vary_reprateratio,
                      '_frr_',fixed_reprateratio,
                      suffix)
  print(filesuffix)
  
  casedata_firstyear = min(as.numeric(colnames(case_data)))
  casedata_lastyear = max(as.numeric(colnames(case_data)))
  
  earliest_yob <- min(firstyear_exposure)
  serosurv_lastyear = max(firstyear_exposure)
  
  lastyear_annualfoi = max(serosurv_lastyear,casedata_lastyear)
  
  # Input data for rstan
  alpha_priormean=0
  alpha_priorvar=1
  
  # Input data for rstan
  data_prov <-
    list(
      N = nrow(serosurv_data), # the number of individuals in the serosurvey
      Y = serosurv_data$denpos, #serostatus of each individual
      firstyear_exposure = firstyear_exposure, # first year of dengue exposure based on age
      serosurv_year = serosurv_year, # calendar year of the serosurvey
      n_serotypes = n_serotypes,
      
      AG=AG,
      max_A=max_A,
      casedata_lastyear=casedata_lastyear,
      casedata_firstyear=casedata_firstyear,
      reported_cases=case_data,
      pop=pop_data,
      lr_bound=lr_bound,
      ur_bound=ur_bound,
      foi_mean=foi_mean_prior,
      foi_sd_prior=foi_sd_prior,
      reporting_priora=reporting_priora,
      reporting_priorb=reporting_priorb,
      vary_reprateratio=vary_reprateratio,
      fixed_reprateratio=fixed_reprateratio,
      
      clust = serosurv_data$clusterid,
      num_clust = max(serosurv_data$clusterid),
      alpha_priormean=alpha_priormean,
      alpha_priorvar=alpha_priorvar
    )
  pars_to_include = c('invroot_phi','reporting_rate',
                      'lambda0_logit',
                      'exp_sp','lambda_i','foi_sd',
                      'reporting_rate','reprate_pvss',
                      'lambda',
                      'frailty_alpha','frailty',
                      'susc','mono','exp_reported_cases','foi_sd')
  pars_to_plot = c('invroot_phi','reporting_rate',
                   'lambda0_logit','reprate_pvss','frailty_alpha','foi_sd')
  
  ### Fit model
  if (runmod) {
    system.time(
      dengue_stanmod <-
        sampling(
          compiled_model,
          data = data_prov,
          chains = 4,
          iter = numiter,
          thin = 1,
          control = list(adapt_delta = 0.80),
          pars=pars_to_include,
          include = T
        )
    )
    dengue_stanmod@stanmodel@dso <- new("cxxdso")
    saveRDS(dengue_stanmod, file = paste0('./dengue_stanmod',filesuffix,'.rds'))
    
  } else {
    dengue_stanmod = readRDS(paste0('./dengue_stanmod',filesuffix,'.rds'))
  }
  
  ## The rest of the code is processing, making plots, etc.
  
  png(paste0('./posteriortrace',filesuffix,'.png'),width=15,height=9,units='in',res=180)
  plot(stan_trace(dengue_stanmod,c(pars_to_plot,'lp__')))
  dev.off()
  
  exp_sp = rstan::extract(dengue_stanmod,pars="exp_sp", inc_warmup=FALSE, permute=T)$exp_sp
  exp_sp_mean <- apply(exp_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
  
  exp_sp_byage = sapply(1:(numiter*2),function(i) {
    serosurv_data$exp_sptemp=exp_sp[i,]
    exp_sp_byage = serosurv_data %>% group_by(age) %>% summarise(esp=mean(exp_sptemp))
    return(exp_sp_byage$esp)
  })
  exp_sp_byage_summ = apply(exp_sp_byage,1,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  exp_sp_bycluster = sapply(1:(numiter*2),function(i) {
    serosurv_data$exp_sptemp=exp_sp[i,]
    exp_sp_bycluster = serosurv_data %>% group_by(clusterid) %>% summarise(esp=mean(exp_sptemp))
    return(exp_sp_bycluster$esp)
  })
  exp_sp_bycluster_summ = apply(exp_sp_bycluster,1,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  lambda = rstan::extract(dengue_stanmod,pars="lambda",inc_warmup=F,permute=T)$lambda
  lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  colnames(lambda)=(casedata_firstyear - max_A):lastyear_annualfoi
  
  bb_data = serosurv_data %>% group_by(age) %>% summarise(n=n(),x=sum(denpos))
  bin_cis = as.data.frame(t(sapply(1:nrow(bb_data),
                                   function(x) {
                                     d=binom.confint(bb_data$x[x],bb_data$n[x])
                                     return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
                                   })))
  colnames(bin_cis)=c("SP","lwr","upr")
  bin_cis$age = bb_data$age
  
  d_forplot = data.frame('Age'=rep(6:11,2),
                         'Model'=rep(c('Observed','Model estimated'),each=nrow(bin_cis)),
                         'SP'=c(bin_cis[,1],
                                exp_sp_byage_summ[2,]),
                         'lwr'=c(bin_cis[,2],exp_sp_byage_summ[1,]),
                         'upr'=c(bin_cis[,3],exp_sp_byage_summ[3,]))
  
  psero=ggplot() + 
    geom_errorbar(data = d_forplot %>% filter(Model=="Observed"),
                  aes(x=Age,ymin=lwr,ymax=upr),width=0) + 
    geom_line(data = d_forplot,aes(x=Age,y=SP,col=factor(Model,levels=c('Observed','Model estimated'))))+
    geom_ribbon(data = d_forplot,aes(x=Age,ymin=lwr,ymax=upr,fill=factor(Model,levels=c('Observed','Model estimated'))),alpha=0.3)+
    geom_point(data = d_forplot, aes(x=Age,y=SP,shape=factor(Model,levels=c('Observed','Model estimated'))))+
    scale_color_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_fill_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_shape_manual(name=NULL,values=c(16,NA),labels=c('Observed','Model estimated'))+
    ylab('Seroprevalence')+ylim(c(0,1))+
    scale_x_continuous(name='Age',breaks=6:11,labels=6:11)+
    theme(axis.text=element_text(size=6),
          axis.title = element_text(size=6),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position.inside = c(0.1,0.89),
          strip.text = element_text(size=6))
  
  ggsave(paste0('./ModelFitBirthYear',filesuffix,'.png'),
         psero,
         device='png',height=6,width=6)
  
  obs_sp_by_cluster = all_df %>% group_by(clusterid) %>% summarise(sp=mean(denpos))
  
  d_forplot_cluster = data.frame('School'=rep(unique(all_df$clusterid),2),
                                 'Model'=rep(c('Observed','Model estimated'),each=length(unique(all_df$clusterid))),
                                 'SP'=c(obs_sp_by_cluster$sp,
                                        exp_sp_bycluster_summ[2,]))
  ggsave(paste0('./ModelFitSeroByCluster',filesuffix,'.png'),
         d_forplot_cluster %>% ggplot() + geom_density(aes(x=SP,colour=Model,fill=Model),alpha=0.4),
         device='png',height=6,width=6)
  
  reprate = rstan::extract(dengue_stanmod,'reporting_rate',inc_warmup=F,permute=T)$reporting_rate
  reprate = quantile(reprate,c(0.025,0.5,0.975))
  
  rel_reprate = rstan::extract(dengue_stanmod,'reprate_pvss',inc_warmup=F,permute=T)$reprate_pvss
  rel_reprate = quantile(rel_reprate,c(0.025,0.5,0.975))
  
  
  
  susc = rstan::extract(dengue_stanmod,'susc',inc_warmup=F,permute=T)$susc
  susc = apply(susc,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
  
  mono = rstan::extract(dengue_stanmod,'mono',inc_warmup=F,permute=T)$mono
  mono = apply(mono,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
  
  repcases = rstan::extract(dengue_stanmod,'exp_reported_cases')$exp_reported_cases
  repcases = apply(repcases,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
  repcases_forplot = as.data.frame(apply(repcases,1,cbind))
  colnames(repcases_forplot)=c("LCI","Median","UCI")
  repcases_forplot$AgeGroup=agegroups
  repcases_forplot$Year=rep(2000:2020,each=10)
  
  repcasestotal_byyear = apply(repcases,c(1,3),sum)
  repcasestotal_byyear = as.data.frame(t(apply(repcasestotal_byyear,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
  colnames(repcasestotal_byyear)=c("LCI","Median","UCI")
  repcasestotal_byyear$Year=2000:2020
  repcasestotal_byyear$Model='Estimated'
  
  obscasestotal = data.frame('Median'=colSums(case_data),
                             'LCI'=NA,
                             'UCI'=NA,
                             'Year'=2000:2020,
                             'Model'="Observed")
  repcasestotal_byyear = rbind(repcasestotal_byyear,obscasestotal)
  
  
  p1=ggplot() + geom_col(data=repcasestotal_byyear %>% filter(Model=="Observed"),aes(x=Year,y=Median),col='grey',alpha=.3)+
    geom_point(data=repcasestotal_byyear %>% filter(Model=="Estimated"),aes(x = Year, y = Median),col = "steelblue") + 
    geom_errorbar(data=repcasestotal_byyear %>% filter(Model=="Estimated"),
                  aes(x=Year,ymin=LCI,ymax=UCI),col = "steelblue")+
    scale_x_continuous(name="Year",breaks=seq(2000,2020,by=2))+ylab("Reported cases")+
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=6))
  
  repcasestotal_byage = apply(repcases,c(1,2),sum)
  repcasestotal_byage = as.data.frame(t(apply(repcasestotal_byage,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
  colnames(repcasestotal_byage)=c("LCI","Median","UCI")
  repcasestotal_byage$AgeGroup=agegroups
  repcasestotal_byage$Model='Estimated'
  
  obscasestotal = data.frame('Median'=rowSums(case_data),
                             'LCI'=NA,
                             'UCI'=NA,
                             'AgeGroup'=agegroups,
                             'Model'="Observed")
  repcasestotal_byage = rbind(repcasestotal_byage,obscasestotal)
  
  repcasestotal_byage$AgeGroup = factor(repcasestotal_byage$AgeGroup,
                                        levels=agegroups)
  
  agegroups_plot = c("1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49")
  p2=ggplot() + geom_col(data=repcasestotal_byage %>% filter(Model=="Observed"),aes(x=AgeGroup,y=Median),col='grey',alpha=.3)+
    geom_point(data=repcasestotal_byage %>% filter(Model=="Estimated"),aes(x = AgeGroup, y = Median),col = "steelblue") + 
    geom_errorbar(data=repcasestotal_byage %>% filter(Model=="Estimated"),
                  aes(x=AgeGroup,ymin=LCI,ymax=UCI),col = "steelblue")+
    scale_x_discrete(name="Age Group",labels=agegroups_plot,guide = guide_axis(n.dodge = 2))+
    ylab("Reported cases")+
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=6))
  ggsave(paste0('./ModelFitSeroprev_Cases',filesuffix,'.png'),
         grid.arrange(psero,grid.arrange(p1,p2,ncol=2),nrow=2,heights=c(4,2)),
         device='png',height=6,width=6,units='in')
  
  foi_table = data.frame('Year'=earliest_yob:serosurv_year)
  foi_table$FOI = c(paste0(round(n_serotypes*lambda[2,colnames(lambda) %in% foi_table$Year],2),' (',
                           round(n_serotypes*lambda[1,colnames(lambda) %in% foi_table$Year],2),',',
                           round(n_serotypes*lambda[3,colnames(lambda) %in% foi_table$Year],2),')')
  )
  
  write.csv(foi_table,paste0('./FOIByYear',filesuffix,'.csv'),row.names=F)
  
  reprate_table = data.frame('Parameter'=c('Secondary infections','Relative Primary vs. Secondary'))
  reprate_table$RepRate = c(paste0(sprintf('%.2f',(reprate[2])),' (',
                                   sprintf('%.2f',(reprate[1])),',',
                                   sprintf('%.2f',(reprate[3])),')'),
                            paste0(sprintf('%.2f',(rel_reprate[2])),' (',
                                   sprintf('%.2f',(rel_reprate[1])),',',
                                   sprintf('%.2f',(rel_reprate[3])),')')
  )
  
  write.csv(reprate_table,paste0('./RepRate',filesuffix,'.csv'),row.names=F)
  
  ## All other parameters
  pars = rstan::extract(dengue_stanmod,c('foi_sd','frailty_alpha','invroot_phi'))
  pars$foi_sd = sqrt(pars$foi_sd^2/2)
  pars$phi = 1/pars$invroot_phi^2
  pars$spatial_sd = sqrt(1/pars$frailty_alpha)
  pars_table = data.frame('Parameter'=c('SD of annual FOI (log-odds scale)',
                                        'SD of spatial RE',
                                        'Overdispersion of case distribution'),
                          'Estimate'= unlist(lapply(pars[c('foi_sd','spatial_sd','phi')],function(x)
                            paste0(sprintf('%.2f',median(x))," (",
                                   sprintf('%.2f',quantile(x,0.025)),",",
                                   sprintf('%.2f',quantile(x,0.975)),")")
                          ))
  )
  
  write.csv(pars_table,paste0('./ParEstsREPhi',filesuffix,'.csv'),row.names=F)
  
  return(dengue_stanmod)
  
}

model_code_re = stanc(file="./denguefoi_case_ss_frailty.stan")
compiled_model_re = stan_model(stanc_ret = model_code_re)

## Example stan fit:
# Confirmed cases and serosurvey
# Reporting rate ratio fixed to 1
denguemod = run_regression_model_case_ss(compiled_model_re,
                                            all_df,
                                            case_data_conf,
                                            pop_data,
                                            lr_bound,
                                            ur_bound,
                                            foi_mean_prior = foi_mean_conf,
                                            foi_sd_prior =foi_sd_conf,
                                            reporting_priora = 1,
                                            reporting_priorb = 1,
                                            vary_reprateratio = 0,
                                            fixed_reprateratio = 1,
                                            n_serotypes,
                                            runmod=T,
                                            numiter=500,
                                            '_confcases')
