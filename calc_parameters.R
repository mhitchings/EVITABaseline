## Code to define prior parameters for the yearly FOI
AG=nrow(case_data) # Number of age groups considered
max_A=49 # Max age considered in the case (note, does not have to be max age. The oldest cases may provide little info about transmission, as we
## are primarily interested in FOI over the period that serosurvey participants were alive)
agegroups = c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")

T_case=ncol(case_data)

# For all cases
incidence <- calc_incidence(case_data, pop_data, 0.1)
yearly_foi_sd <- 0.6
foi_mean <- -0.37
foi_sd <- 0.39

# For confirmed cases
incidence <- calc_incidence(case_data_conf, pop_data, 0.1)
foi_mean_conf <- -1.97
foi_sd_conf <- 1.22

# Lower and upper bound of age classes
lr_bound = c(1,seq(5,45,by=5))
ur_bound = c(seq(4,49,by=5))

## Input parameters to model
T_lambda=12 ## How many years of FOI are we trying to estimate? FOI estimates for years before the earliest birth year in the serosurvey will be highly uncertain
n_serotypes=4 ## Assume that n_serotypes of dengue circulate every year. 

serosurv_year = 2020
