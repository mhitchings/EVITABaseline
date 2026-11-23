# EVITABaseline
Code and data to estimate dengue transmission using baseline serological data from the EVITA Dengue Trial population in Belo Horizonte, Brazil

The following input data are provided:
- Serosurvey_Clean.csv: individual-level, anonymized data including age, de-identified study cluster ID, and dengue seropositivity of 2,479 included individuals enrolled into the EVITA Dengue Trial in 2019
- case_data.csv: age- and year-stratified reported dengue cases in Belo Horizonte
- case_data_labconf.csv: age- and year-stratified reported, lab-confirmed dengue cases in Belo Horizonte
- popdata.csv: population by age of Belo Horizonte from 2010 IBGE census in Brazil

The primary R script to run the catalytic model to estimate annual force of infection from serosurvey and case data is evita_baseline_serology_analysis.R. This script sources calc_parameters and denguefoi_case_ss_frailty.stan, which contains code for the stan model. An example model is fit using lab-confirmed cases only.

The second R script runs the power simulations presented in the paper. These simulations require a stan output to inform annual FOI and variation at the cluster level. This script estimates trial power, cumulative incidence, and ICC for trials starting in given calendar years, with given trial length.
