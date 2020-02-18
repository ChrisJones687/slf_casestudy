library(PoPS)
infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2016_infested.tif"
number_of_observations <- 47174
prior_number_of_observations <- 0
prior_means <- c(0,0,0,0)
prior_cov_matrix <- matrix(ncol = 4, nrow = 4, 0)
params_to_estimate <- c(T, T, T, T)
number_of_generations <- 8
generation_size <- 1000
checks = c(180,250000, 1000, 10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2015_infested.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/tree_of_heaven_0.50.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/total_hosts.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/temp_coefficient_slf_2015_2017.tif"
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "month"
season_month_start <- 5
season_month_end <- 11
start_date <- '2014-01-01'
end_date <- '2014-12-31'
use_lethal_temperature <- TRUE
temperature_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/crit_temp_slf_2015_2017.tif"
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2014-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "number of locations and total distance"
output_frequency <- "year"
movements_file = ""
use_movements = FALSE

slf_cal_2015_2016 <- abc_calibration(infected_years_file, 
                                     number_of_observations, prior_number_of_observations,
                                     prior_means, prior_cov_matrix, params_to_estimate,
                                     number_of_generations,
                                     generation_size,
                                     checks,
                                     infected_file, host_file, total_plants_file, 
                                     temp, temperature_coefficient_file, 
                                     precip, precipitation_coefficient_file, 
                                     time_step, 
                                     season_month_start, season_month_end, 
                                     start_date, end_date, 
                                     use_lethal_temperature, temperature_file,
                                     lethal_temperature, lethal_temperature_month,
                                     mortality_on, mortality_rate, mortality_time_lag, 
                                     management, treatment_dates, treatments_file,
                                     treatment_method,
                                     natural_kernel_type, anthropogenic_kernel_type,
                                     natural_dir, natural_kappa, 
                                     anthropogenic_dir, anthropogenic_kappa,
                                     pesticide_duration, pesticide_efficacy,
                                     mask, success_metric, output_frequency,
                                     movements_file, use_movements)

## Write data to files
write.csv(slf_cal_2015_2016$posterior_means, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_means.csv", row.names = FALSE)
write.csv(slf_cal_2015_2016$posterior_cov_matrix, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_cov_matrix.csv", row.names = FALSE)
write.csv(slf_cal_2015_2016$total_number_of_observations, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_total_observations.csv", row.names = FALSE)
write.csv(slf_cal_2015_2016$raw_calibration_data, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_raw_calibration.csv", row.names = FALSE)




infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2017_infested.tif"
number_of_observations <- 33305
prior_number_of_observations <- slf_cal_2015_2016$total_number_of_observations
prior_means <- slf_cal_2015_2016$posterior_means
prior_cov_matrix <- slf_cal_2015_2016$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T)
number_of_generations <- 8
generation_size <- 1000
checks = c(180,250000, 1000, 10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2016_infested.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/tree_of_heaven_0.50.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/total_hosts.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/temp_coefficient_slf_2015_2017.tif"
precip <- FALSEC
precipitation_coefficient_file <- ""
time_step <- "month"
season_month_start <- 5
season_month_end <- 11
start_date <- '2016-01-01'
end_date <- '2016-12-31'
use_lethal_temperature <- TRUE
temperature_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/crit_temp_slf_2015_2017.tif"
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2016-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "number of locations and total distance"
output_frequency <- "year"
movements_file = ""
use_movements = FALSE


slf_cal_2016_2017 <- abc_calibration(infected_years_file, 
                                     number_of_observations, prior_number_of_observations,
                                     prior_means, prior_cov_matrix, params_to_estimate,
                                     number_of_generations,
                                     generation_size,
                                     checks,
                                     infected_file, host_file, total_plants_file, 
                                     temp, temperature_coefficient_file, 
                                     precip, precipitation_coefficient_file, 
                                     time_step, 
                                     season_month_start, season_month_end, 
                                     start_date, end_date, 
                                     use_lethal_temperature, temperature_file,
                                     lethal_temperature, lethal_temperature_month,
                                     mortality_on, mortality_rate, mortality_time_lag, 
                                     management, treatment_dates, treatments_file,
                                     treatment_method,
                                     natural_kernel_type, anthropogenic_kernel_type,
                                     natural_dir, natural_kappa, 
                                     anthropogenic_dir, anthropogenic_kappa,
                                     pesticide_duration, pesticide_efficacy,
                                     mask, success_metric, output_frequency,
                                     movements_file, use_movements)


write.csv(slf_cal_2016_2017$posterior_means, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2016_2017_means.csv", row.names = FALSE)
write.csv(slf_cal_2016_2017$posterior_cov_matrix, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2016_2017_cov_matrix.csv", row.names = FALSE)
write.csv(slf_cal_2016_2017$total_number_of_observations, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2016_2017_total_observations.csv", row.names = FALSE)
write.csv(slf_cal_2016_2017$raw_calibration_data, "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2016_2017_raw_calibration.csv", row.names = FALSE)


load("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/slf_cal.RData")


infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2018_infested.tif"
number_of_observations <- 33305
prior_number_of_observations <- slf_cal_2016_2017$total_number_of_observations
prior_means <- slf_cal_2016_2017$posterior_means
prior_cov_matrix <- slf_cal_2016_2017$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T)
number_of_generations <- 8
generation_size <- 1000
checks = c(300,1000000, 1000, 10000)
infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2017_infested.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/tree_of_heaven_0.50.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/total_hosts.tif"
temp <- TRUE
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/temp_coefficient_slf_2015_2017.tif"
precip <- FALSEC
precipitation_coefficient_file <- ""
time_step <- "month"
season_month_start <- 5
season_month_end <- 11
start_date <- '2017-01-01'
end_date <- '2017-12-31'
use_lethal_temperature <- TRUE
temperature_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/crit_temp_slf_2015_2017.tif"
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2017-12-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "number of locations and total distance"
output_frequency <- "year"
movements_file = ""
use_movements = FALSE


slf_cal_2017_2018 <- abc_calibration(infected_years_file, 
                                     number_of_observations, prior_number_of_observations,
                                     prior_means, prior_cov_matrix, params_to_estimate,
                                     number_of_generations,
                                     generation_size,
                                     checks,
                                     infected_file, host_file, total_plants_file, 
                                     temp, temperature_coefficient_file, 
                                     precip, precipitation_coefficient_file, 
                                     time_step, 
                                     season_month_start, season_month_end, 
                                     start_date, end_date, 
                                     use_lethal_temperature, temperature_file,
                                     lethal_temperature, lethal_temperature_month,
                                     mortality_on, mortality_rate, mortality_time_lag, 
                                     management, treatment_dates, treatments_file,
                                     treatment_method,
                                     natural_kernel_type, anthropogenic_kernel_type,
                                     natural_dir, natural_kappa, 
                                     anthropogenic_dir, anthropogenic_kappa,
                                     pesticide_duration, pesticide_efficacy,
                                     mask, success_metric, output_frequency,
                                     movements_file, use_movements)
