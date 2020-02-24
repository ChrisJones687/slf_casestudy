## Setup validation comparison across years to see how model improves over time.
library(PoPS)
library(raster)





infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2016_infested.tif"
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
output_frequency <- "year"
movements_file = ""
use_movements = FALSE
num_iterations = 1000
number_of_cores = 30
success_metric = "quantity and configuration"

means <- read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_means.csv")


reproductive_rate <- means[1,1]
percent_natural_dispersal <- means[3,1]
natural_distance_scale <- means[2,1]
anthropogenic_distance_scale <-means[4,1]
natural_kappa <- 0
anthropogenic_kappa <- 0

slf_val_2015_2016 <- PoPS::validate(infected_years_file, num_iterations, number_of_cores,
                                    infected_file, host_file, total_plants_file, 
                                    temp, temperature_coefficient_file, 
                                    precip, precipitation_coefficient_file, 
                                    time_step, reproductive_rate,
                                    season_month_start, season_month_end, 
                                    start_date, end_date, 
                                    use_lethal_temperature, temperature_file,
                                    lethal_temperature, lethal_temperature_month,
                                    mortality_on, mortality_rate, mortality_time_lag, 
                                    management, treatment_dates, treatments_file,
                                    treatment_method,
                                    percent_natural_dispersal,
                                    natural_kernel_type, anthropogenic_kernel_type,
                                    natural_distance_scale, anthropogenic_distance_scale,
                                    natural_dir, natural_kappa, 
                                    anthropogenic_dir, anthropogenic_kappa,
                                    pesticide_duration, pesticide_efficacy,
                                    mask, success_metric, output_frequency,
                                    movements_file, use_movements)




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
output_frequency <- "year"
movements_file = ""
use_movements = FALSE
num_iterations = 1000
number_of_cores = 30
success_metric = "quantity and configuration"

means <- read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_means.csv")


reproductive_rate <- means[1,1]
percent_natural_dispersal <- means[3,1]
natural_distance_scale <- means[2,1]
anthropogenic_distance_scale <-means[4,1]
natural_kappa <- 0
anthropogenic_kappa <- 0


slf_val_2016_2017 <- PoPS::validate(infected_years_file, num_iterations, number_of_cores,
                                    infected_file, host_file, total_plants_file, 
                                    temp, temperature_coefficient_file, 
                                    precip, precipitation_coefficient_file, 
                                    time_step, reproductive_rate,
                                    season_month_start, season_month_end, 
                                    start_date, end_date, 
                                    use_lethal_temperature, temperature_file,
                                    lethal_temperature, lethal_temperature_month,
                                    mortality_on, mortality_rate, mortality_time_lag, 
                                    management, treatment_dates, treatments_file,
                                    treatment_method,
                                    percent_natural_dispersal,
                                    natural_kernel_type, anthropogenic_kernel_type,
                                    natural_distance_scale, anthropogenic_distance_scale,
                                    natural_dir, natural_kappa, 
                                    anthropogenic_dir, anthropogenic_kappa,
                                    pesticide_duration, pesticide_efficacy,
                                    mask, success_metric, output_frequency,
                                    movements_file, use_movements)


infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2018_infested.tif"
number_of_observations <- 33305
prior_number_of_observations <- slf_cal_2016_2017$total_number_of_observations
prior_means <- slf_cal_2016_2017$posterior_means
prior_cov_matrix <- slf_cal_2016_2017$posterior_cov_matrix
params_to_estimate <- c(T, T, T, T)
number_of_generations <- 8
generation_size <- 1000
checks = c(300, 1000000, 1000, 10000)
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
num_iterations = 1000
number_of_cores = 30
success_metric = "quantity and configuration"

means <- read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/2015_2016_means.csv")


reproductive_rate <- means[1,1]
percent_natural_dispersal <- means[3,1]
natural_distance_scale <- means[2,1]
anthropogenic_distance_scale <-means[4,1]
natural_kappa <- 0
anthropogenic_kappa <- 0


slf_val_2017_2018 <- PoPS::validate(infected_years_file, num_iterations, number_of_cores,
                                    infected_file, host_file, total_plants_file, 
                                    temp, temperature_coefficient_file, 
                                    precip, precipitation_coefficient_file, 
                                    time_step, reproductive_rate,
                                    season_month_start, season_month_end, 
                                    start_date, end_date, 
                                    use_lethal_temperature, temperature_file,
                                    lethal_temperature, lethal_temperature_month,
                                    mortality_on, mortality_rate, mortality_time_lag, 
                                    management, treatment_dates, treatments_file,
                                    treatment_method,
                                    percent_natural_dispersal,
                                    natural_kernel_type, anthropogenic_kernel_type,
                                    natural_distance_scale, anthropogenic_distance_scale,
                                    natural_dir, natural_kappa, 
                                    anthropogenic_dir, anthropogenic_kappa,
                                    pesticide_duration, pesticide_efficacy,
                                    mask, success_metric, output_frequency,
                                    movements_file, use_movements)


