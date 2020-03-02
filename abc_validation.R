library(PoPS)
library(raster)
library(ggplot2)
library(scales)
## set parameters that don't change from year to year
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/tree_of_heaven_0.50.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_6_state_region_psuedo_mercator/total_hosts.tif"
temp <- TRUE
precip <- FALSE
precipitation_coefficient_file <- ""
time_step <- "month"
season_month_start <- 5
season_month_end <- 11
use_lethal_temperature <- TRUE
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
output_frequency <- "year"
movements_file = ""
use_movements = FALSE
num_iterations = 10000
number_of_cores = 30
success_metric = "quantity and configuration"
natural_kappa <- 0
anthropogenic_kappa <- 0
params_to_estimate <- c(T, T, T, T, F, F)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4th: anthropogenic distance

## set data files for individual years
infected_file_2015 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2015_infested.tif"
infected_file_2016 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2016_infested.tif"
infected_file_2017 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2017_infested.tif"
infected_file_2018 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2018_infested.tif"
infected_file_2019 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2019_infested.tif"

# infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2016_infested.tif"
# infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/slf2015_infested.tif"
temperature_coefficient_file_2016 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/temp_coef_2016.tif"
temperature_coefficient_file_2017 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/temp_coef_2017.tif"
temperature_coefficient_file_2018 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/temp_coef_2018.tif"
temperature_coefficient_file_2019 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/temp_coef_2019.tif"

temperature_file_2016 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/crit_temp_2016.tif"
temperature_file_2017 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/crit_temp_2017.tif"
temperature_file_2018 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/crit_temp_2018.tif"
temperature_file_2019 <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/crit_temp_2019.tif"

posterior_means_2016 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2016_means.csv")))[1,1:4]
posterior_cov_matrix_2016 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2016_cov_matrix.csv"))

posterior_means_2017 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2017_means.csv")))[1,1:4]
posterior_cov_matrix_2017 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2017_cov_matrix.csv"))

posterior_means_2018 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2018_means.csv")))[1,1:4]
posterior_cov_matrix_2018 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2018_cov_matrix.csv"))

posterior_means_2019 <- as.matrix(t(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2019_means.csv")))[1,1:4]
posterior_cov_matrix_2019 <- as.matrix(read.csv("H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/slf_data_redone_with_all_data_sources/calibration and assessment/actual_weights/2019_cov_matrix.csv"))

## compare all calibrations to hindcast in 2016
start_date <- '2016-01-01'
end_date <- '2016-12-31'
treatment_dates <- c('2016-12-24')

val_2016_params2016 <- abc_validate(infected_years_file = infected_file_2016, 
                                  num_iterations, 
                                  number_of_cores,
                                  posterior_means = posterior_means_2016,
                                  posterior_cov_matrix= posterior_cov_matrix_2016,
                                  params_to_estimate,
                                  infected_file = infected_file_2015, 
                                  host_file, 
                                  total_plants_file, 
                                  temp, 
                                  temperature_coefficient_file = temperature_coefficient_file_2016, 
                                  precip, 
                                  precipitation_coefficient_file, 
                                  time_step,
                                  season_month_start, 
                                  season_month_end, 
                                  start_date, 
                                  end_date,  
                                  use_lethal_temperature, 
                                  temperature_file = temperature_file_2016,
                                  lethal_temperature, 
                                  lethal_temperature_month,
                                  mortality_on, 
                                  mortality_rate, 
                                  mortality_time_lag, 
                                  management, 
                                  treatment_dates, 
                                  treatments_file,
                                  treatment_method,
                                  natural_kernel_type,
                                  anthropogenic_kernel_type,
                                  natural_dir, 
                                  natural_kappa, 
                                  anthropogenic_dir, 
                                  anthropogenic_kappa, 
                                  pesticide_duration, 
                                  pesticide_efficacy,
                                  mask, 
                                  success_metric, 
                                  output_frequency,
                                  movements_file, 
                                  use_movements)

val_2016_params2017 <- abc_validate(infected_years_file = infected_file_2016, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2015, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2016, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2016,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2016_params2018 <- abc_validate(infected_years_file = infected_file_2016, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2015, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2016, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2016,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2016_params2016$year <- "2016"
val_2016_params2017$year <- "2016"
val_2016_params2018$year <- "2016"

val_2016_params2016$cal_year <- "2016"
val_2016_params2017$cal_year <- "2017"
val_2016_params2018$cal_year <- "2018"

# v_16_16_sd <- apply(val_2016_params2016, 2, sd)
# v_16_16_mean <- apply(val_2016_params2016, 2, mean)
# 
# v_16_17_sd <- apply(val_2016_params2017, 2, sd)
# v_16_17_mean <- apply(val_2016_params2017, 2, mean)
# 
# v_16_18_sd <- apply(val_2016_params2018, 2, sd)
# v_16_18_mean <- apply(val_2016_params2018, 2, mean)

## compare all calibrations to hindcast in 2017
start_date <- '2017-01-01'
end_date <- '2017-12-31'
treatment_dates <- c('2017-12-24')

val_2017_params2016 <- abc_validate(infected_years_file = infected_file_2017, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2016,
                                    posterior_cov_matrix= posterior_cov_matrix_2016,
                                    params_to_estimate,
                                    infected_file = infected_file_2016, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2017, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2017,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2017_params2017 <- abc_validate(infected_years_file = infected_file_2017, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2016, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2017, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2017,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2017_params2018 <- abc_validate(infected_years_file = infected_file_2017, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2016, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2017, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2017,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2017_params2016$year <- "2017"
val_2017_params2017$year <- "2017"
val_2017_params2018$year <- "2017"

val_2017_params2016$cal_year <- "2016"
val_2017_params2017$cal_year <- "2017"
val_2017_params2018$cal_year <- "2018"

## compare all calibrations to hindcast in 2018
start_date <- '2018-01-01'
end_date <- '2018-12-31'
treatment_dates <- c('2018-12-24')

val_2018_params2016 <- abc_validate(infected_years_file = infected_file_2018, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2016,
                                    posterior_cov_matrix= posterior_cov_matrix_2016,
                                    params_to_estimate,
                                    infected_file = infected_file_2017, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2018,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2018_params2017 <- abc_validate(infected_years_file = infected_file_2018, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2017, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2018,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2018_params2018 <- abc_validate(infected_years_file = infected_file_2018, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2017, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2018, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2018,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2018_params2016$year <- "2018"
val_2018_params2017$year <- "2018"
val_2018_params2018$year <- "2018"

val_2018_params2016$cal_year <- "2016"
val_2018_params2017$cal_year <- "2017"
val_2018_params2018$cal_year <- "2018"


## compare all calibrations to hindcast in 2019
start_date <- '2019-01-01'
end_date <- '2019-12-31'
treatment_dates <- c('2019-12-24')

val_2019_params2016 <- abc_validate(infected_years_file = infected_file_2019, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2016,
                                    posterior_cov_matrix= posterior_cov_matrix_2016,
                                    params_to_estimate,
                                    infected_file = infected_file_2018, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2019, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2019,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2019_params2017 <- abc_validate(infected_years_file = infected_file_2019, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2017,
                                    posterior_cov_matrix= posterior_cov_matrix_2017,
                                    params_to_estimate,
                                    infected_file = infected_file_2018, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2019, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2019,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2019_params2018 <- abc_validate(infected_years_file = infected_file_2019, 
                                    num_iterations, 
                                    number_of_cores,
                                    posterior_means = posterior_means_2018,
                                    posterior_cov_matrix= posterior_cov_matrix_2018,
                                    params_to_estimate,
                                    infected_file = infected_file_2018, 
                                    host_file, 
                                    total_plants_file, 
                                    temp, 
                                    temperature_coefficient_file = temperature_coefficient_file_2019, 
                                    precip, 
                                    precipitation_coefficient_file, 
                                    time_step,
                                    season_month_start, 
                                    season_month_end, 
                                    start_date, 
                                    end_date,  
                                    use_lethal_temperature, 
                                    temperature_file = temperature_file_2019,
                                    lethal_temperature, 
                                    lethal_temperature_month,
                                    mortality_on, 
                                    mortality_rate, 
                                    mortality_time_lag, 
                                    management, 
                                    treatment_dates, 
                                    treatments_file,
                                    treatment_method,
                                    natural_kernel_type,
                                    anthropogenic_kernel_type,
                                    natural_dir, 
                                    natural_kappa, 
                                    anthropogenic_dir, 
                                    anthropogenic_kappa, 
                                    pesticide_duration, 
                                    pesticide_efficacy,
                                    mask, 
                                    success_metric, 
                                    output_frequency,
                                    movements_file, 
                                    use_movements)

val_2019_params2016$year <- "2019"
val_2019_params2017$year <- "2019"
val_2019_params2018$year <- "2019"

val_2019_params2016$cal_year <- "2016"
val_2019_params2017$cal_year <- "2017"
val_2019_params2018$cal_year <- "2018"


val_total_params2016 <- val_2016_params2016[,1:10] + val_2017_params2016[,1:10] + val_2018_params2016[,1:10] + val_2019_params2016[,1:10]
val_total_params2016$year <- "All"
val_total_params2016$cal_year <- "2016"

val_total_params2017 <- val_2016_params2017[,1:10] + val_2017_params2017[,1:10] + val_2018_params2017[,1:10] + val_2019_params2017[,1:10]
val_total_params2017$year <- "All"
val_total_params2017$cal_year <- "2017"

val_total_params2018 <- val_2016_params2018[,1:10] + val_2017_params2018[,1:10] + val_2018_params2018[,1:10] + val_2019_params2018[,1:10]
val_total_params2018$year <- "All"
val_total_params2018$cal_year <- "2018"

vals <- rbind(val_2016_params2016, val_2016_params2017, val_2016_params2018, val_2017_params2016, val_2017_params2017, val_2017_params2018, val_2018_params2016, val_2018_params2017, val_2018_params2018, val_2019_params2016, val_2019_params2017, val_2019_params2018)
vals_total <- rbind(val_total_params2016, val_total_params2017, val_total_params2018)
vals_total$configuration_disagreement <- vals_total$configuration_disagreement/4
vals_total$configuration_disagreement_weighted[vals_total$cal_year == "2016"] <- val_2016_params2016$configuration_disagreement * (334/total_infs) + val_2017_params2016$configuration_disagreement * (676/total_infs) + val_2018_params2016$configuration_disagreement * (1740/total_infs) + val_2019_params2016$configuration_disagreement * (4861/total_infs)
vals_total$configuration_disagreement_weighted[vals_total$cal_year == "2017"] <- val_2016_params2017$configuration_disagreement * (334/total_infs) + val_2017_params2017$configuration_disagreement * (676/total_infs) + val_2018_params2017$configuration_disagreement * (1740/total_infs) + val_2019_params2017$configuration_disagreement * (4861/total_infs)
vals_total$configuration_disagreement_weighted[vals_total$cal_year == "2018"] <- val_2016_params2018$configuration_disagreement * (334/total_infs) + val_2017_params2018$configuration_disagreement * (676/total_infs) + val_2018_params2018$configuration_disagreement * (1740/total_infs) + val_2019_params2018$configuration_disagreement * (4861/total_infs)


ggplot(vals, aes(x=year, y=quantity_disagreement, fill=cal_year)) + 
  geom_boxplot()

ggplot(vals, aes(x=year, y=configuration_disagreement, fill=cal_year)) + 
  geom_boxplot()

ggplot(vals_total, aes(x=year, y=quantity_disagreement, fill=cal_year)) + 
  geom_boxplot()

ggplot(vals_total, aes(x=year, y=configuration_disagreement, fill=cal_year)) + 
  geom_boxplot()

ggplot(vals_total, aes(x=year, y=configuration_disagreement_weighted, fill=cal_year)) + 
  geom_boxplot() + xlab("All years combined") + ylab("Configuration Disagreement") + 
  scale_y_continuous(labels = percent) + labs(fill ="Calibration Year")

#check total values
inf_2019 <- raster(infected_file_2019)
sum(inf_2019[inf_2019 > 0] > 0)

inf_2018 <- raster(infected_file_2018)
sum(inf_2018[inf_2018 > 0] > 0)

inf_2017 <- raster(infected_file_2017)
sum(inf_2017[inf_2017 > 0] > 0)

inf_2016 <- raster(infected_file_2016)
sum(inf_2016[inf_2016 > 0] > 0)

inf_2015 <- raster(infected_file_2015)
sum(inf_2015[inf_2015 > 0] > 0)
