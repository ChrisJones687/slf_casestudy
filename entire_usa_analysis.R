## Analysis for the entire USA
library(raster)
library(sp)
library(lubridate)
library(PoPS)


infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/whole_usa/slf_infested.tif"
host_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/whole_usa/host.tif"
total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/whole_usa/total_plants.tif"
temperature_file <- ""
temperature_coefficient_file <- "H:/Shared drives/APHIS  Projects/PoPS/Case Studies/spotted_latternfly/whole_usa/slf_infested.tif"
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 5
season_month_end <- 11
time_step <- "month"
start_date <- '2019-01-01'
end_date <- '2020-12-31'
lethal_temperature <- -35
lethal_temperature_month <- 1
random_seed <- 42
reproductive_rate <- 1.7
treatments_file <- ""
treatment_dates <- c('2019-12-24')
treatment_method <- "ratio"
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
percent_natural_dispersal <- .995
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_distance_scale <- 25
anthropogenic_distance_scale <- 8300
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
random_seed = NULL
output_frequency = "year"
movements_file <- ""
use_movements <- FALSE
num_iterations <- 30
number_of_cores <- 30

# data <- PoPS::pops(infected_file, host_file, total_plants_file, 
#                    temp, temperature_coefficient_file, 
#                    precip, precipitation_coefficient_file, 
#                    time_step, reproductive_rate,
#                    season_month_start, season_month_end, 
#                    start_date, end_date, 
#                    use_lethal_temperature, temperature_file,
#                    lethal_temperature, lethal_temperature_month,
#                    mortality_on, mortality_rate, mortality_time_lag, 
#                    management, treatment_dates, treatments_file,
#                    treatment_method,
#                    percent_natural_dispersal,
#                    natural_kernel_type, anthropogenic_kernel_type,
#                    natural_distance_scale, anthropogenic_distance_scale,
#                    natural_dir, natural_kappa, 
#                    anthropogenic_dir, anthropogenic_kappa,
#                    pesticide_duration, pesticide_efficacy,
#                    random_seed, output_frequency,
#                    movements_file, use_movements)

data <- PoPS::pops_multirun(infected_file, host_file, total_plants_file, 
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
                            num_iterations, number_of_cores,
                            pesticide_duration, pesticide_efficacy,
                            random_seed = NULL, output_frequency,
                            movements_file, use_movements)
