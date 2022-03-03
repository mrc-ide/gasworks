
##' @name helium_age_groups
##' @title Helium model age groups
##' @description Helium model age groups
##' @return A list of named model parameters
##' @export
helium_age_groups <- function() {
  age <- c(seq(0, 75, 5), 100)
  list(age_start = age[-length(age)],
       age_end = age[-1],
       n_group = length(age) - 1)
}

##' @title Prepare particle filter data for the helium model
##' @param data The data set to be used for the particle filter,
##' This is essentially a [data.frame()] with at least columns `model_week`,
##' along with the data used in [helium_compare()].
##' @return A a [data.frame()] with columns `step_start`, `step_end`,
##' along with the data used in [helium_compare()].
##' @importFrom mcstate particle_filter_data
helium_prepare_data <- function(data) {
  time <- "model_week"
  states <- c("scarlet_fever_inc", "igas_inc", "daily_pharyngitis_rate",
              "daily_pharyngitis_rate_04", "daily_pharyngitis_rate_05_14",
              "daily_pharyngitis_rate_15_44", "daily_pharyngitis_rate_45_64",
              "daily_pharyngitis_rate_65_74", "daily_pharyngitis_rate_75",
              "daily_scarlet_fever_rate_04", "daily_scarlet_fever_rate_05_14",
              "daily_scarlet_fever_rate_15_44", "daily_scarlet_fever_rate_45_64",
              "daily_scarlet_fever_rate_65_74", "daily_scarlet_fever_rate_75")
  fitted <- data[, c(time, states)]
  mcstate::particle_filter_data(fitted, time = time, rate = 7)
}
