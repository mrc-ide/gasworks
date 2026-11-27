## Create a system that has been run to the end of a given week
system_at_week <- function(pars, week = 1, n_particles = 5, seed = 1L,
                           deterministic = FALSE) {
  sys <- dust2::dust_system_create(model, pars, n_particles = n_particles,
                                   seed = seed, deterministic = deterministic,
                                   preserve_particle_dimension = TRUE)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_system_run_to_time(sys, week * 7)
  sys
}

system_state <- function(sys) {
  state <- dust2::dust_system_state(sys)
  rownames(state) <- model_state_names(sys)
  state
}

## Evaluate the comparison function compiled into the model for a single
## week of observed data, optionally replacing some of the model state first.
compare_compiled <- function(sys, observed, prepare, state = NULL) {
  if (!is.null(state)) {
    full_state <- system_state(sys)
    full_state[rownames(state), ] <- state
    dust2::dust_system_set_state(sys, unname(full_state))
  }
  data <- prepare(data.frame(model_week = 1, observed))
  dust2::dust_system_compare_data(sys, as.list(data[names(data) != "time"]))
}

## Simulate, returning a list of arrays with one element per model variable
simulate_unpacked <- function(pars, times, n_particles = 5, seed = 1L) {
  sys <- dust2::dust_system_create(model, pars, n_particles = n_particles,
                                   seed = seed)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_unpack_state(sys, dust2::dust_system_simulate(sys, times))
}
