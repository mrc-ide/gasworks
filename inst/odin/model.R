## Stochastic model of group A Streptococcus transmission.
##
## Time and durations are measured in days. Transmission and demographic rates
## are given per week (as in the data, which are weekly), so we convert the time
## step into weeks here:
dt_week <- dt / 7
# i: number of age groups, j: number of Erlang compartments
n_group <- parameter(1, type = "integer", constant = TRUE)

## Daily time-varying inputs; values are held at their last entry after the
## end of the series
alpha_t <- alpha[min(as.integer(time) + 1, length(alpha))]
## Probability of reporting a scarlet fever case (time-varying)
q_F_t <- q_F[min(as.integer(time) + 1, length(q_F))]

## Core equations for transitions between compartments:
update(U[])    <- U[i]    + dem_U[i]    + gas_U[i]
update(A[, ])  <- A[i, j] + dem_A[i, j] + gas_A[i, j]
update(E[, ])  <- E[i, j] + dem_E[i, j] + gas_E[i, j]
update(S[, ])  <- S[i, j] + dem_S[i, j] + gas_S[i, j]
update(P[, ])  <- P[i, j] + dem_P[i, j] + gas_P[i, j]
update(F[, ])  <- F[i, j] + dem_F[i, j] + gas_F[i, j]
update(R[, ])  <- R[i, j] + dem_R[i, j] + gas_R[i, j]
update(N[])    <- N[i] + dem_N[i]

n_Nx[] <- n_Ux[i] + sum(n_Ex[i, ]) + sum(n_Ax[i, ]) + sum(n_Sx[i, ]) +
  sum(n_Px[i, ]) + sum(n_Fx[i, ]) + sum(n_Rx[i, ])

dem_N[] <- dem_U[i] + sum(dem_E[i, ]) + sum(dem_A[i, ]) + sum(dem_S[i, ]) +
  sum(dem_P[i, ]) + sum(dem_F[i, ]) + sum(dem_R[i, ])

## Output prevalence
update(prev_A[]) <- sum(A[i, ]) / N[i]
update(prev_R[]) <- sum(R[i, ]) / N[i]

## Output incidence flows:
update(births_inc) <- births_inc + sum(n_xU[])

update(net_leavers_inc) <- net_leavers_inc + sum(n_Nx[])

update(infections_inc) <- infections_inc + sum(n_UE[]) + sum(n_UA[])

new_gas_pharyngitis_inc_by_group[] <- gas_pharyngitis_inc_by_group[i] + n_ES[i]
update(gas_pharyngitis_inc_by_group[]) <- new_gas_pharyngitis_inc_by_group[i]

new_scarlet_fever_inc_by_group[] <- scarlet_fever_inc_by_group[i] + n_PF[i]
update(scarlet_fever_inc_by_group[]) <- new_scarlet_fever_inc_by_group[i]

# output all-group incidence flows
update(gas_pharyngitis_inc) <- sum(new_gas_pharyngitis_inc_by_group[])
update(scarlet_fever_inc) <- sum(new_scarlet_fever_inc_by_group[])
update(scarlet_fever_cases) <- round(sum(new_scarlet_fever_inc_by_group[]) * q_F_t)
update(igas_inc) <- igas_inc + sum(n_I[])

## Output daily incidence rate per 100,000 population
 w[] <- N[i] / 1e5 * 7
update(daily_gas_pharyngitis_rate) <-
sum(new_gas_pharyngitis_inc_by_group[]) * p_T / sum(w[])
update(daily_scarlet_fever_rate) <-
  sum(new_scarlet_fever_inc_by_group[]) / sum(w[])

## When using the age-structured model (helium) output incidence rates for age
## groups used by UK GP surveillance:
##  04,  5-14, 15-44,   45-64,   65-74, 75+
## [1], [2:3], [4:9], [10:13], [14:15], [16]

## Output daily GAS pharyngitis rates by group
update(daily_gas_pharyngitis_rate_04) <- (
  if (n_group == 16)
    new_gas_pharyngitis_inc_by_group[1] * p_T / w[1]
  else 0)
update(daily_gas_pharyngitis_rate_05_14) <- (
  if (n_group == 16)
    sum(new_gas_pharyngitis_inc_by_group[2:3]) * p_T / sum(w[2:3])
  else 0)
update(daily_gas_pharyngitis_rate_15_44) <- (
  if (n_group == 16)
    sum(new_gas_pharyngitis_inc_by_group[4:9]) * p_T / sum(w[4:9])
  else 0)
update(daily_gas_pharyngitis_rate_45_64) <- (
  if (n_group == 16)
    sum(new_gas_pharyngitis_inc_by_group[10:13]) * p_T / sum(w[10:13])
  else 0)
update(daily_gas_pharyngitis_rate_65_74) <- (
  if (n_group == 16)
    sum(new_gas_pharyngitis_inc_by_group[14:15]) * p_T / sum(w[14:15])
  else 0)
update(daily_gas_pharyngitis_rate_75) <- (
  if (n_group == 16)
    new_gas_pharyngitis_inc_by_group[16] * p_T / w[16]
  else 0)

## Output daily scarlet fever rates by group
update(daily_scarlet_fever_rate_04) <- (
  if (n_group == 16) new_scarlet_fever_inc_by_group[1] / w[1]
  else 0)
update(daily_scarlet_fever_rate_05_14) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[2:3]) / sum(w[2:3])
  else 0)
update(daily_scarlet_fever_rate_15_44) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[4:9]) / sum(w[4:9])
  else 0)
update(daily_scarlet_fever_rate_45_64) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[10:13]) / sum(w[10:13])
  else 0)
update(daily_scarlet_fever_rate_65_74) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[14:15]) / sum(w[14:15])
  else 0)
update(daily_scarlet_fever_rate_75) <- (
  if (n_group == 16) new_scarlet_fever_inc_by_group[16] / w[16]
  else 0)

## Output weekly scarlet fever cases by group
update(scarlet_fever_inc_04) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[1]) else 0)
update(scarlet_fever_inc_05_14) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[2:3]) else 0)
update(scarlet_fever_inc_15_44) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[4:9]) else 0)
update(scarlet_fever_inc_45_64) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[10:13]) else 0)
update(scarlet_fever_inc_65_74) <- (
  if (n_group == 16) sum(new_scarlet_fever_inc_by_group[14:15]) else 0)
update(scarlet_fever_inc_75) <- (
  if (n_group == 16) new_scarlet_fever_inc_by_group[16] else 0)

# weight phi_S by population size
weighted_phi_S[] <- phi_S[i] * N[i]

update(etiologic_fraction) <- sum(weighted_phi_S[]) / sum(N[])

## Output etiologic fraction by age group
update(etiologic_fraction_04) <- (
  if (n_group == 16) weighted_phi_S[1] / N[1] else 0)
update(etiologic_fraction_05_14) <- (
  if (n_group == 16) sum(weighted_phi_S[2:3]) / sum(N[2:3]) else 0)
update(etiologic_fraction_15_44) <- (
  if (n_group == 16) sum(weighted_phi_S[4:9]) / sum(N[4:9]) else 0)
update(etiologic_fraction_45_64) <- (
  if (n_group == 16) sum(weighted_phi_S[10:13]) / sum(N[10:13]) else 0)
update(etiologic_fraction_65_74) <- (
  if (n_group == 16) sum(weighted_phi_S[14:15]) / sum(N[14:15]) else 0)
update(etiologic_fraction_75) <- (
  if (n_group == 16) weighted_phi_S[16] / N[16] else 0)

## Force of infection
update(beta_t) <- beta * (1 + sigma * cos(2 * pi * (t0 + time - t_s) / 365.25))
## (there is no transmission from an age group with no population)
lambda[, ] <- if (N[j] > 0) beta_t * m[i, j] *
  (sum(A[j, ]) * theta_A + sum(S[j, ]) + sum(P[j, ])) / N[j] else 0
foi[] <- sum(lambda[i, ])

## Total rates of transmission out of each (sub)compartment. The force of
## infection is a rate per week; durations (delta) are in days, so the rates
## of progression are per day.
r_U[]  <- foi[i]
r_A[]  <- k_A / delta_A
r_E[]  <- k_E / delta_E
r_S[]  <- k_S / delta_S
r_P[]  <- k_P / delta_P
r_F[]  <- k_F / delta_F
r_R[]  <- k_R / delta_R
# rate of developing iGAS
r_I[]  <- foi[i] * p_I

## Calculate number of births
n_xU[1] <- round(alpha_t * dt_week)

## Calculate net number of leavers from each compartment - deterministic
## note that this can be negative when there is a net population increase due
## to immigration > emigration + death. entrants are distributed across disease
## states in proportion to the general population.
n_Ux[]   <- round(U[i] * omega[i] * dt_week)
n_Ax[, ] <- round(A[i, j] * omega[i] * dt_week)
n_Ex[, ] <- round(E[i, j] * omega[i] * dt_week)
n_Sx[, ] <- round(S[i, j] * omega[i] * dt_week)
n_Px[, ] <- round(P[i, j] * omega[i] * dt_week)
n_Fx[, ] <- round(F[i, j] * omega[i] * dt_week)
n_Rx[, ] <- round(R[i, j] * omega[i] * dt_week)

## calculate net aging
n_Ui[] <- (if (i > 1) U[i - 1] else 0) - (if (i < n_group) U[i] else 0)
n_Ai[, ] <- (if (i > 1) A[i - 1, j] else 0) - (if (i < n_group) A[i, j] else 0)
n_Ei[, ] <- (if (i > 1) E[i - 1, j] else 0) - (if (i < n_group) E[i, j] else 0)
n_Si[, ] <- (if (i > 1) S[i - 1, j] else 0) - (if (i < n_group) S[i, j] else 0)
n_Pi[, ] <- (if (i > 1) P[i - 1, j] else 0) - (if (i < n_group) P[i, j] else 0)
n_Fi[, ] <- (if (i > 1) F[i - 1, j] else 0) - (if (i < n_group) F[i, j] else 0)
n_Ri[, ] <- (if (i > 1) R[i - 1, j] else 0) - (if (i < n_group) R[i, j] else 0)

## Calculate all demographic changes
dem_U[]   <- n_xU[i] + round(n_Ui[i] * r_age * dt_week) - n_Ux[i]
dem_A[, ] <- round(n_Ai[i, j] * r_age * dt_week) - n_Ax[i, j]
dem_E[, ] <- round(n_Ei[i, j] * r_age * dt_week) - n_Ex[i, j]
dem_S[, ] <- round(n_Si[i, j] * r_age * dt_week) - n_Sx[i, j]
dem_P[, ] <- round(n_Pi[i, j] * r_age * dt_week) - n_Px[i, j]
dem_F[, ] <- round(n_Fi[i, j] * r_age * dt_week) - n_Fx[i, j]
dem_R[, ] <- round(n_Ri[i, j] * r_age * dt_week) - n_Rx[i, j]

## Draws from binomial distributions for numbers leaving each compartment
## all demographic transitions are done first
n_U[] <- Binomial(U[i] + dem_U[i], 1 - exp(-r_U[i] * dt_week))
n_A[, ] <- Binomial(A[i, j] + dem_A[i, j], 1 - exp(-r_A[i] * dt))
n_E[, ] <- Binomial(E[i, j] + dem_E[i, j], 1 - exp(-r_E[i] * dt))
n_S[, ] <- Binomial(S[i, j] + dem_S[i, j], 1 - exp(-r_S[i] * dt))
n_P[, ] <- Binomial(P[i, j] + dem_P[i, j], 1 - exp(-r_P[i] * dt))
n_F[, ] <- Binomial(F[i, j] + dem_F[i, j], 1 - exp(-r_F[i] * dt))
n_R[, ] <- Binomial(R[i, j] + dem_R[i, j], 1 - exp(-r_R[i] * dt))

# Number developing iGAS, applies to whole population
n_I[] <- Binomial(N[i] + dem_N[i], 1 - exp(-r_I[i] * dt_week))

## Draw the numbers of transitions between compartments
n_UE[] <- Binomial(n_U[i], p_S)
n_UA[] <- n_U[i] - n_UE[i]
n_AR[] <- Binomial(n_A[i, k_A], p_R)
n_AU[] <- n_A[i, k_A] - n_AR[i]
n_EP[] <- Binomial(n_E[i, k_E], p_F)
n_ES[] <- n_E[i, k_E] - n_EP[i]
n_SR[] <- n_S[i, k_S]
n_PF[] <- n_P[i, k_P]
n_FR[] <- n_F[i, k_F]
n_RU[] <- n_R[i, k_R]

## GAS transitions
gas_U[]   <- n_AU[i] + n_RU[i] - n_U[i]
gas_E[, ] <- (if (j == 1) n_UE[i] else n_E[i, j - 1]) - n_E[i, j]
gas_A[, ] <- (if (j == 1) n_UA[i] else n_A[i, j - 1]) - n_A[i, j]
gas_S[, ] <- (if (j == 1) n_ES[i] else n_S[i, j - 1]) - n_S[i, j]
gas_P[, ] <- (if (j == 1) n_EP[i] else n_P[i, j - 1]) - n_P[i, j]
gas_F[, ] <- (if (j == 1) n_PF[i] else n_F[i, j - 1]) - n_F[i, j]
gas_R[, ] <- (if (j == 1) n_AR[i] + n_SR[i] + n_FR[i]
              else n_R[i, j - 1]) - n_R[i, j]

## Initial states:
initial(U[]) <- U0[i]
initial(A[, ]) <- A0[i, j]
initial(E[, ]) <- E0[i, j]
initial(S[, ]) <- S0[i, j]
initial(P[, ]) <- P0[i, j]
initial(F[, ]) <- F0[i, j]
initial(R[, ]) <- R0[i, j]
initial(N[]) <- U0[i] + sum(A0[i, ]) + sum(E0[i, ]) + sum(S0[i, ]) +
  sum(P0[i, ]) + sum(F0[i, ]) + sum(R0[i, ])
initial(prev_A[]) <- sum(A0[i, ]) /
  (U0[i] + sum(A0[i, ]) + sum(E0[i, ]) + sum(S0[i, ]) + sum(P0[i, ]) +
     sum(F0[i, ]) + sum(R0[i, ]))
initial(prev_R[]) <- sum(R0[i, ]) /
  (U0[i] + sum(A0[i, ]) + sum(E0[i, ]) + sum(S0[i, ]) + sum(P0[i, ]) +
     sum(F0[i, ]) + sum(R0[i, ]))

initial(gas_pharyngitis_inc_by_group[], zero_every = 7) <- 0
initial(scarlet_fever_inc_by_group[], zero_every = 7) <- 0
initial(infections_inc, zero_every = 7) <- 0
initial(gas_pharyngitis_inc) <- 0
initial(scarlet_fever_inc) <- 0
initial(scarlet_fever_cases) <- 0
initial(igas_inc, zero_every = 7) <- 0
initial(births_inc, zero_every = 7) <- 0
initial(net_leavers_inc, zero_every = 7) <- 0
initial(beta_t) <- 0
initial(daily_gas_pharyngitis_rate) <- 0
initial(daily_scarlet_fever_rate) <- 0
initial(etiologic_fraction) <- 0

initial(daily_gas_pharyngitis_rate_04)   <-  0
initial(daily_gas_pharyngitis_rate_05_14) <- 0
initial(daily_gas_pharyngitis_rate_15_44) <- 0
initial(daily_gas_pharyngitis_rate_45_64) <- 0
initial(daily_gas_pharyngitis_rate_65_74) <- 0
initial(daily_gas_pharyngitis_rate_75)    <- 0

initial(daily_scarlet_fever_rate_04)   <-  0
initial(daily_scarlet_fever_rate_05_14) <- 0
initial(daily_scarlet_fever_rate_15_44) <- 0
initial(daily_scarlet_fever_rate_45_64) <- 0
initial(daily_scarlet_fever_rate_65_74) <- 0
initial(daily_scarlet_fever_rate_75)    <- 0

initial(scarlet_fever_inc_04)    <- 0
initial(scarlet_fever_inc_05_14) <- 0
initial(scarlet_fever_inc_15_44) <- 0
initial(scarlet_fever_inc_45_64) <- 0
initial(scarlet_fever_inc_65_74) <- 0
initial(scarlet_fever_inc_75)    <- 0

initial(etiologic_fraction_04)    <- 0
initial(etiologic_fraction_05_14) <- 0
initial(etiologic_fraction_15_44) <- 0
initial(etiologic_fraction_45_64) <- 0
initial(etiologic_fraction_65_74) <- 0
initial(etiologic_fraction_75)    <- 0

## User defined parameters - default in parentheses:
## Initial number in each state
U0 <- parameter()
A0 <- parameter()
E0 <- parameter()
S0 <- parameter()
P0 <- parameter()
F0 <- parameter()
R0 <- parameter()

beta <- parameter()
m <- parameter()
sigma <- parameter()
t0 <- parameter(0) # day of model initialisation
t_s <- parameter()
p_S <- parameter()
p_R <- parameter()
p_I <- parameter()
p_F <- parameter()
p_T <- parameter()
delta_A <- parameter()
delta_E <- parameter()
delta_S <- parameter()
delta_P <- parameter()
delta_F <- parameter()
delta_R <- parameter()
k_A <- parameter(type = "integer", constant = TRUE)
k_E <- parameter(type = "integer", constant = TRUE)
k_S <- parameter(type = "integer", constant = TRUE)
k_P <- parameter(type = "integer", constant = TRUE)
k_F <- parameter(type = "integer", constant = TRUE)
k_R <- parameter(type = "integer", constant = TRUE)
theta_A <- parameter()

alpha <- parameter() # time-varying number of births
dim(alpha) <- parameter(rank = 1)
omega <- parameter()
r_age <- parameter(0) # rate of aging - determined by group size

q_F <- parameter()
dim(q_F) <- parameter(rank = 1)

phi_S <- parameter()

## Object dimensions

dim(m)      <- c(n_group, n_group)
dim(lambda) <- c(n_group, n_group)
dim(foi)   <- n_group
dim(omega) <- n_group
dim(phi_S) <- n_group

dim(U)  <- n_group
dim(A)  <- c(n_group, k_A)
dim(E)  <- c(n_group, k_E)
dim(S)  <- c(n_group, k_S)
dim(P)  <- c(n_group, k_P)
dim(F)  <- c(n_group, k_F)
dim(R)  <- c(n_group, k_R)
dim(N)  <- n_group

dim(prev_A) <- n_group
dim(prev_R) <- n_group

dim(U0)  <- n_group
dim(A0)  <- c(n_group, k_A)
dim(E0)  <- c(n_group, k_E)
dim(S0)  <- c(n_group, k_S)
dim(P0)  <- c(n_group, k_P)
dim(F0)  <- c(n_group, k_F)
dim(R0)  <- c(n_group, k_R)

dim(r_U) <- n_group
dim(r_A) <- n_group
dim(r_E) <- n_group
dim(r_S) <- n_group
dim(r_P) <- n_group
dim(r_F) <- n_group
dim(r_R) <- n_group
dim(r_I) <- n_group

dim(n_xU) <- n_group

dim(n_U) <- n_group
dim(n_A) <- c(n_group, k_A)
dim(n_E) <- c(n_group, k_E)
dim(n_S) <- c(n_group, k_S)
dim(n_P) <- c(n_group, k_P)
dim(n_F) <- c(n_group, k_F)
dim(n_R) <- c(n_group, k_R)
dim(n_I) <- n_group

dim(n_Ux) <- n_group
dim(n_Ax) <- c(n_group, k_A)
dim(n_Ex) <- c(n_group, k_E)
dim(n_Sx) <- c(n_group, k_S)
dim(n_Px) <- c(n_group, k_P)
dim(n_Fx) <- c(n_group, k_F)
dim(n_Rx) <- c(n_group, k_R)
dim(n_Nx) <- n_group

dim(n_Ui) <- n_group
dim(n_Ai) <- c(n_group, k_A)
dim(n_Ei) <- c(n_group, k_E)
dim(n_Si) <- c(n_group, k_S)
dim(n_Pi) <- c(n_group, k_P)
dim(n_Fi) <- c(n_group, k_F)
dim(n_Ri) <- c(n_group, k_R)

dim(dem_U) <- n_group
dim(dem_A) <- c(n_group, k_A)
dim(dem_E) <- c(n_group, k_E)
dim(dem_S) <- c(n_group, k_S)
dim(dem_P) <- c(n_group, k_P)
dim(dem_F) <- c(n_group, k_F)
dim(dem_R) <- c(n_group, k_R)
dim(dem_N) <- n_group

dim(gas_U) <- n_group
dim(gas_A) <- c(n_group, k_A)
dim(gas_E) <- c(n_group, k_E)
dim(gas_S) <- c(n_group, k_S)
dim(gas_P) <- c(n_group, k_P)
dim(gas_F) <- c(n_group, k_F)
dim(gas_R) <- c(n_group, k_R)

dim(n_UE) <- n_group
dim(n_UA) <- n_group
dim(n_AR) <- n_group
dim(n_AU) <- n_group
dim(n_EP) <- n_group
dim(n_ES) <- n_group
dim(n_SR) <- n_group
dim(n_PF) <- n_group
dim(n_FR) <- n_group
dim(n_RU) <- n_group

dim(w) <- n_group
dim(gas_pharyngitis_inc_by_group) <- n_group
dim(new_gas_pharyngitis_inc_by_group) <- n_group
dim(scarlet_fever_inc_by_group) <- n_group
dim(new_scarlet_fever_inc_by_group) <- n_group
dim(weighted_phi_S)   <- n_group

## Comparison to data
##
## Data are weekly. Any missing observation is skipped, so that the same
## comparison serves both the single-group (hydrogen) and age-structured
## (helium) models: the former supplies the all-age scarlet fever rate, the
## latter supplies rates by UKHSA age group.
exp_noise <- parameter(1e6) # rate of exponential noise added to case counts
k_gp <- parameter() # standard deviation of GP surveillance rates
k_hpr <- parameter() # overdispersion of notified cases

## Notified scarlet fever and iGAS cases
obs_scarlet_fever_cases <- data()
obs_igas_inc <- data()
noise_scarlet_fever_cases <- Exponential(rate = exp_noise)
noise_igas_inc <- Exponential(rate = exp_noise)
obs_scarlet_fever_cases ~ NegativeBinomial(
  size = 1 / k_hpr, mu = scarlet_fever_cases + noise_scarlet_fever_cases)
obs_igas_inc ~ NegativeBinomial(
  size = 1 / k_hpr, mu = igas_inc + noise_igas_inc)

## GP consultation rates for scarlet fever
obs_daily_scarlet_fever_rate <- data()
obs_daily_scarlet_fever_rate ~ Normal(daily_scarlet_fever_rate, k_gp)
obs_daily_scarlet_fever_rate_04 <- data()
obs_daily_scarlet_fever_rate_05_14 <- data()
obs_daily_scarlet_fever_rate_15_44 <- data()
obs_daily_scarlet_fever_rate_45_64 <- data()
obs_daily_scarlet_fever_rate_65_74 <- data()
obs_daily_scarlet_fever_rate_75 <- data()
obs_daily_scarlet_fever_rate_04 ~
  Normal(daily_scarlet_fever_rate_04, k_gp)
obs_daily_scarlet_fever_rate_05_14 ~
  Normal(daily_scarlet_fever_rate_05_14, k_gp)
obs_daily_scarlet_fever_rate_15_44 ~
  Normal(daily_scarlet_fever_rate_15_44, k_gp)
obs_daily_scarlet_fever_rate_45_64 ~
  Normal(daily_scarlet_fever_rate_45_64, k_gp)
obs_daily_scarlet_fever_rate_65_74 ~
  Normal(daily_scarlet_fever_rate_65_74, k_gp)
obs_daily_scarlet_fever_rate_75 ~
  Normal(daily_scarlet_fever_rate_75, k_gp)

## GP consultation rates for pharyngitis (all causes). The observed rate is
## scaled by the etiologic fraction and compared with the modelled rate of GAS
## pharyngitis, i.e. obs * etiologic_fraction ~ Normal(modelled, k_gp); this is
## written in terms of the observation as odin requires data on the left.
obs_daily_pharyngitis_rate <- data()
obs_daily_pharyngitis_rate ~ Normal(
  obs_daily_pharyngitis_rate + daily_gas_pharyngitis_rate -
    obs_daily_pharyngitis_rate * etiologic_fraction, k_gp)
obs_daily_pharyngitis_rate_04 <- data()
obs_daily_pharyngitis_rate_05_14 <- data()
obs_daily_pharyngitis_rate_15_44 <- data()
obs_daily_pharyngitis_rate_45_64 <- data()
obs_daily_pharyngitis_rate_65_74 <- data()
obs_daily_pharyngitis_rate_75 <- data()
obs_daily_pharyngitis_rate_04 ~ Normal(
  obs_daily_pharyngitis_rate_04 + daily_gas_pharyngitis_rate_04 -
    obs_daily_pharyngitis_rate_04 * etiologic_fraction_04, k_gp)
obs_daily_pharyngitis_rate_05_14 ~ Normal(
  obs_daily_pharyngitis_rate_05_14 + daily_gas_pharyngitis_rate_05_14 -
    obs_daily_pharyngitis_rate_05_14 * etiologic_fraction_05_14, k_gp)
obs_daily_pharyngitis_rate_15_44 ~ Normal(
  obs_daily_pharyngitis_rate_15_44 + daily_gas_pharyngitis_rate_15_44 -
    obs_daily_pharyngitis_rate_15_44 * etiologic_fraction_15_44, k_gp)
obs_daily_pharyngitis_rate_45_64 ~ Normal(
  obs_daily_pharyngitis_rate_45_64 + daily_gas_pharyngitis_rate_45_64 -
    obs_daily_pharyngitis_rate_45_64 * etiologic_fraction_45_64, k_gp)
obs_daily_pharyngitis_rate_65_74 ~ Normal(
  obs_daily_pharyngitis_rate_65_74 + daily_gas_pharyngitis_rate_65_74 -
    obs_daily_pharyngitis_rate_65_74 * etiologic_fraction_65_74, k_gp)
obs_daily_pharyngitis_rate_75 ~ Normal(
  obs_daily_pharyngitis_rate_75 + daily_gas_pharyngitis_rate_75 -
    obs_daily_pharyngitis_rate_75 * etiologic_fraction_75, k_gp)
