## Definition of the time-step and output as "time"
steps_per_week <- 7
dt <- 1 / steps_per_week
initial(time) <- 0
update(time) <- (step + 1) * dt
# i: number of age groups, j: number of Erlang compartments
n_group <- user(1)

## What we really want is min(step + 1, length(alpha_t)) but that's not
## supported by odin (it could be made to support this).
alpha_t <- if (as.integer(step) >= length(alpha))
  alpha[length(alpha)] else alpha[step + 1]
## Probability of reporting a scarlet fever case (time-varying)
q_F_t <- if (as.integer(step) >= length(q_F))
  q_F[length(q_F)] else q_F[step + 1]

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
update(births_inc) <- (
  if (step %% steps_per_week == 0) sum(n_xU[])
  else births_inc + sum(n_xU[])
  )

update(net_leavers_inc) <- (
  if (step %% steps_per_week == 0) sum(n_Nx[])
  else net_leavers_inc + sum(n_Nx[])
  )

update(infections_inc) <- (
  if (step %% steps_per_week == 0) sum(n_UE[]) + sum(n_UA[])
  else infections_inc + sum(n_UE[]) + sum(n_UA[])
  )

gas_pharyngitis_inc_by_group[] <- (
  if (step %% steps_per_week == 0) n_ES[i]
  else gas_pharyngitis_inc_by_group[i] + n_ES[i]
  )

scarlet_fever_inc_by_group[] <- (
  if (step %% steps_per_week == 0) n_PF[i]
  else scarlet_fever_inc_by_group[i] + n_PF[i])

# output all-group incidence flows
update(gas_pharyngitis_inc) <- sum(gas_pharyngitis_inc_by_group[])
update(scarlet_fever_inc) <- sum(scarlet_fever_inc_by_group[])
update(scarlet_fever_cases) <- round(sum(scarlet_fever_inc_by_group[]) * q_F_t)
update(igas_inc) <- (
  if (step %% steps_per_week == 0) sum(n_I[]) else igas_inc + sum(n_I[]))

## Output daily incidence rate per 100,000 population
 w[] <- N[i] / 1e5 * 7
update(daily_gas_pharyngitis_rate) <-
  sum(gas_pharyngitis_inc_by_group[]) / sum(w[])

## When using the age-structured model (helium) output incidence rates for age
## groups used by UK GP surveillance:
##  04,  5-14, 15-44,   45-64,   65-74, 75+
## [1], [2:3], [4:9], [10:13], [14:15], [16]

## Output daily GAS pharyngitis rates by group
update(daily_gas_pharyngitis_rate_04) <- (
  if (n_group == 16) gas_pharyngitis_inc_by_group[1] / w[1]
  else 0)
update(daily_gas_pharyngitis_rate_05_14) <- (
  if (n_group == 16) sum(gas_pharyngitis_inc_by_group[2:3]) / sum(w[2:3])
  else 0)
update(daily_gas_pharyngitis_rate_15_44) <- (
  if (n_group == 16) sum(gas_pharyngitis_inc_by_group[4:9]) / sum(w[4:9])
  else 0)
update(daily_gas_pharyngitis_rate_45_64) <- (
  if (n_group == 16) sum(gas_pharyngitis_inc_by_group[10:13]) / sum(w[10:13])
  else 0)
update(daily_gas_pharyngitis_rate_65_74) <- (
  if (n_group == 16) sum(gas_pharyngitis_inc_by_group[14:15]) / sum(w[14:15])
  else 0)
update(daily_gas_pharyngitis_rate_75) <- (
  if (n_group == 16) gas_pharyngitis_inc_by_group[16] / w[16]
  else 0)

## Output daily scarlet fever rates by group
update(daily_scarlet_fever_rate_04) <- (
  if (n_group == 16) scarlet_fever_inc_by_group[1] / w[1]
  else 0)
update(daily_scarlet_fever_rate_05_14) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[2:3]) / sum(w[2:3])
  else 0)
update(daily_scarlet_fever_rate_15_44) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[4:9]) / sum(w[4:9])
  else 0)
update(daily_scarlet_fever_rate_45_64) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[10:13]) / sum(w[10:13])
  else 0)
update(daily_scarlet_fever_rate_65_74) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[14:15]) / sum(w[14:15])
  else 0)
update(daily_scarlet_fever_rate_75) <- (
  if (n_group == 16) scarlet_fever_inc_by_group[16] / w[16]
  else 0)

## Output weekly scarlet fever cases by group
update(scarlet_fever_inc_04) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[1]) else 0)
update(scarlet_fever_inc_05_14) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[2:3]) else 0)
update(scarlet_fever_inc_15_44) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[4:9]) else 0)
update(scarlet_fever_inc_45_64) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[10:13]) else 0)
update(scarlet_fever_inc_65_74) <- (
  if (n_group == 16) sum(scarlet_fever_inc_by_group[14:15]) else 0)
update(scarlet_fever_inc_75) <- (
  if (n_group == 16) scarlet_fever_inc_by_group[16] else 0)

## Force of infection
pi <- 3.14159265358979
update(beta_t) <- beta * (1 + sigma * cos(2 * pi * (t0 + step - t_s) / 365.25))
lambda[, ] <- beta_t * m[i, j] *
  (sum(A[j, ]) * theta_A + sum(S[j, ]) + sum(P[j, ])) / N[j]
foi[] <- sum(lambda[i, ])

## Total rates of transmission out of each (sub)compartment
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
n_xU[1] <- round(alpha_t * dt)

## Calculate net number of leavers from each compartment - deterministic
## note that this can be negative when there is a net population increase due
## to immigration > emigration + death. entrants are distributed across disease
## states in proportion to the general population.
n_Ux[]   <- round(U[i] * omega[i] * dt)
n_Ax[, ] <- round(A[i, j] * omega[i] * dt)
n_Ex[, ] <- round(E[i, j] * omega[i] * dt)
n_Sx[, ] <- round(S[i, j] * omega[i] * dt)
n_Px[, ] <- round(P[i, j] * omega[i] * dt)
n_Fx[, ] <- round(F[i, j] * omega[i] * dt)
n_Rx[, ] <- round(R[i, j] * omega[i] * dt)

## calculate net aging
n_Ui[] <- (if (i > 1) U[i - 1] else 0) - (if (i < n_group) U[i] else 0)
n_Ai[, ] <- (if (i > 1) A[i - 1, j] else 0) - (if (i < n_group) A[i, j] else 0)
n_Ei[, ] <- (if (i > 1) E[i - 1, j] else 0) - (if (i < n_group) E[i, j] else 0)
n_Si[, ] <- (if (i > 1) S[i - 1, j] else 0) - (if (i < n_group) S[i, j] else 0)
n_Pi[, ] <- (if (i > 1) P[i - 1, j] else 0) - (if (i < n_group) P[i, j] else 0)
n_Fi[, ] <- (if (i > 1) F[i - 1, j] else 0) - (if (i < n_group) F[i, j] else 0)
n_Ri[, ] <- (if (i > 1) R[i - 1, j] else 0) - (if (i < n_group) R[i, j] else 0)

## Calculate all demographic changes
dem_U[]   <- n_xU[i] + round(n_Ui[i] * r_age * dt) - n_Ux[i]
dem_A[, ] <- round(n_Ai[i, j] * r_age * dt) - n_Ax[i, j]
dem_E[, ] <- round(n_Ei[i, j] * r_age * dt) - n_Ex[i, j]
dem_S[, ] <- round(n_Si[i, j] * r_age * dt) - n_Sx[i, j]
dem_P[, ] <- round(n_Pi[i, j] * r_age * dt) - n_Px[i, j]
dem_F[, ] <- round(n_Fi[i, j] * r_age * dt) - n_Fx[i, j]
dem_R[, ] <- round(n_Ri[i, j] * r_age * dt) - n_Rx[i, j]

## Draws from binomial distributions for numbers leaving each compartment
## all demographic transitions are done first
n_U[]   <- rbinom(U[i]    + dem_U[i],    1 - exp(-r_U[i] * dt))
n_A[, ] <- rbinom(A[i, j] + dem_A[i, j], 1 - exp(-r_A[i] * dt))
n_E[, ] <- rbinom(E[i, j] + dem_E[i, j], 1 - exp(-r_E[i] * dt))
n_S[, ] <- rbinom(S[i, j] + dem_S[i, j], 1 - exp(-r_S[i] * dt))
n_P[, ] <- rbinom(P[i, j] + dem_P[i, j], 1 - exp(-r_P[i] * dt))
n_F[, ] <- rbinom(F[i, j] + dem_F[i, j], 1 - exp(-r_F[i] * dt))
n_R[, ] <- rbinom(R[i, j] + dem_R[i, j], 1 - exp(-r_R[i] * dt))

# Number developing iGAS, applies to whole population
n_I[] <- rbinom(N[i] + dem_N[i], 1 - exp(-r_I[i] * dt))

## Draw the numbers of transitions between compartments
n_UE[] <- rbinom(n_U[i], p_S)
n_UA[] <- n_U[i] - n_UE[i]
n_AR[] <- rbinom(n_A[i, k_A], p_R)
n_AU[] <- n_A[i, k_A] - n_AR[i]
n_EP[] <- rbinom(n_E[i, k_E], p_F)
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

initial(infections_inc) <- 0
initial(gas_pharyngitis_inc) <- 0
initial(scarlet_fever_inc) <- 0
initial(scarlet_fever_cases) <- 0
initial(igas_inc) <- 0
initial(births_inc) <- 0
initial(net_leavers_inc) <- 0
initial(beta_t) <- 0
initial(daily_gas_pharyngitis_rate) <- 0

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

## User defined parameters - default in parentheses:
## Initial number in each state
U0[] <- user()
A0[, ] <- user()
E0[, ] <- user()
S0[, ] <- user()
P0[, ] <- user()
F0[, ] <- user()
R0[, ] <- user()

beta <- user() # rate of transmission
m[, ] <- user()
sigma <- user() # amplitude of seasonal effect
t0 <- user(0)  # day of model initialisation
t_s <- user() # day of peak seasonal transmission
p_S <- user() # probability of pharyngitis symptoms after infection
p_R <- user() # probability of immunity after carriage
p_I <- user() # probability of invasive disease after infection
p_F <- user() # probability of scarlet fever after pharyngitis
delta_A <- user() # mean duration of carriage
delta_E <- user() # mean duration of incubation period
delta_S <- user() # mean duration of pharyngitis symptoms
delta_P <- user() # mean duration from pharyngitis to scarlet fever rash
delta_F <- user() # mean duration of scarlet fever rash
delta_R <- user() # mean duration of natural immunity
k_A <- user() # number of sub-compartments for carriage
k_E <- user() # number of sub-compartments for incubation period
k_S <- user() # number of sub-compartments for pharyngitis symptoms
k_P <- user() # number of sub-compartments for pre-scarlet fever
k_F <- user() # number of sub-compartments for scarlet fever rash
k_R <- user() # number of sub-compartments for natural immunity
theta_A <- user() # infectiousness of carriers relative to symptomatics

alpha[]  <- user() # time-varying number of births
dim(alpha) <- user()
omega[] <- user() # rate of population exit, can be negative
r_age <- user(0)   # rate of aging - determined by group size

q_F[] <- user() # time-varying probability of reporting a scarlet fever case
dim(q_F) <- user()

## Object dimensions

dim(m)      <- c(n_group, n_group)
dim(lambda) <- c(n_group, n_group)
dim(foi)   <- n_group
dim(omega) <- n_group

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
dim(scarlet_fever_inc_by_group)   <- n_group
