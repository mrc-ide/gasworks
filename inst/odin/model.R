## Definition of the time-step and output as "time"
steps_per_week <- 7
dt <- 1 / steps_per_week
initial(time) <- 0
update(time) <- (step + 1) * dt
# i: number of age groups
n_group <- user(1)

## What we really want is min(step + 1, length(alpha_t)) but that's not
## supported by odin (it could be made to support this).
alpha_t <- if (as.integer(step) >= length(alpha))
  alpha[length(alpha)] else alpha[step + 1]

## Core equations for transitions between compartments:
update(U[])  <- U[i]  + dem_U[i]  - n_UE[i] - n_UA[i] + n_AU[i] + n_RU[i]
update(A[])  <- A[i]  + dem_A[i]  + n_UA[i] - n_AU[i] - n_AR[i]
update(E[])  <- E[i]  + dem_E[i]  + n_UE[i] - n_ES[i] - n_EP[i]
update(S1[]) <- S1[i] + dem_S1[i] + n_ES[i] - n_SS[i]
update(S2[]) <- S2[i] + dem_S2[i] + n_SS[i] - n_SR[i]
update(P[])  <- P[i]  + dem_P[i]  + n_EP[i] - n_PF[i]
update(F[])  <- F[i]  + dem_F[i]  + n_PF[i] - n_FR[i]
update(R[])  <- R[i]  + dem_R[i]  + n_AR[i] + n_SR[i] + n_FR[i] - n_RU[i]
update(N[])  <- N[i]  + dem_N[i]

n_Nx[] <- n_Ux[i] + n_Ex[i] + n_Ax[i] + n_S1x[i] + n_S2x[i] + n_Px[i] +
  n_Fx[i] + n_Rx[i]

dem_N[] <- dem_U[i] + dem_E[i] + dem_A[i] + dem_S1[i] + dem_S2[i] + dem_P[i] +
  dem_F[i] + dem_R[i]

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

pharyngitis_inc_by_group[] <- (
  if (step %% steps_per_week == 0) n_SS[i]
  else pharyngitis_inc_by_group[i] + n_SS[i]
  )
update(pharyngitis_inc) <- sum(pharyngitis_inc_by_group[])

scarlet_fever_inc_by_group[] <- (
  if (step %% steps_per_week == 0) n_PF[i]
  else scarlet_fever_inc_by_group[i] + n_PF[i])
update(scarlet_fever_inc) <- sum(scarlet_fever_inc_by_group[])

update(igas_inc) <- (
  if (step %% steps_per_week == 0) sum(n_I[])
  else igas_inc + sum(n_I[])
  )

## Output incidence rates per 100,000 population
## 1. pharyngitis or scarlet fever (prob of reporting / etiologic fraction)
## 2. scarlet fever
pharyngitis_scarlet_fever[] <- pharyngitis_inc_by_group[i] * p_T / phi_S[i] +
  scarlet_fever_inc_by_group[i]

update(pharyngitis_scarlet_fever_rate) <-
  sum(pharyngitis_scarlet_fever[]) / sum(N[]) * 1e5
update(scarlet_fever_rate) <- sum(scarlet_fever_inc_by_group[]) / sum(N[]) * 1e5

## Force of infection
pi <- 3.14159265358979
update(beta_t) <- beta * (1 + sigma * cos(2 * pi * (t0 + step - t_s) / 365.25))
lambda[, ] <- beta_t * m[i, j] *  (A[j] * theta_A + S1[j] + S2[j] + P[j]) / N[j]
foi[] <- sum(lambda[i, ])

## Total rates of transmission out of each compartment
r_U[]  <- foi[i]
r_A[]  <- 1 / delta_A
r_E[]  <- 1 / delta_E
r_S1[] <- 1 / delta_S
r_S2[] <- 1 / delta_S
r_P[]  <- 1 / delta_P
r_F[]  <- 1 / delta_F
r_R[]  <- 1 / delta_R
# rate of developing iGAS
r_I[]  <- foi[i] * p_I

## Calculate number of births
n_xU[1] <- round(alpha_t * dt)

## Calculate net number of leavers from each compartment - deterministic
## note that this can be negative when there is a net population increase due
## to immigration > emigration + death. entrants are distributed across disease
## states in proportion to the general population.
n_Ux[]  <- round(U[i]  * omega[i] * dt)
n_Ax[]  <- round(A[i]  * omega[i] * dt)
n_Ex[]  <- round(E[i]  * omega[i] * dt)
n_S1x[] <- round(S1[i] * omega[i] * dt)
n_S2x[] <- round(S2[i] * omega[i] * dt)
n_Px[]  <- round(P[i]  * omega[i] * dt)
n_Fx[]  <- round(F[i]  * omega[i] * dt)
n_Rx[]  <- round(R[i]  * omega[i] * dt)

## calculate net aging
n_Ui[]  <- (if (i > 1) U[i - 1] else 0) - (if (i < n_group) U[i] else 0)
n_Ai[]  <- (if (i > 1) A[i - 1] else 0) - (if (i < n_group) A[i] else 0)
n_Ei[]  <- (if (i > 1) E[i - 1] else 0) - (if (i < n_group) E[i] else 0)
n_S1i[] <- (if (i > 1) S1[i - 1] else 0) - (if (i < n_group) S1[i] else 0)
n_S2i[] <- (if (i > 1) S2[i - 1] else 0) - (if (i < n_group) S2[i] else 0)
n_Pi[]  <- (if (i > 1) P[i - 1] else 0) - (if (i < n_group) P[i] else 0)
n_Fi[]  <- (if (i > 1) F[i - 1] else 0) - (if (i < n_group) F[i] else 0)
n_Ri[]  <- (if (i > 1) R[i - 1] else 0) - (if (i < n_group) R[i] else 0)

## Calculate all demographic changes
dem_U[]  <- n_xU[i] + round(n_Ui[i] * r_age * dt) - n_Ux[i]
dem_A[]  <- round(n_Ai[i] * r_age * dt) - n_Ax[i]
dem_E[]  <- round(n_Ei[i] * r_age * dt) - n_Ex[i]
dem_S1[] <- round(n_S1i[i] * r_age * dt) - n_S1x[i]
dem_S2[] <- round(n_S2i[i] * r_age * dt) - n_S2x[i]
dem_P[]  <- round(n_Pi[i] * r_age * dt) - n_Px[i]
dem_F[]  <- round(n_Fi[i] * r_age * dt) - n_Fx[i]
dem_R[]  <- round(n_Ri[i] * r_age * dt) - n_Rx[i]

## Draws from binomial distributions for numbers leaving each compartment
## all demographic transitions are done first
n_U[]  <- rbinom(U[i] + dem_U[i],  1 - exp(-r_U[i]  * dt))
n_A[]  <- rbinom(A[i] + dem_A[i],  1 - exp(-r_A[i]  * dt))
n_E[]  <- rbinom(E[i] + dem_E[i],  1 - exp(-r_E[i]  * dt))
n_S1[] <- rbinom(S1[i] + dem_S1[i], 1 - exp(-r_S1[i] * dt))
n_S2[] <- rbinom(S2[i] + dem_S2[i], 1 - exp(-r_S2[i] * dt))
n_P[]  <- rbinom(P[i] + dem_P[i],  1 - exp(-r_P[i]  * dt))
n_F[]  <- rbinom(F[i] + dem_F[i],  1 - exp(-r_F[i]  * dt))
n_R[]  <- rbinom(R[i] + dem_R[i],  1 - exp(-r_R[i]  * dt))

# Number developing iGAS, applies to whole population
n_I[]  <- rbinom(N[i] + dem_N[i],  1 - exp(-r_I[i]  * dt))


## Draw the numbers of transitions between compartments
n_UE[] <- rbinom(n_U[i], p_S)
n_UA[] <- n_U[i] - n_UE[i]
n_AR[] <- rbinom(n_A[i], p_R)
n_AU[] <- n_A[i] - n_AR[i]
n_EP[] <- rbinom(n_E[i], p_F)
n_ES[] <- n_E[i] - n_EP[i]
n_SS[] <- n_S1[i]
n_SR[] <- n_S2[i]
n_PF[] <- n_P[i]
n_FR[] <- n_F[i]
n_RU[] <- n_R[i]

## Initial states:
initial(U[])  <- U0[i]
initial(A[])  <- A0[i]
initial(E[])  <- E0[i]
initial(S1[]) <- S10[i]
initial(S2[]) <- S20[i]
initial(P[])  <- P0[i]
initial(F[])  <- F0[i]
initial(R[])  <- R0[i]
initial(N[])  <- U0[i] + A0[i] + E0[i] + S10[i] + S20[i] + F0[i] + P0[i] + R0[i]
initial(infections_inc) <- 0
initial(pharyngitis_inc) <- 0
initial(scarlet_fever_inc) <- 0
initial(igas_inc) <- 0
initial(births_inc) <- 0
initial(net_leavers_inc) <- 0
initial(beta_t) <- 0
initial(pharyngitis_scarlet_fever_rate) <- 0
initial(scarlet_fever_rate) <- 0

## User defined parameters - default in parentheses:
## Initial number in each state
U0[] <- user()
A0[] <- user()
E0[] <- user()
S10[] <- user()
S20[] <- user()
P0[] <- user()
F0[] <- user()
R0[] <- user()

beta <- user() # rate of transmission
m[, ] <- user()
sigma <- user() # amplitude of seasonal effect
t0 <- user(0)  # day of model initialisation
t_s <- user() # day of peak seasonal transmission
p_S <- user() # probability of pharyngitis symptoms after infection
p_R <- user() # probability of immunity after carriage
p_I <- user() # probability of invasive disease after infection
p_F <- user() # probability of scarlet fever after pharyngitis
p_T <- user() # probability of seeking treatment for pharyngitis
delta_A <- user() # mean duration of carriage
delta_E <- user() # mean duration of incubation period
delta_S <- user() # mean duration of pharyngitis symptoms (x 2)
delta_P <- user() # mean duration from pharyngitis to scarlet fever rash
delta_F <- user() # mean duration of scarlet fever rash
delta_R <- user() # mean duration of natural immunity
theta_A <- user() # infectiousness of carriers relative to symptomatics
phi_S[] <- user() # proportion of all pharyngitis attributable to GAS

alpha[]  <- user() # time-varying number of births
dim(alpha) <- user()
omega[] <- user() # rate of population exit, can be negative
r_age <- user(0)   # rate of aging - determined by group size

## Object dimensions

dim(m)      <- c(n_group, n_group)
dim(lambda) <- c(n_group, n_group)
dim(foi)   <- n_group
dim(omega) <- n_group
dim(phi_S) <- n_group

dim(U)  <- n_group
dim(A)  <- n_group
dim(E)  <- n_group
dim(S1) <- n_group
dim(S2) <- n_group
dim(P)  <- n_group
dim(F)  <- n_group
dim(R)  <- n_group
dim(N)  <- n_group

dim(U0)  <- n_group
dim(A0)  <- n_group
dim(E0)  <- n_group
dim(S10) <- n_group
dim(S20) <- n_group
dim(P0)  <- n_group
dim(F0)  <- n_group
dim(R0)  <- n_group


dim(r_U) <- n_group
dim(r_A) <- n_group
dim(r_E) <- n_group
dim(r_S1) <- n_group
dim(r_S2) <- n_group
dim(r_P) <- n_group
dim(r_F) <- n_group
dim(r_R) <- n_group
dim(r_I) <- n_group

dim(n_xU) <- n_group

dim(n_U) <- n_group
dim(n_A) <- n_group
dim(n_E) <- n_group
dim(n_S1) <- n_group
dim(n_S2) <- n_group
dim(n_P) <- n_group
dim(n_F) <- n_group
dim(n_R) <- n_group
dim(n_I) <- n_group

dim(n_Ux) <- n_group
dim(n_Ax) <- n_group
dim(n_Ex) <- n_group
dim(n_S1x) <- n_group
dim(n_S2x) <- n_group
dim(n_Px) <- n_group
dim(n_Fx) <- n_group
dim(n_Rx) <- n_group
dim(n_Nx) <- n_group

dim(n_Ui) <- n_group
dim(n_Ai) <- n_group
dim(n_Ei) <- n_group
dim(n_S1i) <- n_group
dim(n_S2i) <- n_group
dim(n_Pi) <- n_group
dim(n_Fi) <- n_group
dim(n_Ri) <- n_group

dim(dem_U) <- n_group
dim(dem_A) <- n_group
dim(dem_E) <- n_group
dim(dem_S1) <- n_group
dim(dem_S2) <- n_group
dim(dem_P) <- n_group
dim(dem_F) <- n_group
dim(dem_R) <- n_group
dim(dem_N) <- n_group

dim(n_UE) <- n_group
dim(n_UA) <- n_group
dim(n_AR) <- n_group
dim(n_AU) <- n_group
dim(n_EP) <- n_group
dim(n_ES) <- n_group
dim(n_SS) <- n_group
dim(n_SR) <- n_group
dim(n_PF) <- n_group
dim(n_FR) <- n_group
dim(n_RU) <- n_group

dim(pharyngitis_inc_by_group) <- n_group
dim(scarlet_fever_inc_by_group) <- n_group
dim(pharyngitis_scarlet_fever) <- n_group
