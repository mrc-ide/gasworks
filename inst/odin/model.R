## Definition of the time-step and output as "time"
steps_per_week <- 7
dt <- 1 / steps_per_week
initial(time) <- 0
update(time) <- (step + 1) * dt
# i: number of age groups
n_group <- user(1)

## Core equations for transitions between compartments:
update(U[])  <- U[i]  + n_xU[i] - n_UE[i] - n_UA[i] + n_AU[i] + n_RU[i] -
  n_Ux[i]
update(E[])  <- E[i]  + n_UE[i] - n_ES[i] - n_EI[i] - n_Ex[i]
update(A[])  <- A[i]  + n_UA[i] - n_AU[i] - n_AR[i] - n_Ax[i]
update(S1[]) <- S1[i] + n_ES[i] - n_SF[i] - n_SS[i] - n_S1x[i]
update(S2[]) <- S2[i] + n_SS[i] - n_SR[i] - n_S2x[i]
update(F[])  <- F[i]  + n_SF[i] - n_FR[i] - n_Fx[i]
update(I[])  <- I[i]  + n_EI[i] - n_IR[i] - n_Ix[i]
update(R[])  <- R[i]  + n_AR[i] + n_SR[i] + n_FR[i] + n_IR[i] - n_RU[i] -
  n_Rx[i]
update(N[])  <- N[i]  + n_xU[i] - n_Nx[i]

n_Nx[] <- n_Ux[i] + n_Ex[i] + n_Ax[i] + n_S1x[i] + n_S2x[i] + n_Fx[i] +
  n_Ix[i] + n_Rx[i]

## Output incidence flows:
update(entrants_inc) <- (
  if (step %% steps_per_week == 0) sum(n_xU[])
  else entrants_inc + sum(n_xU[]))

update(leavers_inc) <- (
  if (step %% steps_per_week == 0) sum(n_Nx[])
  else leavers_inc + sum(n_Nx[]))

update(infections_inc) <- (
  if (step %% steps_per_week == 0) sum(n_UE[]) + sum(n_UA[])
  else infections_inc + sum(n_UE[]) + sum(n_UA[]))

pharyngitis_inc_by_group[] <- (
  if (step %% steps_per_week == 0) n_SS[i] + n_SF[i]
  else pharyngitis_inc_by_group[i] + n_SS[i] + n_SF[i])
update(pharyngitis_inc) <- sum(pharyngitis_inc_by_group[])

scarlet_fever_inc_by_group[] <- (
  if (step %% steps_per_week == 0) sum(n_SF[i])
  else scarlet_fever_inc_by_group[i] + n_SF[i])
update(scarlet_fever_inc) <- sum(scarlet_fever_inc_by_group[])

update(igas_inc) <- (
  if (step %% steps_per_week == 0) sum(n_EI[])
  else igas_inc + sum(n_EI[]))

## Output incidence rates per 100,000 population
all_pharyngitis[] <- pharyngitis_inc_by_group[i] / phi_S[i]
update(pharyngitis_rate) <- sum(all_pharyngitis[]) / sum(N[]) * 1e5
update(scarlet_fever_rate) <- sum(scarlet_fever_inc_by_group[]) / sum(N[]) * 1e5

## Force of infection
pi <- 3.14159265358979
update(beta_t) <- beta * (1 + sigma * cos(2 * pi * (t0 + step - t_s) / 365.25))
lambda[, ] <- beta_t * m[i, j] *  (A[j] * theta_A + S1[j] + S2[j]) / N[j]
foi[] <- sum(lambda[i, ])


## Rates of transition between compartments
r_UE[] <- p_S * foi[i]
r_UA[] <- (1 - p_S) * foi[i]
r_AR[] <- p_R / delta_A
r_AU[] <- (1 - p_R) / delta_A
r_EI[] <- p_I / delta_E
r_ES[] <- (1 - p_I) / delta_E
r_IR[] <- 1 / delta_I
r_SF[] <- p_F / delta_S
r_SS[] <- (1 - p_F) / delta_S
r_SR[] <- 1 / delta_S
r_FR[] <- 1 / delta_F
r_RU[] <- 1 / delta_R

## Total rates of transmission out of each compartment
r_U[]  <- r_UE[i] + r_UA[i] + omega[i]
r_A[]  <- r_AU[i] + r_AR[i] + omega[i]
r_E[]  <- r_EI[i] + r_ES[i] + omega[i]
r_I[]  <- r_IR[i] + omega[i]
r_S1[] <- r_SF[i] + r_SS[i] + omega[i]
r_S2[] <- r_SR[i] + omega[i]
r_F[]  <- r_FR[i] + omega[i]
r_R[]  <- r_RU[i] + omega[i]

##Draw number of entrants
n_xU[] <- rpois(alpha[i] * dt)

## Draws from binomial distributions for numbers leaving each compartments
n_U[]  <- rbinom(U[i],  1 - exp(-r_U[i]  * dt))
n_A[]  <- rbinom(A[i],  1 - exp(-r_A[i]  * dt))
n_E[]  <- rbinom(E[i],  1 - exp(-r_E[i]  * dt))
n_I[]  <- rbinom(I[i],  1 - exp(-r_I[i]  * dt))
n_S1[] <- rbinom(S1[i], 1 - exp(-r_S1[i] * dt))
n_S2[] <- rbinom(S2[i], 1 - exp(-r_S2[i] * dt))
n_F[]  <- rbinom(F[i],  1 - exp(-r_F[i]  * dt))
n_R[]  <- rbinom(R[i],  1 - exp(-r_R[i]  * dt))

# Draw the number of leavers from each compartment
n_Ux[]  <- rbinom(n_U[i],  omega[i] / r_U[i])
n_Ax[]  <- rbinom(n_A[i],  omega[i] / r_A[i])
n_Ex[]  <- rbinom(n_E[i],  omega[i] / r_E[i])
n_Ix[]  <- rbinom(n_I[i],  omega[i] / r_I[i])
n_S1x[] <- rbinom(n_S1[i], omega[i] / r_S1[i])
n_S2x[] <- rbinom(n_S2[i], omega[i] / r_S2[i])
n_Fx[]  <- rbinom(n_F[i],  omega[i] / r_F[i])
n_Rx[]  <- rbinom(n_R[i],  omega[i] / r_R[i])

## Draw the numbers of transitions between compartments
n_UE[] <- rbinom(n_U[i] - n_Ux[i], p_S)
n_UA[] <- n_U[i] - n_Ux[i] - n_UE[i]
n_AR[] <- rbinom(n_A[i] - n_Ax[i], p_R)
n_AU[] <- n_A[i] - n_Ax[i] - n_AR[i]
n_EI[] <- rbinom(n_E[i] - n_Ex[i], p_I)
n_ES[] <- n_E[i] - n_Ex[i] - n_EI[i]
n_IR[] <- n_I[i] - n_Ix[i]
n_SF[] <- rbinom(n_S1[i] - n_S1x[i], p_F)
n_SS[] <- n_S1[i] - n_S1x[i] - n_SF[i]
n_SR[] <- n_S2[i] - n_S2x[i]
n_FR[] <- n_F[i] - n_Fx[i]
n_RU[] <- n_R[i] - n_Rx[i]

## Initial states:
initial(U[])  <- U0[i]
initial(A[])  <- A0[i]
initial(E[])  <- E0[i]
initial(I[])  <- I0[i]
initial(S1[]) <- S10[i]
initial(S2[]) <- S20[i]
initial(F[])  <- F0[i]
initial(R[])  <- R0[i]
initial(N[])  <- U0[i] + A0[i] + E0[i] + I0[i] + S10[i] + S20[i] + F0[i] + R0[i]
initial(infections_inc) <- 0
initial(pharyngitis_inc) <- 0
initial(scarlet_fever_inc) <- 0
initial(igas_inc) <- 0
initial(entrants_inc) <- 0
initial(leavers_inc) <- 0
initial(beta_t) <- 0
initial(pharyngitis_rate) <- 0
initial(scarlet_fever_rate) <- 0

## User defined parameters - default in parentheses:
## Initial number in each state
U0[] <- user()
A0[] <- user()
E0[] <- user()
I0[] <- user()
S10[] <- user()
S20[] <- user()
F0[] <- user()
R0[] <- user()

beta <- user() # rate of transmission
m[, ] <- user()
sigma <- user() # amplitude of seasonal effect
t0 <- user(0)  # day of model initialisation
t_s <- user() # day of peak seasonal transmission
p_S <- user() # probability of any symptoms after infection
p_R <- user() # probability of immunity after carriage
p_I <- user() # probability of invasive disease after infection
p_F <- user() # probability of scarlet fever after pharyngitis
delta_A <- user() # mean duration of carriage
delta_E <- user() # mean duration of incubation period
delta_I <- user() # mean duration of invasive disease
delta_S <- user() # mean duration of pharyngitis symptoms (x 2)
delta_F <- user() # mean duration of scarlet fever
delta_R <- user() # mean duration of natural immunity
theta_A <- user() # infectiousness of carriers relative to symptomatics
phi_S[] <- user() # proportion of all pharyngitis attributable to GAS

alpha[] <- user() # number of population entrants
omega[] <- user() # rate of population exit

## Object dimensions

dim(m)      <- c(n_group, n_group)
dim(lambda) <- c(n_group, n_group)
dim(foi)   <- n_group
dim(alpha) <- n_group
dim(omega) <- n_group
dim(phi_S) <- n_group

dim(U)  <- n_group
dim(A)  <- n_group
dim(E)  <- n_group
dim(I)  <- n_group
dim(S1) <- n_group
dim(S2) <- n_group
dim(F)  <- n_group
dim(R)  <- n_group
dim(N)  <- n_group

dim(U0)  <- n_group
dim(A0)  <- n_group
dim(E0)  <- n_group
dim(I0)  <- n_group
dim(S10) <- n_group
dim(S20) <- n_group
dim(F0)  <- n_group
dim(R0)  <- n_group

dim(r_UE) <- n_group
dim(r_UA) <- n_group
dim(r_AR) <- n_group
dim(r_AU) <- n_group
dim(r_EI) <- n_group
dim(r_ES) <- n_group
dim(r_IR) <- n_group
dim(r_SF) <- n_group
dim(r_SS) <- n_group
dim(r_SR) <- n_group
dim(r_FR) <- n_group
dim(r_RU) <- n_group

dim(r_U) <- n_group
dim(r_A) <- n_group
dim(r_E) <- n_group
dim(r_I) <- n_group
dim(r_S1) <- n_group
dim(r_S2) <- n_group
dim(r_F) <- n_group
dim(r_R) <- n_group

dim(n_xU) <- n_group

dim(n_U) <- n_group
dim(n_A) <- n_group
dim(n_E) <- n_group
dim(n_I) <- n_group
dim(n_S1) <- n_group
dim(n_S2) <- n_group
dim(n_F) <- n_group
dim(n_R) <- n_group

dim(n_Ux) <- n_group
dim(n_Ax) <- n_group
dim(n_Ex) <- n_group
dim(n_Ix) <- n_group
dim(n_S1x) <- n_group
dim(n_S2x) <- n_group
dim(n_Fx) <- n_group
dim(n_Rx) <- n_group
dim(n_Nx) <- n_group

dim(n_UE) <- n_group
dim(n_UA) <- n_group
dim(n_AR) <- n_group
dim(n_AU) <- n_group
dim(n_EI) <- n_group
dim(n_ES) <- n_group
dim(n_IR) <- n_group
dim(n_SF) <- n_group
dim(n_SS) <- n_group
dim(n_SR) <- n_group
dim(n_FR) <- n_group
dim(n_RU) <- n_group

dim(pharyngitis_inc_by_group) <- n_group
dim(scarlet_fever_inc_by_group) <- n_group
dim(all_pharyngitis) <- n_group
