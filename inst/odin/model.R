## Definition of the time-step and output as "time"
steps_per_week <- 7
dt <- 1 / steps_per_week
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(U)  <- U  + n_xU - n_UE - n_UA + n_AU + n_RU - n_Ux
update(E)  <- E  + n_UE - n_ES - n_EI - n_Ex
update(A)  <- A  + n_UA - n_AU - n_AR - n_Ax
update(S1) <- S1 + n_ES - n_SF - n_SS - n_S1x
update(S2) <- S2 + n_SS - n_SR - n_S2x
update(F)  <- F  + n_SF - n_FR - n_Fx
update(I)  <- I  + n_EI - n_IR - n_Ix
update(R)  <- R  + n_AR + n_SR + n_FR + n_IR - n_RU - n_Rx
update(N)  <- N  + n_xU - n_Nx

n_Nx <- n_Ux + n_Ex + n_Ax + n_S1x + n_S2x + n_Fx + n_Ix + n_Rx


## Output incidence flows:
update(infections_inc) <- (if (step %% steps_per_week == 0) n_UE + n_UA
                           else infections_inc + n_UE + n_UA)
update(pharyngitis_inc) <- (if (step %% steps_per_week == 0) n_SS + n_SF
                            else pharyngitis_inc + n_SS + n_SF)
update(scarlet_fever_inc) <- (if (step %% steps_per_week == 0) n_SF
                              else scarlet_fever_inc + n_SF)
update(igas_inc) <- (if (step %% steps_per_week == 0) n_EI
                     else igas_inc + n_EI)
update(entrants_inc) <- (if (step %% steps_per_week == 0) n_xU
                         else entrants_inc + n_xU)
update(leavers_inc) <- (if (step %% steps_per_week == 0) n_Nx
                        else leavers_inc + n_Nx)

## Force of infection
pi <- 3.14159265358979
seasonality <-  1 + sigma * cos(2 * pi * (t0 + step - t_s) / 365.25)
lambda <- beta * seasonality * (A + S1 + S2) / N
update(foi) <- lambda

## Rates of transition between compartments
r_UE <- p_S * lambda
r_UA <- (1 - p_S) * lambda
r_AR <- p_R / delta_A
r_AU <- (1 - p_R) / delta_A
r_EI <- p_I / delta_E
r_ES <- (1 - p_I) / delta_E
r_IR <- 1 / delta_I
r_SF <- p_F / delta_S
r_SS <- (1 - p_F) / delta_S
r_SR <- 1 / delta_S
r_FR <- 1 / delta_F
r_RU <- 1 / delta_R

## Total rates of transmission out of each compartment
r_U  <- r_UE + r_UA + omega
r_A  <- r_AU + r_AR + omega
r_E  <- r_EI + r_ES + omega
r_I  <- r_IR + omega
r_S1 <- r_SF + r_SS + omega
r_S2 <- r_SR + omega
r_F  <- r_FR + omega
r_R  <- r_RU + omega

##Draw number of entrants
n_xU <- rpois(alpha * dt)

## Draws from binomial distributions for numbers leaving each compartments
n_U  <- rbinom(U,  1 - exp(-r_U  * dt))
n_A  <- rbinom(A,  1 - exp(-r_A  * dt))
n_E  <- rbinom(E,  1 - exp(-r_E  * dt))
n_I  <- rbinom(I,  1 - exp(-r_I  * dt))
n_S1 <- rbinom(S1, 1 - exp(-r_S1 * dt))
n_S2 <- rbinom(S2, 1 - exp(-r_S2 * dt))
n_F  <- rbinom(F,  1 - exp(-r_F  * dt))
n_R  <- rbinom(R,  1 - exp(-r_R  * dt))

# Draw the number of leavers from each compartment
n_Ux  <- rbinom(n_U,  omega / r_U)
n_Ax  <- rbinom(n_A,  omega / r_A)
n_Ex  <- rbinom(n_E,  omega / r_E)
n_Ix  <- rbinom(n_I,  omega / r_I)
n_S1x <- rbinom(n_S1, omega / r_S1)
n_S2x <- rbinom(n_S2, omega / r_S2)
n_Fx  <- rbinom(n_F,  omega / r_F)
n_Rx  <- rbinom(n_R,  omega / r_R)

## Draw the numbers of transitions between compartments
n_UE <- rbinom(n_U - n_Ux, p_S)
n_UA <- n_U - n_Ux - n_UE
n_AR <- rbinom(n_A - n_Ax, p_R)
n_AU <- n_A - n_Ax - n_AR
n_EI <- rbinom(n_E - n_Ex, p_I)
n_ES <- n_E - n_Ex - n_EI
n_IR <- n_I - n_Ix
n_SF <- rbinom(n_S1 - n_S1x, p_F)
n_SS <- n_S1 - n_S1x - n_SF
n_SR <- n_S2 - n_S2x
n_FR <- n_F - n_Fx
n_RU <- n_R - n_Rx

## Initial states:
initial(U)  <- U0
initial(A)  <- A0
initial(E)  <- E0
initial(I)  <- I0
initial(S1) <- S10
initial(S2) <- S20
initial(F)  <- F0
initial(R)  <- R0
initial(N)  <- U0 + A0 + E0 + I0 + S10 + S20 + F0 + R0
initial(infections_inc) <- 0
initial(pharyngitis_inc) <- 0
initial(scarlet_fever_inc) <- 0
initial(igas_inc) <- 0
initial(entrants_inc) <- 0
initial(leavers_inc) <- 0
initial(foi) <- 0

## User defined parameters - default in parentheses:
## Initial number in each state
U0  <- user(0)
A0  <- user(0)
E0  <- user(0)
I0  <- user(0)
S10 <- user(0)
S20 <- user(0)
F0  <- user(0)
R0  <- user(0)

beta <- user() # rate of transmission
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

alpha <- user() # number of population entrants
omega <- user() # rate of population exit
