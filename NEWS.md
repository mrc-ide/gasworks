# gasworks 0.1.9

* update model date functions to align with UKHSA epi weeks, which end on Sunday
i.e. model week 1 ends on Sun 5 Jan 2014.

# gasworks 0.1.8

* output daily_pharyngitis_scarlet_fever_rate and daily_scarlet_fever_rate
by age groups used in UKHSA data
* ensure incidence rates are daily 7-day moving averages, and clarify in naming

# gasworks 0.1.7

* allow time in each compartment to be Erlang distributed

# gasworks 0.1.6

* add pre scarlet fever compartment P (pharyngitis symptoms but no rash)

# gasworks 0.1.5

* remove iGAS compartment - output iGAS cases at observation level

# gasworks 0.1.4

* add time-varying births
* add under-reporting of pharyngitis (p_T)

# gasworks 0.1.3

* allow immigration via negative omega

# gasworks 0.1.2

* add aging between groups
* make all demographic processes deterministic

# gasworks 0.1.1

* add age-structure to model

# gasworks 0.1.0

* add stochastic compartmental model of GAS transmission and pmcmc fitting
capabilities via mrc-ide/mcstate.git

# gasworks 0.0.1

* initial setup
