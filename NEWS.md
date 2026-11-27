# gasworks 0.4.0

* Port the model to `odin2` and `dust2`, and the fitting functions to `monty`, replacing `odin.dust`, `dust` and `mcstate`. In deterministic mode the model reproduces the trajectories and likelihood of version 0.3.6.
* Time is now measured in days (with a default time step of one day) and durations (`delta_*`) are passed to the model in days. Transmission and demographic rates are still given per week.
* The comparison to data is now written in odin, and is always compiled. `hydrogen_compare()` and `helium_compare()` remain as reference implementations in R.
* `hydrogen_filter()` and `helium_filter()` return `dust2` filters. The likelihood of the time-invariant data used by the helium model is now created with `helium_constant_likelihood()`.
* New functions `model_simulate()`, `hydrogen_packer()`, `helium_packer()`, `hydrogen_saved_states()` and `helium_saved_states()`; the last two replace `hydrogen_index()` and `helium_index()`.
* `model_parameters()` no longer replaces supplied values of `k_gp` and `k_hpr` with 1, which meant that the fitted values were not used in the likelihood.
* There is no transmission from an age group with no population (previously this gave a force of infection of `NaN`).

# gasworks 0.3.5

* Compiled compare function allowing running on a GPU

# gasworks 0.3.3

* helium_compare: use total daily pharyngitis rate when age-disaggregated unavaliable
* outputt UKHSA-aggregated etiologic fraction from model

# gasworks 0.3.2

* update helium and hydrogen compare to fit to GP scarlet fever incidence,
as well as notifiable cases (subject to time-varying under-reporting)

# gasworks 0.3.1

* update helium compare to use Normal for pharyngitis rates and Dirichlet for scarlet fever proportions

# gasworks 0.3.0

* add running and fitting functions for both hydrogen and helium models

# gasworks 0.2.0

* reduce helium age groups to 16, so max group is 75+

# gasworks 0.1.12

* output new model states: prev_A and prev_R, will scale by age

# gasworks 0.1.11

* output new model states: gas_pharyngitis_inc, pharyngitis_inc (all-cause
cases)
* when n_group == 19, additionally output: 
proportion of scarlet fever cases, proportion of all-cause pharyngitis cases,
proportion of all-cause pharyngitis or scarlet fever cases in each  UKHSA age
group

# gasworks 0.1.10

* update model date functions to align with UKHSA epi weeks, which end on Sunday
i.e. model week 1 ends on Sun 5 Jan 2014.

# gasworks 0.1.9

* output daily incidence rates (based on 7-day average) for scarlet fever and
pharyngitis-or-scarlet-fever
* allow fitting of daily incidence rates by age-groups reported in GP
surveillance


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
