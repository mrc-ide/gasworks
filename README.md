# gasworks

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/gasworks/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/gasworks/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/gasworks/badge)](https://www.codefactor.io/repository/github/mrc-ide/gasworks)
[![codecov](https://codecov.io/gh/mrc-ide/gasworks/branch/main/graph/badge.svg?token=XpGjm3l0ru)](https://codecov.io/gh/mrc-ide/gasworks)
<!-- badges: end -->

This package implements models of Group A Streptococcal (GAS) transmission.


## Installation

You will need a compiler to install dependencies for the package, and to build
the models. Use `pkgbuild::check_build_tools()` to see if your system is usable.

The package is built on [`odin2`](https://mrc-ide.github.io/odin2),
[`dust2`](https://mrc-ide.github.io/dust2) and
[`monty`](https://mrc-ide.github.io/monty), which are available from the
mrc-ide r-universe:

```r
install.packages(
  c("odin2", "dust2", "monty"),
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

The package can then be installed directly from GitHub with:

```r
remotes::install_github("mrc-ide/gasworks", upgrade = FALSE)
```

## Usage

The model is a `dust2` system, with time measured in days:

```r
library(gasworks)
pars <- example_parameters()
y <- model_simulate(pars, times = seq(0, 364, by = 7), n_particles = 10)
matplot(t(y["igas_inc", , ]), type = "l", lty = 1, col = "#00000055",
        xlab = "Week", ylab = "iGAS cases")
```

To fit the single-group (hydrogen) model to weekly surveillance data, create a
particle filter and a packer for the fitted parameters, and combine these with
a prior using `monty`:

```r
filter <- hydrogen_filter(data, n_particles = 100)
packer <- hydrogen_packer(c("beta", "sigma"), demographic_pars, fixed = fixed)
likelihood <- dust2::dust_likelihood_monty(
  filter, packer, save_trajectories = hydrogen_saved_states())
```

The age-structured (helium) model works in the same way, using
`helium_filter()`, `helium_packer()` and `helium_constant_likelihood()`.

## Development

The model is written in odin, in `inst/odin/model.R`. After changing it,
regenerate the compiled code by running `./scripts/generate_odin` from the
package root.

## License

MIT © Imperial College of Science, Technology and Medicine

