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

You will need the packages `odin` and `mcstate`, which can be installed using:

```r
remotes::install_github("mrc-ide/odin", upgrade = FALSE)
remotes::install_github("mrc-ide/mcstate", upgrade = FALSE)
```


The package can then be installed directly from GitHub with:

```r
remotes::install_github("mrc-ide/gasworks", upgrade = FALSE)
```

## License

MIT © Imperial College of Science, Technology and Medicine

