# Pre-workshop

Pre-workshop instructions for participants. Let’s get ready to rock! 🚀

## R and RStudio 💻

[R](https://en.wikipedia.org/wiki/R_(programming_language)) is a fantastic software for statistical analyses. 📊 [RStudio](https://posit.co/products/open-source/rstudio/) is your trusty sidekick, helping you navigate the R universe with ease. It’s like a cozy integrated development environment (IDE) for R. 🌟

There are plenty of guides available to help you obtain or update R and RStudio. Here are a couple of them to get you started:

- [R Basics for Paleoecologists](https://ckiahtipes.github.io/) by C.A. Kiahtipes, a previous part of the APD series of workshops.
- [Install or Update R tutorial](https://jennhuck.github.io/workshops/install_update_R.html) by Jennifer Huck.

## Packages 📦

Packages are like magic toolboxes 🧰 that contain a collection of functions for specific needs. We want to make sure that everyone has the necessary packages installed for this workshop.

### Install packages

Let’s create a list of packages that we’ll need from CRAN. Here they are:

``` r
package_list <-
  c(
    "tidyverse", # general data wrangling and visualisation ✨
    "pander", # nice tables 😍
    "Bchron", # age-depth modelling 🕰️
    "janitor", # string cleaning 🧹
    "remotes", # installing packages from GitHub 🚀
    "neotoma2", # access to the Neotoma database 🌿
    "here" # for working directory 🗺️
  )
```

Now, let’s install all these amazing packages from CRAN:

``` r
lapply(
  package_list, utils::install.packages
)
```

#### Install the RRatepol package 🌼

The {RRatepol} is an R package for estimating rate of change (RoC) from community data in time series. Take a peek at its [website](https://hope-uib-bio.github.io/R-Ratepol-package/) for more information.

Let’s try installing the package from GitHub:

``` r
# Install R-Ratepol
remotes::install_github("HOPE-UIB-BIO/R-Ratepol-package")
```

### Test if everything is set up ✅

Let’s do a quick test to make sure everything is in order. Running the following code should produce `"Everything is good to go"` instead of an error message saying `"All required packages are not installed"`.

``` r
if (
  isTRUE(
    all(
      c(package_list, "RRatepol") %in%
        as.data.frame(
          utils::installed.packages()
        )[, 1]
    )
  )
) {
  cat("Everything is good to go")
} else {
  warning("All required packages are not installed")
}
```
