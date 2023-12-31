# Part 2 - XRF data

This workflow should show the full strength of the [*RRatepol package*](https://hope-uib-bio.github.io/R-Ratepol-package/) for working with data types **other** than fossil pollen, specifically for **XRF data**. It should serve as step-by-step guidance starting from downloading datasets from the Neotoma and building age-depth models, to estimating rate-of-change using age uncertainty.

For a workflow using **Geochemical data**, see [Part 1 - Geochemical data](https://ondrejmottl.github.io/OCCR_R-Ratepol_workshop/part_2_xrf_data.html).

⚠️**This workflow is only meant as an example**: There may be several additional steps for data preparation, which should be done to properly implement RRatepol and assess the rate of change for your specific project.

## Prerequisites

Please follow the [pre-workshop instructions](https://ondrejmottl.github.io/OCCR_R-Ratepol_workshop/pre_workshop.html) to make sure all packages are installed.

This is the second part of the workshop. Therefore, some steps will be omitted and/or not explained. For all details about individual steps, see [Part 1 - Geochemical data](https://ondrejmottl.github.io/OCCR_R-Ratepol_workshop/part_2_xrf_data.html).

## Attach packages

``` r
library(tidyverse) # general data wrangling and visualisation ✨
library(pander) # nice tables 😍
library(RRatepol) # rate-of-vegetation change ! > v1.2.0 ! 📈
library(neotoma2) # access to the Neotoma database 🌿
library(Bchron) # age-depth modelingng 🕰️
library(janitor) # string cleaning 🧹
library(here) # for working directory 🗺️
```

## Download a dataset from Neotoma

Here we have selected the **Dářko peat** record (ID = 48306) by Kletetschka, Günther.

Reference paper: Roleček, J., H. Svitavská Svobodová, E. Jamrichová, L. Dudová, P. Hájková, G. Kletetschka, P. Kuneš, and V. Abraham. 2020. Conservation targets from the perspective of a palaeoecological reconstruction: the case study of Dářko peat bog in the Czech Republic. Preslia 92(2):87-114. DOI: 10.23855/preslia.2020.087

``` r
sel_dataset_download <-
  neotoma2::get_downloads(48306)
```

## Prepare the geochemical data

``` r
# get samples
data_samples <-
  neotoma2::samples(sel_dataset_download)  %>% 
  tibble::as_tibble()

# prepare taxa table
data_community <-
  data_samples %>%
  dplyr::mutate(sample_id = as.character(sampleid)) %>%
  dplyr::select("sample_id", "value", "variablename") %>%
  # turn into the wider format
  tidyr::pivot_wider(
    names_from = "variablename",
    values_from = "value",
    values_fill = 0
  ) %>%
  # clean names
  janitor::clean_names()  

# make table with units for later
data_units <-
  data_samples %>%
  dplyr::distinct(variablename, units) 

head(data_community)[, 1:5]
```

| sample_id | strontium | yttrium | zirconium | niobium |
|:---------:|:---------:|:-------:|:---------:|:-------:|
|  445127   |   4e-04   |  9e-04  |  0.0013   | 0.0032  |
|  445128   |   7e-04   | 0.0015  |  0.0024   | 0.0035  |
|  445129   |   8e-04   | 0.0017  |  0.0021   | 0.0037  |
|  445130   |   7e-04   | 0.0015  |  0.0026   | 0.0037  |
|  445131   |   7e-04   | 0.0017  |  0.0023   | 0.0039  |
|  445132   |   8e-04   | 0.0016  |  0.0022   | 0.0041  |

Here, we strongly advocate that careful preparation of the datasets (with additional steps) may be needed before using R-Ratepol!

## Preparation of the levels

### Sample depth

Extract depth for each level

``` r
data_levels <-
  neotoma2::samples(sel_dataset_download) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(sample_id = as.character(sampleid)) %>%
  dplyr::distinct(sample_id, depth) %>%
  dplyr::relocate(sample_id)

head(data_levels)
```

| sample_id | depth |
|:---------:|:-----:|
|  445127   |   0   |
|  445128   |  0.2  |
|  445129   |  0.4  |
|  445130   |  0.6  |
|  445131   |  0.8  |
|  445132   |   1   |

### Age-depth modelling

We will recalculate the age-depth model ‘de novo’ using the [*Bchron* package](http://andrewcparnell.github.io/Bchron/).

#### Prepare chron.control table and run Bchron

The chronology control table contains all the dates (mostly radiocarbon) to create the age-depth model.

Here we only present a few of the important steps of preparation of the chronology control table. There are many more potential issues, but solving those is not the focus of this workflow.

``` r
# First, get the chronologies and check which we want to use used
sel_chron_control_table_download <-
  neotoma2::chroncontrols(sel_dataset_download)

print(sel_chron_control_table_download)
```

| siteid | chronologyid | depth | thickness | agelimitolder | chroncontrolid |
|:------:|:------------:|:-----:|:---------:|:-------------:|:--------------:|
| 27171  |    33682     |   0   |    NA     |      -60      |     104411     |
| 27171  |    33682     |  25   |     1     |      645      |     104412     |
| 27171  |    33682     |  50   |     1     |     1795      |     104413     |
| 27171  |    33682     |  100  |     1     |     2400      |     104414     |
| 27171  |    33682     |  155  |     1     |     2810      |     104415     |
| 27171  |    33682     |  200  |     1     |     3590      |     104416     |

Table continues below

| agelimityounger | chroncontrolage | chroncontroltype |
|:---------------:|:---------------:|:----------------:|
|       -70       |       -65       |     Core top     |
|       585       |       615       |   Radiocarbon    |
|      1735       |      1765       |   Radiocarbon    |
|      2340       |      2370       |   Radiocarbon    |
|      2740       |      2775       |   Radiocarbon    |
|      3530       |      3560       |   Radiocarbon    |

``` r
# prepare the table
data_chron_control_table <-
  sel_chron_control_table_download %>%
  # Here select the ID of one of the chronology
  dplyr::filter(chronologyid == 33682) %>%
  tibble::as_tibble() %>%
  # Here we calculate the error as the average of the age `limitolder` and
  #   `agelimityounger`
  dplyr::mutate(
    error = round((agelimitolder - agelimityounger) / 2)
  ) %>%
  # As Bchron cannot accept a error of 0, we need to replace the value with 1
  dplyr::mutate(
    error = replace(error, error == 0, 1),
    error = ifelse(is.na(error), 1, error),
    thickness = ifelse(is.na(thickness), 1, thickness)
  ) %>%
  # We need to specify which calibration curve should be used for what point
  dplyr::mutate(
    curve = ifelse(as.data.frame(sel_dataset_download)["lat"] > 0, "intcal20", "shcal20"),
    curve = ifelse(chroncontroltype != "Radiocarbon", "normal", curve)
  ) %>%
  tibble::column_to_rownames("chroncontrolid") %>%
  dplyr::arrange(depth) %>%
  dplyr::select(
    chroncontrolage, error, depth, thickness, chroncontroltype, curve
  )

head(data_chron_control_table)
```

| chroncontrolage | error | depth | thickness | chroncontroltype |  curve   |
|:---------------:|:-----:|:-----:|:---------:|:----------------:|:--------:|
|       -65       |   5   |   0   |     1     |     Core top     |  normal  |
|       615       |  30   |  25   |     1     |   Radiocarbon    | intcal20 |
|      1765       |  30   |  50   |     1     |   Radiocarbon    | intcal20 |
|      2370       |  30   |  100  |     1     |   Radiocarbon    | intcal20 |
|      2775       |  35   |  155  |     1     |   Radiocarbon    | intcal20 |
|      3560       |  30   |  200  |     1     |   Radiocarbon    | intcal20 |

As this is just a toy example, we will use only the iteration multiplier (`i_multiplier`) of `0.1` to reduce the computation time. However, we strongly recommend increasing it to 5 for any normal age-depth model construction.

``` r
i_multiplier <- 0.1 # increase to 5

# Those are default values suggested by the Bchron package
n_iteration_default <- 10e3
n_burn_default <- 2e3
n_thin_default <- 8

# Let's multiply them by our i_multiplier
n_iteration <- n_iteration_default * i_multiplier
n_burn <- n_burn_default * i_multiplier
n_thin <- max(c(1, n_thin_default * i_multiplier))

# run Bchron
sel_bchron <-
  Bchron::Bchronology(
    ages = data_chron_control_table$chroncontrolage,
    ageSds = data_chron_control_table$error,
    positions = data_chron_control_table$depth,
    calCurves = data_chron_control_table$curve,
    positionThicknesses = data_chron_control_table$thickness,
    iterations = n_iteration,
    burn = n_burn,
    thin = n_thin
  )
```

Visually check the age-depth models

``` r
Bchron:::plot.BchronologyRun(sel_bchron) + # or just simple plot(sel_bchron)
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 25),
    axis.text = ggplot2::element_text(size = 15),
    strip.text = ggplot2::element_text(size = 15),
    panel.grid = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = "age (cal yr BP)",
    y = "depth"
  )
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20Bchron-1.png" data-fig-align="center" />

#### Predict ages

Let’s first extract posterior ages (i.e. possible ages) from the age-depth model.

``` r
age_position <-
  Bchron:::predict.BchronologyRun( # or just simple predict(sel_bchron)
    object = sel_bchron,
    newPositions = data_levels$depth
  )

age_uncertainties <-
  age_position %>%
  as.data.frame() %>%
  dplyr::mutate_all(., as.integer) %>%
  as.matrix()

colnames(age_uncertainties) <- data_levels$sample_id
```

``` r
data_levels_predicted <-
  data_levels %>%
  dplyr::mutate(
    age = apply(
      age_uncertainties, 2,
      stats::quantile,
      probs = 0.5
    )
  )

head(data_levels_predicted)
```

| sample_id | depth | age |
|:---------:|:-----:|:---:|
|  445127   |   0   | -64 |
|  445128   |  0.2  | -59 |
|  445129   |  0.4  | -53 |
|  445130   |  0.6  | -48 |
|  445131   |  0.8  | -43 |
|  445132   |   1   | -37 |

### Visualisation of our data

Let’s now make a simple diagram with proportions of the concentrations (x-axis) against our estimated ages along the depth (y-axis).

``` r
data_community %>%
  dplyr::inner_join(
    data_levels_predicted,
    by = dplyr::join_by(sample_id)
  ) %>%
  tidyr::pivot_longer(
    cols = -c(sample_id, depth, age),
    names_to = "element",
    values_to = "value"
  ) %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = age,
      y = value,
      col = element
    ),
  ) +
  ggplot2::facet_wrap(
    ~element,
    nrow = 1,
    scales = "free_x"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_x_continuous(trans = "reverse") +
  ggplot2::scale_y_continuous(breaks = c(0, 50, 100)) +
  ggplot2::theme(
     panel.grid.major.x = ggplot2::element_line(
      colour = "grey",
      linewidth = 0.1
    ),
    legend.position = "none",
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = "age (cal yr BP)",
    y = "Percentage"
  ) +
  ggplot2::geom_line()
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20data%20with%20ages-1.png" data-fig-align="center" />

## Estimation Rate-of-Change

Now we will use our prepared XRF data and age-depth model to estimate the rate of change. We will present several scenarios (i.e. approaches) to calculate RoC.

### Selection of dissimilarity coefficient

We can check the units of individual measured values:

``` r
data_units %>%
  dplyr::pull("units") %>%
  table()
```

| percent |
|:-------:|
|   29    |

As we can see, all measured values are in percentage units. Therefore, for all scenarios, we will be using the `chisq` dissimilarity coefficient (works with closed data), and `time_standardisation` == 500 (this means that all ROC values are ‘change per 500 yr’).

### Scenario 1 - Estimating RoC for each level

This is the “Classic” approach that uses each sampled depth in a record (i.e. individual level) to estimate RoC.

``` r
scenario_1 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "chisq",
    time_standardisation = 500,
    working_units = "levels" # here is set to use individual levels
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_1)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20scenario%201-1.png" data-fig-align="center" />

### Scenario 2 - Estimating RoC for each level with smoothing of data

We do the same as in Scenario 1 but now we smooth the community data before calculating RoC. This may be useful to mitigate the error of measurement. Specifically, we will add `smooth_method` = “shep” (i.e. Shepard’s 5-term filter).

``` r
scenario_2 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "chisq",
    time_standardisation = 500,
    working_units = "levels",
    smooth_method = "shep" # Shepard's 5-term filter
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_2)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20scenario%202-1.png" data-fig-align="center" />

We see that the absolute RoC scores are decreased and the pattern changed slightly (x-axis).

### Scenario 3 - Using uncertainty matrix

For RoC analysis, it is important to consider age uncertainties. For each iteration, RRatepol will randomly select one age-sequence from the uncertainty matrix (see the age-depth modeling section for more info).

Because of that, we need to increase the number of randomizations. This is again a toy example for a quick computation and therefore we only do 100 randomizations. We would recommend increasing the *set_randomisations* to 10.000 for any real estimation.

``` r
set_randomisations <- 100
```

To speed the process up, you can also set `use_parallel` == `TRUE`, which will use all cores of your computer.

``` r
scenario_3 <-
  RRatepol::estimate_roc(
  data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "chisq",
    time_standardisation = 500,
    smooth_method = "shep",
    age_uncertainty = age_uncertainties, # Add the uncertainty matrix
    rand = set_randomisations,  # set number of randomisations
    use_parallel = TRUE # do use parallel computing
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_3)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20scenario%203-1.png" data-fig-align="center" />

We will now also visualize uncertainty around the RoC scores shown by a grey shadow. We see a substantial increase in RoC value in certain regions, this is caused by the extremely small numbers of age differences and the low number of randomizations.

### Scenario 4 - Estimating RoC per bin

In order to get rid of the effect of uneven distribution of sampled depths (i.e. levels) in a record, we can bin the data.

Specifically, we will change the `working_units` from single levels to `"bins"`. Here we select bins of 500 years each instead of the individual levels.

Note that one level is randomly selected as a representation of that time bin. Therefore, it is important to increase the number of randomizations.

``` r
scenario_4 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "chisq",
    working_units = "bins", # change the "bins"
    bin_size = 500, # size of a time bin
    time_standardisation = 500,
    smooth_method = "shep",
    age_uncertainty = age_uncertainties,
    rand = set_randomisations, 
    use_parallel = TRUE 
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_4)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20scenario%204-1.png" data-fig-align="center" />

Here we can see a drastic change in the shape of RoC but a large loss of temporal precision.

### Scenario 5 - Estimating RoC with the new “Moving-window” approach

In order to reduce the temporal uncertainty and improve temporal precision, we can apply a novel approach in RRATEPOL called “moving window”.

``` r
scenario_5 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "chisq",
    working_units = "MW", # change the "MW" to apply the "moving window"
    bin_size = 500,
    number_of_shifts = 5, # number of shifts
    time_standardisation = 500,
    smooth_method = "shep",
    rand = set_randomisations,
    use_parallel = TRUE,
    age_uncertainty = age_uncertainties
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_5)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20scenario%205-1.png" data-fig-align="center" />

### Scenario 6 - Detecting peak points

Throughout the record, there can be periods when the RoC will substantially change. We can detect RoC increases that are significant by identifying so-called *peak-points*. Here, we will use the “Non-linear” method, which will detect a significant change from a non-linear trend of RoC.

``` r
scenario_5_peak <-
  RRatepol::detect_peak_points(
    data_source = scenario_5,
    sel_method = "trend_non_linear"
  )
```

Now we will plot the RoC estimates showing the peak points. So here we can see that there were rates of vegetation change throughout the record but only at certain moments in time (green dots - peak points) these changes were significant. There you go!

``` r
RRatepol::plot_roc(
  data_source = scenario_5_peak,
  peaks = TRUE
)
```

<img src="part_2_xrf_data_files/figure-commonmark/plot%20peak%20points-1.png" data-fig-align="center" />
