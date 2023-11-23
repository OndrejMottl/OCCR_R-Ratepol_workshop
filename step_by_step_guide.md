# Step-by-step guide

This workflow should show the full strength of the [*RRatepol package*](https://hope-uib-bio.github.io/R-Ratepol-package/) for working with data types **other** than fossil pollen. It should serve as step-by-step guidance starting from downloading dataset from Neotoma, building age-depth models, to estimating rate-of-change using age uncertainty.

⚠️**This workflow is only meant as an example**: There may be several additional steps for data preparation which should be done to properly implement RRatepol and asssess rate of change for your specific project.

For a step-by-spet workflow with fossil pollen data, please see other package materials, such as [African Polled Database workshop](https://ondrejmottl.github.io/APD_R-Ratepol_workshop/index.html).

Additionally, please see [**FOSSILPOL**](https://hope-uib-bio.github.io/FOSSILPOL-website/), an R-based modular workflow to process multiple fossil pollen records to create a comprehensive, standardised dataset compilation, ready for multi-record and multi-proxy analyses at various spatial and temporal scales.

## Install packages

Please follow the [pre-workshop instructions](https://ondrejmottl.github.io/OCCR_R-Ratepol_workshop/pre_workshop.html) to make sure all packages are installed.

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

Here we have selected the **Chickaree Lake** record (ID = 47613) by Higuera, Philip E. and Dunnette, Paul V.

Reference paper: Dunnette, P.V., P.E. Higuera, K.K. McLauchlan, K.M. Derr, C.E. Briles, and M.H. Keefe. 2014. Biogeochemical impacts of wildfires over four millennia in a Rocky Mountain subalpine watershed. New Phytologist 203(3):900-912. \[DOI: 10.1111/nph.12828\]

``` r
sel_dataset_download <-
  neotoma2::get_downloads(47613)
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
  janitor::clean_names()  %>% 
  # remove the ratio as it is just calcualtion from other variables
  dplyr::select(-carbon_nitrogen)

# make table with units for later
data_units <-
  data_samples %>%
  dplyr::distinct(variablename, units)  %>% 
  dplyr::filter(
    variablename != "Carbon:Nitrogen"
  )

head(data_community)[, 1:5]
```

| sample_id | nitrogen | carbon | d15n  |  d13c  |
|:---------:|:--------:|:------:|:-----:|:------:|
|  439811   |    2     | 19.07  | 0.08  | -27.33 |
|  439812   |   1.4    | 15.15  | -0.13 | -26.89 |
|  439813   |   1.43   | 16.67  | -0.23 | -25.73 |
|  439814   |   1.28   | 14.01  | -0.3  | -26.03 |
|  439815   |   1.19   | 12.72  | -0.19 | -25.91 |
|  439816   |   1.15   | 12.94  | 0.04  | -26.78 |

Here, we strongly advocate that careful preparation may of the datasets (with additional steps) may be needed before using R-Ratepol!

We can now try to visualise the taxa per sample_id

``` r
data_community %>%
  tibble::rowid_to_column("ID") %>%
  tidyr::pivot_longer(
    cols = -c(sample_id, ID),
    names_to = "element",
    values_to = "value"
  ) %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = ID,
      y = value,
      col = element
    ),
  ) +
  ggplot2::facet_wrap(
    ~ element,
    nrow = 1,
    scales = "free_x"
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = "grey",
      linewidth = 0.1
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(), 
    legend.position = "none"
  ) +
  ggplot2::labs(
    x = "sample id",
    y = "value"
  )  +
  ggplot2::geom_line()  
```

<img src="step_by_step_guide_files/figure-commonmark/geochemistry%20visualise-1.png" data-fig-align="center" />

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
|  439811   |  0.8  |
|  439812   |  2.9  |
|  439813   |  5.1  |
|  439814   |  7.8  |
|  439815   | 10.5  |
|  439816   | 12.7  |

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
| 27028  |    33168     |  0.8  |    1.6    |    -59.62     |     103110     |
| 27028  |    33168     |  2.4  |    1.1    |    -55.24     |     103111     |
| 27028  |    33168     |  3.5  |    1.1    |      -51      |     103112     |
| 27028  |    33168     |  4.6  |    2.1    |    -45.22     |     103113     |
| 27028  |    33168     |  6.7  |    2.2    |    -39.53     |     103114     |
| 27028  |    33168     |  8.9  |    2.1    |    -33.58     |     103115     |

Table continues below

| agelimityounger | chroncontrolage | chroncontroltype |
|:---------------:|:---------------:|:----------------:|
|     -60.38      |       -60       |     Lead-210     |
|     -56.04      |     -55.64      |     Lead-210     |
|     -51.84      |     -51.42      |     Lead-210     |
|     -46.14      |     -45.68      |     Lead-210     |
|     -40.55      |     -40.04      |     Lead-210     |
|     -34.74      |     -34.16      |     Lead-210     |

``` r
# prepare the table
data_chron_control_table <-
  sel_chron_control_table_download %>%
  # Here select the ID of one of the chronology
  dplyr::filter(chronologyid == 33168) %>%
  tibble::as_tibble() %>%
  # Here we calculate the error as the average of the age `limitolder` and
  #   `agelimityounger`
  dplyr::mutate(
    error = round((agelimitolder - agelimityounger) / 2)
  ) %>%
  # As Bchron cannot accept a error of 0, we need to replace the value with 1
  dplyr::mutate(
    error = replace(error, error == 0, 1),
    error = ifelse(is.na(error), 1, error)
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

| chroncontrolage | error | depth | thickness | chroncontroltype | curve  |
|:---------------:|:-----:|:-----:|:---------:|:----------------:|:------:|
|       -60       |   1   |  0.8  |    1.6    |     Lead-210     | normal |
|     -55.64      |   1   |  2.4  |    1.1    |     Lead-210     | normal |
|     -51.42      |   1   |  3.5  |    1.1    |     Lead-210     | normal |
|     -45.68      |   1   |  4.6  |    2.1    |     Lead-210     | normal |
|     -40.04      |   1   |  6.7  |    2.2    |     Lead-210     | normal |
|     -34.16      |   1   |  8.9  |    2.1    |     Lead-210     | normal |

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

<img src="step_by_step_guide_files/figure-commonmark/plot%20Bchron-1.png" data-fig-align="center" />

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

head(age_uncertainties, n = 8)[, 1:8]
```

Here we see samples (e.g., 439811, 439812, 439813,…) and their possible ages (age-sequence) with each model iteration (posterior). Each age-sequence is similar but there are differences of tens or hundreds of years. We will call this *the uncertainty matrix*.

| 439811 | 439812 | 439813 | 439814 | 439815 | 439816 | 439817 | 439818 |
|:------:|:------:|:------:|:------:|:------:|:------:|:------:|:------:|
|  -61   |  -53   |  -46   |  -37   |  -29   |  -23   |  -11   |   -8   |
|  -61   |  -52   |  -42   |  -37   |  -30   |  -20   |  -13   |   -8   |
|  -60   |  -53   |  -42   |  -37   |  -28   |  -19   |  -12   |   -6   |
|  -60   |  -55   |  -44   |  -40   |  -28   |  -20   |  -11   |   0    |
|  -60   |  -54   |  -43   |  -36   |  -28   |  -19   |  -12   |   -4   |
|  -61   |  -53   |  -41   |  -36   |  -28   |  -19   |  -12   |   -7   |
|  -61   |  -53   |  -44   |  -37   |  -28   |  -19   |  -15   |   0    |
|  -61   |  -55   |  -44   |  -36   |  -28   |  -20   |  -17   |   -8   |

We can visualise these “possible ages” (age-sequence) of each iteration.

``` r
# create a data.frame for plotting
data_age_uncertainties <-
  age_uncertainties %>%
  as.data.frame() %>%
  tibble::rowid_to_column("ID") %>%
  tidyr::pivot_longer(
    cols = -ID,
    names_to = "sample_id",
    values_to = "age"
  ) %>%
  dplyr::left_join(
    data_levels,
    by = dplyr::join_by(sample_id)
  )
```

Each line is a single potential age-depth model iteration (age-sequence). Green points represent the radiocarbon dates. Horizontal lines are depths of our samples.

``` r
(
  fig_age_uncertainties <-
    data_age_uncertainties %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        x = age,
        y = depth
      )
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        group = ID
      ),
      alpha = 0.05,
      linewidth = 0.1
    ) +
    ggplot2::geom_hline(
      yintercept = data_levels$depth,
      lty = 2,
      color = "gray50",
      alpha = 0.1,
      linewidth = 0.1
    ) +
    ggplot2::geom_point(
      data = data_chron_control_table,
      mapping = ggplot2::aes(
        x = chroncontrolage
      ),
      color = "green",
      shape = 15,
      size = 3
    ) +
    ggplot2::scale_y_continuous(trans = "reverse") +
    ggplot2::scale_x_continuous(trans = "reverse") +
    ggplot2::labs(
      x = "age (cal yr BP)"
    )
)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20uncertainty%20matrix-1.png" data-fig-align="center" />

We can visualise all age-depth “possible ages” together as the range of values. Here, each line representing one sampled depth in our record.

``` r
data_age_uncertainties %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = age,
      y = depth,
      group = depth
    )
  ) +
  ggplot2::geom_hline(
    yintercept = data_levels$depth,
    lty = 2,
    color = "gray50",
    alpha = 0.1,
    linewidth = 0.1
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA
  ) +
  ggplot2::scale_y_continuous(trans = "reverse") +
  ggplot2::scale_x_continuous(trans = "reverse") +
  ggplot2::labs(
    x = "age (cal yr BP)"
  )
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20uncertainty%20matrix%20boxplots-1.png" data-fig-align="center" />

Let’s take the median age of all possible ages (i.e. the estimated age from each age-depth model run) as our default.

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
|  439811   |  0.8  | -60 |
|  439812   |  2.9  | -54 |
|  439813   |  5.1  | -44 |
|  439814   |  7.8  | -37 |
|  439815   | 10.5  | -28 |
|  439816   | 12.7  | -21 |

We can visualise the median age by drawing a red line. This age is the age that is often reported in publications but in essence it represents multiple age-depth model runs with smaller or larger age uncertainties throughout the record.

``` r
fig_age_uncertainties +
  ggplot2::geom_point(
    data = data_levels_predicted,
    color = "red",
    size = 3
  ) +
  ggplot2::geom_line(
    data = data_levels_predicted,
    color = "red",
    linewidth = 0.5
  )
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20median%20age-1.png" data-fig-align="center" />

### Visualisation of our data

Let’s now make a simple pollen diagram with proportions of the main pollen taxa (x-axis) against our estimated ages along depth (y-axis).

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
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = "grey",
      linewidth = 0.1
    ),
    legend.position = "none"
  ) +
  ggplot2::labs(
    x = "age (cal yr BP)",
    y = "value"
  ) +
  ggplot2::geom_line()
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20data%20with%20ages-1.png" data-fig-align="center" />

## Estimation Rate-of-Change

Now we will use our prepared geochemistry data and age-depth model to estimate the rate of change. We will present several scenarios (i.e. approaches) to calculate RoC.

### Selection of dissimilarity coefficient

We can check the units of individual measured values:

``` r
data_units
```

| variablename |       units        |
|:------------:|:------------------:|
|   Nitrogen   |      percent       |
|    Carbon    |      percent       |
|     δ15N     |  per mille air N   |
|     δ13C     | per mille VPDB/17O |

As we can see, all measured values are in difrent units. Therefore,for all scenarios, we will be using the `gower` dissimilarity coefficient (works with data in various units), and `time_standardisation` == 250 (this means that all ROC values are ‘change per 250 yr’).

### Scenario 1 - Estimating RoC for each level

This is the “Classic” approach that uses each sampled depth in a record (i.e. individual level) to estimate RoC.

``` r
scenario_1 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "gower",
    time_standardisation = 250,
    working_units = "levels" # here is set to use individual levels
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_1)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20scenario%201-1.png" data-fig-align="center" />

### Scenario 2 - Estimating RoC for each level with smoothing of data

We do the same as in Scenario 1 but now we smooth the community data before calculating RoC. This may be usefull to mitigate the error of measumerments. Specifically, we will add `smooth_method` = “shep” (i.e. Shepard’s 5-term filter).

``` r
scenario_2 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "gower",
    time_standardisation = 250,
    working_units = "levels",
    smooth_method = "shep" # Shepard's 5-term filter
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_2)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20scenario%202-1.png" data-fig-align="center" />

We see that the absolute RoC scores are similar ut the pattern changed slightly (x-axis).

### Scenario 3 - Estimating RoC per bin

In order to get rid of the effect of uneven distribution of sampled depths (i.e. levels) in a record, we can bin the data.

Specifically, we will change the `working_units` from single levels to `"bins"`. Here we select bins of 250 years each instead of the individual levels.

Note that one level is randomly selected as a representation of that time bin. Because of that, we need to increase the number of randomisations. This is again a toy example for a quick computation and therefore we only do 100 randomisations. We would recommend increasing the *set_randomisations* to 10.000 for any real estimation.

``` r
set_randomisations <- 100
```

To speed the process up, you can also set `use_parallel` == `TRUE`, which will use all cores of your computer.

``` r
scenario_3 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "gower",
    working_units = "bins", # change the "bins"
    bin_size = 250, # size of a time bin
    time_standardisation = 250,
    smooth_method = "shep",
    rand = set_randomisations,  # set number of randomisations
    use_parallel = TRUE # do use parallel computing
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_3)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20scenario%203-1.png" data-fig-align="center" />

We will now also visualize uncertainty around the RoC scores shown by a grey shadow.We see a substantial increase in temporal uncertainty around the RoC scores (grey shadow), indicating a loss of temporal precision.

### Scenario 4 - Estimating RoC per bin and calculating age uncertainties

For RoC analysis, it is important to consider age uncertainties. For each iteration, RRatepol will randomly select one age-sequence from the uncertainty matrix (see the age-depth modelling section for more info).

``` r
scenario_4 <-
  RRatepol::estimate_roc(
  data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "gower",
    working_units = "bins", 
    bin_size = 250, 
    time_standardisation = 250,
    smooth_method = "shep",
    rand = set_randomisations, 
    use_parallel = TRUE,
    age_uncertainty = age_uncertainties # Add the uncertainty matrix
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_4)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20scenario%204-1.png" data-fig-align="center" />

Here zou can see that the pattern chnage only slightly. This is because we are randomly sampling age with a small number of randomisations.

### Scenario 5 - Estimating RoC with the new “Moving-window” approach

In order to reduce the temporal uncertainty and improve temporal precision, we can apply a novel approach in RRATEPOL called “moving window”.

``` r
scenario_5 <-
  RRatepol::estimate_roc(
    data_source_community = data_community,
    data_source_age = data_levels_predicted,
    dissimilarity_coefficient = "gower",
    working_units = "MW", # change the "MW" to apply the "moving window"
    bin_size = 250,
    number_of_shifts = 5, # number of shifts
    time_standardisation = 250,
    smooth_method = "shep",
    rand = set_randomisations,
    use_parallel = TRUE,
    age_uncertainty = age_uncertainties
  )
```

``` r
RRatepol::plot_roc(data_source = scenario_5)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20scenario%205-1.png" data-fig-align="center" />

### Scenario 6 - Detecting peak points

Throughout the record, there can be periods when the RoC will substantially change. We can detect RoC increases that are significant by identifying so called *peak-points*. Here, we will use “Non-linear” method, which will detect significant change from a non-linear trend of RoC.

``` r
scenario_5_peak <-
  RRatepol::detect_peak_points(
    data_source = scenario_5,
    sel_method = "trend_non_linear"
  )
```

Now we will plot the RoC estimates showing the peak-points. So here we can see that there were rates of vegetation change throughout the record but only at certain moments in time (green dots - peak points) these changes were significant. There you go!

``` r
RRatepol::plot_roc(
  data_source = scenario_5_peak,
  peaks = TRUE
)
```

<img src="step_by_step_guide_files/figure-commonmark/plot%20peak%20points-1.png" data-fig-align="center" />
