if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(readr, here, dplyr, tidyr, lubridate, ggplot2)
pacman::p_load_gh("epiforecasts/EpiNow2@develop") ## get development version

source(here::here("R", "fix_onset_date.R"))
## code to create snapshots -- this will go into EpiNow2 eventually
source(paste0(
  "https://raw.githubusercontent.com/epiforecasts/nowcasting.example/",
  "main/R/utils.r"
))

## plot data by notification date
ll |>
  count(date = notification_date) |>
  ggplot(aes(x = date, y = n)) + geom_col()

## load serioal interval (in Finlay's script)
si <- readRDS(here::here("data-raw", "si.rds")) |>
  na.omit()

## fit gamma distribution to observed serial intervals (for EpiNow2)
si_dist_file <- here::here("data-raw", "si_dist.rds")
if (!file.exists(si_dist_file)) {
  si_dist <- bootstrapped_dist_fit(si, dist = "gamma")
  saveRDS(si_dist, here::here("data-raw", "si_dist.rds"))
} else {
  si_dist <- readRDS(si_dist_file)
}

## estimate right truncation of onset dates
max_trunc <- max(ll$notification_delay_in_days)
snapshots <- ll |>
  create_snapshots(
    max_delay = max_trunc,
    second = "notification_date",
    first = "onset_date",
    continuous = TRUE
  )
trunc <- estimate_truncation(
  snapshots,
  trunc_max = max_trunc,
  chains = 2, iter = 2000,
  verbose = FALSE
)
trunc_dist <- trunc$dist

## generate count data by onset date
df <- ll |>
  count(date = onset_date, name = "confirm")

## estimate R and make forecast
inf <- estimate_infections(
  reported_cases = df,
  generation_time = si_dist,
  truncation = trunc_opts(trunc_dist),
  stan = stan_opts(cores = 2, chains = 2),
  rt = rt_opts(rw = 7, future = "estimate"),
  gp = NULL,
  verbose = FALSE
)

## extract onset to admission delays
onset_to_admissions <- ll |>
  filter(!is.na(hospitalisation_date)) |>
  mutate(onset_to_admissions = as.integer(hospitalisation_date - onset_date)) |>
  filter(onset_to_admissions <= 14) |> ## exclude strange long delays
  pull(onset_to_admissions)

## estimate distribution from onset to admission
onset_to_admissions_dist_file <- here::here("data-raw", "ota_dist.rds")
if (!file.exists(onset_to_admissions_dist_file)) {
  onset_to_admission_dist <- bootstrapped_dist_fit(
    onset_to_admissions, dist = "gamma"
  )
  saveRDS(onset_to_admission_dist, onset_to_admission_dist_file)
} else {
  onset_to_admission_dist <- readRDS(onset_to_admissions_dist_file)
}

## generate admissions time series
admissions <- ll |>
  count(date = hospitalisation_date, name = "hospitalisations")

## estimate admissions truncation
ll <- ll |>
  mutate(
    admission_to_notification =
      as.integer(notification_date - hospitalisation_date)
  )
max_admissions_trunc <- max(ll$admission_to_notification, na.rm = TRUE)
adm_snapshots <- ll |>
  filter(!is.na(hospitalisation_date)) |>
  create_snapshots(
    max_delay = max_admissions_trunc,
    second = "notification_date",
    first = "hospitalisation_date",
    continuous = TRUE
  )
adm_trunc <- estimate_truncation(
  adm_snapshots,
  trunc_max = max_admissions_trunc,
  chains = 2, iter = 2000,
  verbose = FALSE
)

## create data frame for estimate_secondary
df_sec <- df |>
  left_join(admissions, by = "date") |>
  rename(primary = "confirm", secondary = "hospitalisations") |>
  replace_na(list(secondary = 0)) |>
  filter(cumsum(secondary) > 5) ## remove small numbers

## proportion admitted prior
onsets_admitted <- sum(df_sec$secondary) /
  sum(df_sec$primary)

## estimate relationship between onsets and admissions
sec <- estimate_secondary(
  df_sec,
  delays = delay_opts(onset_to_admission_dist),
  truncation = trunc_opts(adm_trunc$dist),
  secondary = secondary_opts(type = "incidence"),
  obs = obs_opts(
    scale = list(mean = onsets_admitted, sd = onsets_admitted),
    family = "negbin",
    week_effect = TRUE
  ),
  verbose = FALSE
)

## forecast admissions
fc <- forecast_secondary(
  estimate = sec,
  primary = inf,
)

## plot forecasts
plot(fc)
