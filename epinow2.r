if(!require(pacman)) install.packages("pacman")
pacman::p_load(readr, here, dplyr, tidyr, lubridate, ggplot2)
pacman::p_load_gh("epiforecasts/EpiNow2@develop") ## get development version

source(here::here("R", "fix_onset_date.R"))
## code to create snapshots
source("https://raw.githubusercontent.com/epiforecasts/nowcasting.example/main/R/utils.r")

ll |>
  count(date = notification_date) |>
  ggplot(aes(x = date, y = n)) + geom_col()

si <- readRDS(here::here("data-raw", "si.rds")) |>
  na.omit()

## si_dist <- bootstrapped_dist_fit(si, dist = "gamma")
## saveRDS(si_dist, here::here("data-raw", "si_dist.rds"))
si_dist <- readRDS(here::here("data-raw", "si_dist.rds"))

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

df <- ll |>
  count(date = onset_date, name = "confirm")

inf <- estimate_infections(
  reported_cases = df,
  generation_time = si_dist,
  truncation = trunc_opts(trunc_dist),
  stan = stan_opts(cores = 2, chains = 2),
  rt = rt_opts(rw = 7, future = "estimate"),
  gp = NULL,
  verbose = FALSE
)

onset_to_admissions <- ll |>
  filter(!is.na(hospitalisation_date)) |>
  mutate(onset_to_admissions = as.integer(hospitalisation_date - onset_date)) |>
  filter(onset_to_admissions <= 14) |>
  pull(onset_to_admissions)
## onset_to_admission_dist <- bootstrapped_dist_fit(onset_to_admissions, dist = "gamma")
## saveRDS(onset_to_admission_dist, here::here("data-raw", "ota_dist.rds"))
onset_to_admission_dist <- readRDS(here::here("data-raw", "ota_dist.rds"))
admissions <- ll |>
  count(date = hospitalisation_date, name = "hospitalisations")

ll <- ll |>
  mutate(admission_to_notification = as.integer(notification_date - hospitalisation_date))
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

df_sec <- df |>
  left_join(admissions, by = "date") |>
  rename(primary = "confirm", secondary = "hospitalisations") |>
  replace_na(list(secondary = 0)) |>
  filter(cumsum(secondary) > 5) ## remove small numbers

onsets_admitted <- sum(df_sec$secondary) /
  sum(df_sec$primary)

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

fc <- forecast_secondary(
  estimate = sec,
  primary = inf,
)
