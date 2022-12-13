library("dplyr")
library("lubridate")

ll <- read_csv(here::here("data-raw", "data_linelist.csv")) %>%
  mutate(across(c(hospitalisation_date, death_date, date_first_contact, date_last_contact), dmy)) %>%
  mutate(notification_date = ymd(notification_date)) %>%
  # There was some data entry issues for onset_date so we re-create it from other columns
  mutate(onset_date = notification_date - notification_delay_in_days)
