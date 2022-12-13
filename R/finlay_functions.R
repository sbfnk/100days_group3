## get cleaned linelist
get_ll_clean <- function() {

  import(here("data-raw/data_linelist.csv")) %>%
    as_tibble() %>%
    mutate(
      across(c(hospitalisation_date, death_date, date_first_contact, date_last_contact), dmy),
      notification_date = ymd(notification_date),
      onset_date = notification_date - notification_delay_in_days,
      gender = replace(gender, gender %in% c("1", "f", "F", "FF"), "female"),
      gender = replace(gender, gender %in% c("0", "man"), "male"),
      rural_exposure = c("yes" = TRUE, "No" = FALSE, "yEs" = TRUE)[rural_exposure],
      is_farmer = grepl("farm", job)
    )

}

## get clean contacts
get_contacts_clean <- function() {

  import(here("data-raw/Contact tracing data.csv")) %>%
    as_tibble() %>%
    mutate(was_case = c("N" = FALSE, "Y" = TRUE)[was_case])

}

## plot saving function
save_plot <- function(p, file,
                      width = 6.4, height = 4.2,
                      dpi = 300,
                      folder = "figures",
                      ...) {

  if(is.null(p)) return(NULL)

  if(!file.exists(here(folder))) dir.create(here(folder))

  ggsave(
    filename = here(folder, file),
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )

}

## show farmer
vis_farmer <- function(ll) {

  ggplot(ll, aes(onset_date, fill = rural_exposure | is_farmer)) +
    geom_histogram(binwidth = 1, color = "black") +
    scale_fill_manual(values = c("orange", "purple")) +
    labs(x = "Date of onset", y = "Number of cases", fill = "Farmer or rural exposure") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = 'bottom',
      plot.background = element_rect(fill = "white", color = NA)
    )

}

## show phylogeny
vis_tree <- function(tree) {

  p_tree <- ggtree(tree) +
    theme_tree2() +
    geom_tiplab()

}

## A vectorised function to count the occurences of x in y
count_in <- function(x, y) {
  .f <- function(i, y) sum(y == i, na.rm = TRUE)
  return(sapply(x, .f, y))
}

## offspring distribution
vis_offspring <- function(epi) {

  infections <- filter(epi$contacts, was_case)
  count_in(infections$to, infections$from) %>%
    table() %>%
    as.data.frame() %>%
    ggplot(aes(., Freq)) +
    geom_col() +
    theme_minimal() +
    labs(x = "Number of infections caused", y = "Count") +
    theme(plot.background = element_rect(fill = "white", color = NA))

}

## vis offspring over time
vis_offspring_time <- function(epi) {

  infections <- filter(epi$contacts)
  counts <- count_in(infections$to, infections$from)

  tibble(id = names(counts), counts = unname(counts)) %>%
    left_join(epi$linelist) %>%
    ggplot(aes(onset_date, counts)) +
    geom_point() +
    geom_smooth(method = "loess", span = 0.25) +
    theme_minimal() +
    labs(
      x = "Date", y = "Individual reproduction number",
      caption = "*Recent reproduction number estimates will be an underestimate in an ongoing outbreak") +
    theme(plot.background = element_rect(fill = "white", color = NA))

}
