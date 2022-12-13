## load packages and functions
if(!require(pacman)) install.packages("pacman")
pacman::p_load_gh("reconhub/epicontacts@timeline")
pacman::p_load(tidyverse, rio, magrittr, here, distcrete, epitrix, treeio, ape,
               epicontacts, ggtree, lubridate)
source(here("R/epicontacts_functions.R"))

## load dna and tree
dna <- read.FASTA(here("data-raw/sequences.fasta"))
tree <- read.tree(here("data-raw/sequences.fasta.treefile"))

## get cleaned linelist
ll <- get_ll_clean()

## get cleaned contacts
contacts <- get_contacts_clean()

## make epicontacts
epi <- make_epicontacts(
  ll, contacts,
  id = "case_name",
  from = "part_name",
  to = "contact_name",
  directed = TRUE
)

## epicontacts with only infection links
ttree <- thin(epi, what = "contacts")
ttree$linelist %<>% mutate(
  R = epicontacts::get_degree(ttree, "out", only_linelist = TRUE)
)

## serial interval
si <- EpiNow2::bootstrapped_dist_fit(get_pairwise(ttree, "onset_date"), dist = "gamma")

## visualise tree
vis_tree(tree) %>%
  save_plot("tree.png")

## visualise proportion farmer
vis_farmer(ll) %>%
  save_plot("farmer.png")

## visualise offspring distribution
vis_offspring(epi) %>%
  save_plot("offspring.png")

## visualise offspring distribution
vis_offspring_time(epi) %>%
  save_plot("offspring_time.png")

## make visnetwork
net <- plot(
  thin(epi, "contacts"),
  x_axis = "onset_date",
  height = 2000,
  arrow_size = 0.1,
  label = FALSE,
  node_size = 8,
  node_color = "is_farmer",
  col_pal = c("TRUE" = "orange", "FALSE" = "purple"),
  ## col_pal = c("female" = "orange", "male" = "purple")
  reverse_root_order = TRUE,
  selector = "id",
  axis_type = "double"
)

## export visnetwork
visNetwork::visSave(net, here("figures/ttree.html"))
