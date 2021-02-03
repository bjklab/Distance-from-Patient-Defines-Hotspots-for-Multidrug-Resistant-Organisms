#' ###########################################
#' load libraries and set seed
#' ###########################################

library(tidyverse)
library(tidybayes)
library(brms)
library(ggtext)
set.seed(16)



#' ###########################################
#' load data
#' ###########################################

read_csv("./tabs/iceman_cx_data_covid_correction.csv") %>%
  identity() -> iceman_corrected
iceman_corrected



vroom::vroom("./tabs/iceman_asv_data.csv.bz2") %>%
  identity() -> meta_asv
meta_asv



iceman_corrected %>%
  mutate(genus_label = case_when(genus == "Acinetobacter" ~ "<i>Acinetobacter</i> (CR)",
                                 genus == "Citrobacter" ~ "<i>Citrobacter</i> (ESBL/CR)",
                                 genus == "Clostridioides" ~ "<i>Clostridioides difficile</i>",
                                 genus == "Enterobacter" ~ "<i>Enterobacter</i> (ESBL/CR)",
                                 genus == "Enterococcus" ~ "<i>Enterococcus</i> (VRE)",
                                 genus == "Escherichia" ~ "<i>Escherichia</i> (ESBL/CR)",
                                 genus == "Klebsiella" ~ "<i>Klebsiella</i> (ESBL/CR)",
                                 genus == "Pseudomonas" ~ "<i>Pseudomonas</i> (CR)",
                                 genus == "Raoultella" ~ "<i>Raoultella</i> (ESBL/CR)",
                                 genus == "Staphylococcus" ~ "<i>Staphylococcus</i> (MRSA)",
                                 genus == "Stenotrophomonas" ~ "<i>Stenotrophomonas</i> (ESBL/CR)")) %>%
  identity() -> dat
dat



#' ###########################################
#' 
#' model culture positivity as a result of GP(mean distance)
#' 
#' ###########################################


#' get prior
dat %>%
  brms::get_prior(data = ., family = bernoulli,
                  cx_positive ~ gp(distance_meters, by = genus)
  ) %>%
  gt::gt()


#' visualize priors
tribble(
  ~ dist,      ~ args,
  "gamma",      list(1.49, 0.05),
  "student_t", list(3, 0, 2.5)
) %>%
  ggplot(aes(y = dist, dist = dist, args = args)) +
  stat_dist_halfeye() +
  ggtitle("Prior Distributions via stat_dist_halfeye()")



#' run model
# dat %>%
#   brm(formula = cx_positive ~ gp(mean_distance_meters, by = genus),
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 18),
#       backend = "cmdstanr",
#       seed = 16) -> m_mean_distance_gp_mix_genus
# 
# m_mean_distance_gp_mix_genus %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus.rds.gz", compress = "gz")
m_mean_distance_gp_mix_genus <- read_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus.rds.gz")

m_mean_distance_gp_mix_genus
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_genus$fit)
m_mean_distance_gp_mix_genus %>% pp_check()

m_mean_distance_gp_mix_genus %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_genus$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         genus = unique(genus)) %>%
  add_fitted_draws(m_mean_distance_gp_mix_genus) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_mean_distance_gp_mix_genus_fitted
m_mean_distance_gp_mix_genus_fitted


m_mean_distance_gp_mix_genus_fitted %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  facet_wrap(facets = ~ genus_label, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_gp
p_cx_gp


# p_cx_gp %>%
#   ggsave(filename = "./figs/p_cx_gp.pdf", height = 5, width = 8, units = "in")
# p_cx_gp %>%
#   ggsave(filename = "./figs/p_cx_gp.png", height = 5, width = 8, units = "in", dpi = 600)
# p_cx_gp %>%
#   ggsave(filename = "./figs/p_cx_gp.svg", height = 5, width = 8, units = "in")






