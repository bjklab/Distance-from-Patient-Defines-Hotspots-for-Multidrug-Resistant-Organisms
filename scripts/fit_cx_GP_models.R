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
dat %>%
  brm(formula = cx_positive ~ gp(mean_distance_meters, by = genus),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_mix_genus

m_mean_distance_gp_mix_genus %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus.rds.gz", compress = "gz")
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





#' ###########################################
#' 
#' ADD RANDOM INTERCEPT FOR SUBJECTS & model culture positivity as a result of GP(mean distance)
#' 
#' ###########################################

#' run model
dat %>%
  brm(formula = cx_positive ~ gp(mean_distance_meters, by = genus) + (1|subject_id),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_mix_genus_subject

m_mean_distance_gp_mix_genus_subject %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus_subject.rds.gz", compress = "gz")
m_mean_distance_gp_mix_genus_subject <- read_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus_subject.rds.gz")

m_mean_distance_gp_mix_genus_subject
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_genus_subject$fit)
m_mean_distance_gp_mix_genus_subject %>% pp_check()

m_mean_distance_gp_mix_genus_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_genus_subject$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         genus = unique(genus),
         subject_id = unique(subject_id)) %>%
  add_fitted_draws(m_mean_distance_gp_mix_genus_subject) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_mean_distance_gp_mix_genus_subject_fitted
m_mean_distance_gp_mix_genus_subject_fitted


m_mean_distance_gp_mix_genus_subject_fitted %>%
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
        strip.background = element_blank()) -> p_cx_gp_subject
p_cx_gp_subject


# p_cx_gp_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_subject.pdf", height = 5, width = 8, units = "in")
# p_cx_gp_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_subject.png", height = 5, width = 8, units = "in", dpi = 600)
# p_cx_gp_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_subject.svg", height = 5, width = 8, units = "in")





#' ###########################################
#' 
#' model aggregate by Gram stain as a result of GP(mean distance)
#' 
#' ###########################################

dat %>%
  #count(genus)
  mutate(Gram = case_when(genus %in% c("Acinetobacter", "Citrobacter", "Enterobacter", "Escherichia", "Klebsiella", "Pseudomonas", "Raoultella", "Stenotrophomonas") ~ "Gram negative",
                          genus %in% c("Clostridioides", "Enterococcus", "Staphylococcus") ~ "Gram positive")) %>%
  identity() -> dat_gram
dat_gram



#' run model
dat_gram %>%
  brm(formula = cx_positive ~ gp(mean_distance_meters, by = Gram),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_mix_gram

m_mean_distance_gp_mix_gram %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_gram.rds.gz", compress = "gz")
m_mean_distance_gp_mix_gram <- read_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_gram.rds.gz")

m_mean_distance_gp_mix_gram
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_gram$fit)
m_mean_distance_gp_mix_gram %>% pp_check()

m_mean_distance_gp_mix_gram %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_gram$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         Gram = unique(Gram)) %>%
  add_fitted_draws(m_mean_distance_gp_mix_gram) %>%
  identity() -> m_mean_distance_gp_mix_gram_fitted
m_mean_distance_gp_mix_gram_fitted


m_mean_distance_gp_mix_gram_fitted %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  facet_wrap(facets = ~ Gram, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.1, 0.75),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_gp_gram
p_cx_gp_gram


# p_cx_gp_gram %>%
#   ggsave(filename = "./figs/p_cx_gp_gram.pdf", height = 4, width = 5, units = "in")
# p_cx_gp_gram %>%
#   ggsave(filename = "./figs/p_cx_gp_gram.png", height = 4, width = 5, units = "in", dpi = 600)
# p_cx_gp_gram %>%
#   ggsave(filename = "./figs/p_cx_gp_gram.svg", height = 4, width = 5, units = "in")




#' ###########################################
#' 
#' ADD RANDOM INTERCEPT to Gram stain as a result of GP(mean distance)
#' 
#' ###########################################

#' run model
dat_gram %>%
  brm(formula = cx_positive ~ gp(mean_distance_meters, by = Gram) + (1|subject_id),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_mix_gram_subject

m_mean_distance_gp_mix_gram_subject %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_gram_subject.rds.gz", compress = "gz")
m_mean_distance_gp_mix_gram_subject <- read_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_gram_subject.rds.gz")

m_mean_distance_gp_mix_gram_subject
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_gram_subject$fit)
m_mean_distance_gp_mix_gram_subject %>% pp_check()

m_mean_distance_gp_mix_gram_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_gram_subject$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         Gram = unique(Gram)) %>%
  add_fitted_draws(m_mean_distance_gp_mix_gram) %>%
  identity() -> m_mean_distance_gp_mix_gram_fitted_subject
m_mean_distance_gp_mix_gram_fitted_subject


m_mean_distance_gp_mix_gram_fitted_subject %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  facet_wrap(facets = ~ Gram, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  scale_y_continuous(limits = c(0,NA)) +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 10),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.1, 0.75),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_gp_gram_subject
p_cx_gp_gram_subject


# p_cx_gp_gram_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_gram_subject.pdf", height = 4, width = 5, units = "in")
# p_cx_gp_gram_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_gram_subject.png", height = 4, width = 5, units = "in", dpi = 600)
# p_cx_gp_gram_subject %>%
#   ggsave(filename = "./figs/p_cx_gp_gram_subject.svg", height = 4, width = 5, units = "in")






