#' ###########################################
#' load libraries and set seed
#' ###########################################

library(tidyverse)
library(tidybayes)
library(brms)
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
  mutate(genus_label = case_when(genus == "Acinetobacter" ~ "<i>Acinetobacter</i> (ESBL/CR)",
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



#' prepare data
meta_asv %>%
  # complete ASVs filtered for 0 reads
  complete(nesting(specimen_id,subject_id,sampleevent,specimen_site), 
           nesting(seqvar_id, taxonomy, kingdom, phylum, class, order, family, genus, species)) %>%
  replace_na(list("read_count" = 0)) %>%
  # recalculate CLR
  mutate(clr_read_count = read_count,
         clr_read_count = replace(clr_read_count, clr_read_count == 0, 1)) %>%
  group_by(specimen_id) %>%
  mutate(clr_read_count = as.vector(compositions::clr(clr_read_count)),
         seq_positive = read_count > 0) %>%
  ungroup() %>%
  # add distance data
  left_join(select(dat, subject_id, specimen_site, distance_meters, mean_distance_meters),
            by = c("subject_id", "specimen_site")) %>%
  distinct() %>%
  identity() -> meta_asv_complete
meta_asv_complete






#' ###########################################
#' 
#' model sequencing detection as a result of GP(mean distance)
#' 
#' ###########################################

#' prepare data
meta_asv %>%
  # standardize genera names
  mutate(genus = gsub("g__|\\]|\\[","",genus)) %>%
  # aggregate reads at genus level
  group_by(specimen_id, subject_id, sampleevent, specimen_site, genus) %>%
  summarise(read_count = sum(read_count, na.rm = TRUE)) %>%
  ungroup() %>%
  # recalculate CLR
  mutate(clr_read_count = read_count,
         clr_read_count = replace(clr_read_count, clr_read_count == 0, 1)) %>%
  group_by(specimen_id) %>%
  mutate(clr_read_count = as.vector(compositions::clr(clr_read_count)),
         seq_positive = read_count > 0) %>%
  ungroup() %>%
  identity() -> meta_asv_genus
meta_asv_genus



meta_asv_genus %>%
  # include Clostridium for Clostridioides difficile
  filter(genus %in% c(unique(iceman_corrected$genus), "Clostridium")) %>%
  #count(genus)
  left_join(select(dat, subject_id, specimen_site, distance_meters, mean_distance_meters),
            by = c("subject_id", "specimen_site")) %>%
  distinct() %>%
  # make sure all genera detected by sequencing
  group_by(genus) %>%
  filter(sum(seq_positive, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  identity() -> dat_asv_genus
dat_asv_genus



#' get prior
dat_asv_genus %>%
  brms::get_prior(data = ., family = bernoulli,
                  seq_positive ~ gp(distance_meters, by = genus)
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
# dat_asv_genus %>%
#   brm(formula = seq_positive ~ gp(mean_distance_meters, by = genus),
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 18),
#       backend = "cmdstanr",
#       seed = 16) -> m_mean_distance_gp_mix_genus_seq
# 
# m_mean_distance_gp_mix_genus_seq %>% write_rds(file = "./models/binomial/m_asv_v_mean_distance_gp_mix_genus.rds.gz", compress = "gz")
m_mean_distance_gp_mix_genus_seq <- read_rds(file = "./models/binomial/m_asv_v_mean_distance_gp_mix_genus.rds.gz")

m_mean_distance_gp_mix_genus_seq
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_genus_seq$fit)
m_mean_distance_gp_mix_genus_seq %>% pp_check()

m_mean_distance_gp_mix_genus_seq %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_genus_seq$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         genus = unique(genus)) %>%
  add_fitted_draws(m_mean_distance_gp_mix_genus_seq) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_mean_distance_gp_mix_genus_seq_fitted
m_mean_distance_gp_mix_genus_seq_fitted


#' fitted
m_mean_distance_gp_mix_genus_seq_fitted %>%
  mutate(genus_italic = glue::glue("<i>{genus}</i>")) %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  facet_wrap(facets = ~ genus_italic, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Probability of Detection by 16S rRNA Sequencing",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_asv_gp
p_asv_gp


# p_asv_gp %>%
#   ggsave(filename = "./figs/p_asv_gp.pdf", height = 5, width = 8, units = "in")
# p_asv_gp %>%
#   ggsave(filename = "./figs/p_asv_gp.png", height = 5, width = 8, units = "in", dpi = 600)
# p_asv_gp %>%
#   ggsave(filename = "./figs/p_asv_gp.svg", height = 5, width = 8, units = "in")







