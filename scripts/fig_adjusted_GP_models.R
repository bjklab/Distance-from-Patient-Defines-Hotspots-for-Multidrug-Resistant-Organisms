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


meta_asv %>%
  # complete ASVs filtered for 0 reads
  complete(nesting(specimen_id,subject_id,sampleevent,specimen_site), 
           nesting(seqvar_id, taxonomy, kingdom, phylum, class, order, family, genus, species)) %>%
  replace_na(list("read_count" = 0))



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













#' ###########################################
#' 
#' model culture detection as a result of GP(mean distance), after adjustment for sequence abundance
#' 
#' ###########################################

dat_asv_genus %>%
  gt::gt_preview()

dat_asv_genus %>%
  mutate(genus = replace(genus, genus == "Clostridium", "Clostridioides")) %>%
  left_join(select(dat, subject_id, specimen_site, sampleevent, genus, genus_label, cx_positive), by = c("subject_id", "specimen_site", "sampleevent", "genus")) %>%
  # make sure all genera detected by culture (Raoultella not in this subset)
  group_by(genus) %>%
  filter(sum(cx_positive, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  identity() -> dat_cx_asv_genus
dat_cx_asv_genus


#' get prior
dat_cx_asv_genus %>%
  brms::get_prior(data = ., family = bernoulli,
                  cx_positive ~ (1 + clr_read_count|genus) + gp(distance_meters, by = genus)
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
# dat_cx_asv_genus %>%
#   brm(formula = cx_positive ~ (1 + clr_read_count|genus) + gp(mean_distance_meters, by = genus),
#       data = .,
#       family = bernoulli,
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 18),
#       backend = "cmdstanr",
#       seed = 16) -> m_mean_distance_gp_mix_genus_seq_adjusted
# 
# m_mean_distance_gp_mix_genus_seq_adjusted %>% write_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus_seq_adjusted.rds.gz", compress = "gz")
m_mean_distance_gp_mix_genus_seq_adjusted <- read_rds(file = "./models/binomial/m_cx_v_mean_distance_gp_mix_genus_seq_adjusted.rds.gz")

m_mean_distance_gp_mix_genus_seq_adjusted
rstan::check_hmc_diagnostics(m_mean_distance_gp_mix_genus_seq_adjusted$fit)
m_mean_distance_gp_mix_genus_seq_adjusted %>% pp_check()

m_mean_distance_gp_mix_genus_seq_adjusted %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_mix_genus_seq_adjusted$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         genus = unique(genus),
         clr_read_count = 0) %>%
  add_fitted_draws(m_mean_distance_gp_mix_genus_seq_adjusted) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_mean_distance_gp_mix_genus_seq_adjusted_fitted
m_mean_distance_gp_mix_genus_seq_adjusted_fitted


#' fitted plot
m_mean_distance_gp_mix_genus_seq_adjusted_fitted %>%
  mutate(genus_italic = glue::glue("<i>{genus}</i>")) %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  facet_wrap(facets = ~ genus_label, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Purples") +
  labs(x = "Distance from Patient & Towards Bathroom (meters)",
       y = "Adjusted Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_gp_seq_adjusted
p_cx_gp_seq_adjusted


# p_cx_gp_seq_adjusted %>%
#   ggsave(filename = "./figs/p_cx_gp_seq_adjusted.pdf", height = 5, width = 8, units = "in")
# p_cx_gp_seq_adjusted %>%
#   ggsave(filename = "./figs/p_cx_gp_seq_adjusted.png", height = 5, width = 8, units = "in", dpi = 600)
# p_cx_gp_seq_adjusted %>%
#   ggsave(filename = "./figs/p_cx_gp_seq_adjusted.svg", height = 5, width = 8, units = "in")





#' ###########################################
#' 
#' mixed effects for location & genus, culture positivity ~ random slope abundance
#' 
#' ###########################################

dat_asv_genus %>%
  gt::gt_preview()

dat_asv_genus %>%
  mutate(genus = replace(genus, genus == "Clostridium", "Clostridioides")) %>%
  left_join(select(dat, subject_id, specimen_site, sampleevent, genus, genus_label, cx_positive), by = c("subject_id", "specimen_site", "sampleevent", "genus")) %>%
  # make sure all genera detected by culture (Raoultella not in this subset)
  group_by(genus) %>%
  filter(sum(cx_positive, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  identity() -> dat_cx_asv_genus
dat_cx_asv_genus


#' get prior
dat_cx_asv_genus %>%
  brms::get_prior(data = ., family = bernoulli,
                  cx_positive ~ 1 + (0 + clr_read_count | genus / specimen_site)
  ) %>%
  gt::gt()


#' run model
# dat_cx_asv_genus %>%
#   brm(formula = cx_positive ~ 1 + (0 + clr_read_count | genus / specimen_site),
#       data = .,
#       family = bernoulli,
#       #prior = c(set_prior(prior = "student_t(3, 0, 1)", class = "sd")),
#       chains = 4,
#       cores = 4,
#       control = list("adapt_delta" = 0.999, max_treedepth = 18),
#       backend = "cmdstanr",
#       seed = 16) -> m_cx_pos_mixed_site_genus
# 
# m_cx_pos_mixed_site_genus %>% write_rds(file = "./models/binomial/m_cx_pos_mixed_site_genus.rds.gz", compress = "gz")
m_cx_pos_mixed_site_genus <- read_rds(file = "./models/binomial/m_cx_pos_mixed_site_genus.rds.gz")

m_cx_pos_mixed_site_genus
rstan::check_hmc_diagnostics(m_cx_pos_mixed_site_genus$fit)
m_cx_pos_mixed_site_genus %>% pp_check()

m_cx_pos_mixed_site_genus %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' posterior samples

# m_cx_pos_mixed_site_genus %>%
#   spread_draws(r_genus[genus]) %>%
#   mutate(r_genus_OR = exp(r_genus)) %>%
#   ggplot(data = ., aes(y = genus, x = r_genus_OR, fill = stat(x > 1))) +
#   tidybayes::stat_halfeye(alpha = 0.9) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   scale_fill_manual(values = c("gray80","skyblue")) +
#   theme_bw() +
#   theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
#         axis.text.x = ggtext::element_markdown(color = "black"),
#         axis.text.y = ggtext::element_markdown(color = "black"),
#         legend.position = "none",
#         legend.background = element_rect(fill = "white", color = "black", size = 0.25),
#         strip.background = element_blank()) +
#   labs(x = "Marginal OR of MDRO Detection per CLR-Transformed 16S rRNA Abundance (all sites)",
#        y = "") -> p_genus_clr_marginal
# p_genus_clr_marginal

  


m_cx_pos_mixed_site_genus %>%
  spread_draws(`r_genus:specimen_site`[genus_site]) %>%
  rename(r_genus_site = `r_genus:specimen_site`) %>%
  mutate(r_genus_site_OR = exp(r_genus_site)) %>%
  separate(col = genus_site, into = c("genus","site"), sep = "_", remove = TRUE) %>%
  left_join(distinct(select(dat,genus,genus_label)), by = "genus") %>% 
  mutate(genus_label = gsub("</i>", "</i><br>", genus_label)) %>%
  mutate(genus_label = gsub(" difficile","<br>difficile", genus_label)) %>%
  mutate(genus = factor(genus_label),
         genus = reorder(genus_label, desc(genus_label)),
         site = factor(site),
         site = fct_recode(site,
                           "Near Patient" = "ES1",
                           "Intermediate Distance" = "ES2",
                           "Far From Patient, Near Wastewater" = "ES3")) %>%
  ggplot(data = ., aes(y = genus, x = r_genus_site_OR, fill = stat(x > 1))) +
  tidybayes::stat_halfeye(alpha = 0.9) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue")) +
  facet_wrap(facets = ~ site, nrow = 1) +
  coord_cartesian(xlim = c(0,3.5)) +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 9),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "none",
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) +
  labs(x = "Odds Ratio of MDRO Detection as 16S rRNA Sequence Abundance of Genus Increases",
       y = "") -> p_genus_clr_marginal_per_site
p_genus_clr_marginal_per_site



# p_genus_clr_marginal_per_site %>%
#   ggsave(filename = "./figs/p_genus_clr_marginal_per_site.pdf", height = 5, width = 8, units = "in")
# p_genus_clr_marginal_per_site %>%
#   ggsave(filename = "./figs/p_genus_clr_marginal_per_site.png", height = 5, width = 8, units = "in", dpi = 600)
# p_genus_clr_marginal_per_site %>%
#   ggsave(filename = "./figs/p_genus_clr_marginal_per_site.svg", height = 5, width = 8, units = "in")


m_cx_pos_mixed_site_genus %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  filter(param != "lp__") %>%
  mutate_if(.predicate = ~ is.numeric(.x), .f = ~ exp(.x)) %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3) %>%
  gt::data_color(
    columns = vars(Estimate, Q2.5, Q97.5),
    apply_to = "fill",
    colors = scales::col_bin(
      bins = c(-Inf,1,Inf),
      # custom defined values - notice that order matters!
      palette = c("#ffffff", "#35b0ab"), #"#f2fbd2", "#c9ecb4", "#93d3ab",
      domain = NULL
    )
  ) %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")




#' ###########################################
#' 
#' model residuals
#' 
#' ###########################################

m_cx_pos_mixed_site_genus$data %>%
  as_tibble() %>%
  expand(clr_read_count = modelr::seq_range(clr_read_count, n = 10),
         genus = unique(genus),
         specimen_site = unique(specimen_site),
         cx_positive = unique(cx_positive)) %>%
  add_residual_draws(m_cx_pos_mixed_site_genus) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  ungroup() %>%
  mutate(site = factor(specimen_site),
         site = fct_recode(site,
                           "Near Patient" = "ES1",
                           "Intermediate Distance" = "ES2",
                           "Far From Patient, Near Bathroom" = "ES3")) %>%
  ungroup() %>%
  identity() -> m_cx_pos_mixed_site_genus_resid
m_cx_pos_mixed_site_genus_resid

#qplot(data = ., x = genus, y = .residual, fill = genus, geom = "boxplot", facets = ~ site)


