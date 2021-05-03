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



#' prepare data
# meta_asv %>%
#   # complete ASVs filtered for 0 reads
#   complete(nesting(specimen_id,subject_id,sampleevent,specimen_site), 
#            nesting(seqvar_id, taxonomy, kingdom, phylum, class, order, family, genus, species)) %>%
#   replace_na(list("read_count" = 0)) %>%
#   # recalculate CLR
#   mutate(clr_read_count = read_count,
#          clr_read_count = replace(clr_read_count, clr_read_count == 0, 1)) %>%
#   group_by(specimen_id) %>%
#   mutate(clr_read_count = as.vector(compositions::clr(clr_read_count)),
#          seq_positive = read_count > 0) %>%
#   ungroup() %>%
#   # add distance data
#   left_join(select(dat, subject_id, specimen_site, distance_meters, mean_distance_meters),
#             by = c("subject_id", "specimen_site")) %>%
#   distinct() %>%
#   identity() -> meta_asv_complete
# meta_asv_complete



#' calculate Shannon diversity and number of ASVs
meta_asv %>%
  group_by(specimen_id) %>%
  summarise(shannon = vegan::diversity(x = read_count),
            num_asvs = sum(read_count > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  identity() -> shannon
shannon



#' add qPCR data
read_csv("./tabs/iceman_qpcr.csv") %>%
  left_join(shannon, by = "specimen_id") %>%
  identity() -> qpcr
qpcr





#' ###########################################
#' 
#' qPCR
#' 
#' ###########################################

#' ###########################################
#' 
#' model total bacterial abundance as a result of GP(mean distance)
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
  # add qPCR bacterial abundance data
  left_join(qpcr, by = "specimen_id") %>%
  identity() -> meta_asv_genus_qpcr
meta_asv_genus_qpcr



meta_asv_genus_qpcr %>%
  # include Clostridium for Clostridioides difficile
  filter(genus %in% c(unique(iceman_corrected$genus), "Clostridium")) %>%
  #count(genus)
  left_join(select(dat, subject_id, specimen_site, distance_meters, mean_distance_meters),
            by = c("subject_id", "specimen_site")) %>%
  distinct() %>%
  identity() -> dat_asv_genus_qpcr
dat_asv_genus_qpcr



#' simple data form
dat_asv_genus_qpcr %>%
  select(specimen_id, subject_id, sampleevent, specimen_site, qpcr_copy_per_ul, shannon, num_asvs, distance_meters, mean_distance_meters) %>%
  distinct() %>%
  rename(qpcr = qpcr_copy_per_ul) %>%
  mutate_at(.vars = vars(qpcr, shannon, num_asvs), .funs = list("scaled" = ~ scale(.x)[,1], "log" = ~ log(.x))) %>%
  identity() -> dat_scaled
dat_scaled



#' get prior
dat_scaled %>%
  brms::get_prior(data = ., family = gaussian,
                  qpcr_scaled ~ gp(mean_distance_meters)
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
dat_scaled %>%
  brm(formula = qpcr_log ~ gp(mean_distance_meters),
      data = .,
      family = gaussian,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_qpcr

m_mean_distance_gp_qpcr %>% write_rds(file = "./models/gaussian/m_qpcr_v_mean_distance_gp.rds.gz", compress = "gz")
m_mean_distance_gp_qpcr <- read_rds(file = "./models/gaussian/m_qpcr_v_mean_distance_gp.rds.gz")

m_mean_distance_gp_qpcr
rstan::check_hmc_diagnostics(m_mean_distance_gp_qpcr$fit)
m_mean_distance_gp_qpcr %>% pp_check()

m_mean_distance_gp_qpcr %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_qpcr$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         ) %>%
  add_fitted_draws(m_mean_distance_gp_qpcr) %>%
  identity() -> m_mean_distance_gp_qpcr_fitted
m_mean_distance_gp_qpcr_fitted


#' fitted
m_mean_distance_gp_qpcr_fitted %>%
  mutate(.value = exp(.value)) %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  #scale_y_log10() +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "16S rRNA copies/μL",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.08, 0.78),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_qpcr_gp
p_qpcr_gp


# p_qpcr_gp %>%
#   ggsave(filename = "./figs/p_qpcr_gp.pdf", height = 4, width = 6, units = "in", device = cairo_pdf)
# p_qpcr_gp %>%
#   ggsave(filename = "./figs/p_qpcr_gp.png", height = 4, width = 6, units = "in", dpi = 600)
# p_qpcr_gp %>%
#   ggsave(filename = "./figs/p_qpcr_gp.svg", height = 4, width = 6, units = "in")





#' ###########################################
#' 
#' ADD RANDOM INTERCEPT FOR SUBJECTS & model as a result of GP(mean distance)
#' 
#' ###########################################

#' run model
dat_scaled %>%
  brm(formula = qpcr_log ~ gp(mean_distance_meters) + (1|subject_id),
      data = .,
      family = gaussian,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_qpcr_subject

m_mean_distance_gp_qpcr_subject %>% write_rds(file = "./models/gaussian/m_qpcr_v_mean_distance_gp_subject.rds.gz", compress = "gz")
m_mean_distance_gp_qpcr_subject <- read_rds(file = "./models/gaussian/m_qpcr_v_mean_distance_gp_subject.rds.gz")

m_mean_distance_gp_qpcr_subject
rstan::check_hmc_diagnostics(m_mean_distance_gp_qpcr_subject$fit)
m_mean_distance_gp_qpcr_subject %>% pp_check()

m_mean_distance_gp_qpcr_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_qpcr_subject$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         subject_id = unique(subject_id)
  ) %>%
  add_fitted_draws(m_mean_distance_gp_qpcr_subject) %>%
  identity() -> m_mean_distance_gp_qpcr_subject_fitted
m_mean_distance_gp_qpcr_subject_fitted


#' fitted
m_mean_distance_gp_qpcr_subject_fitted %>%
  mutate(.value = exp(.value)) %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  #scale_y_log10() +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "16S rRNA copies/μL",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        legend.direction = "horizontal",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_qpcr_gp_subject
p_qpcr_gp_subject


# p_qpcr_gp_subject %>%
#   ggsave(filename = "./figs/p_qpcr_gp_subject.pdf", height = 4, width = 6, units = "in", device = cairo_pdf)
# p_qpcr_gp_subject %>%
#   ggsave(filename = "./figs/p_qpcr_gp_subject.png", height = 4, width = 6, units = "in", dpi = 600)
# p_qpcr_gp_subject %>%
#   ggsave(filename = "./figs/p_qpcr_gp_subject.svg", height = 4, width = 6, units = "in")






#' ###########################################
#' 
#' Shannon Diversity
#' 
#' ###########################################

#' ###########################################
#' 
#' ADD RANDOM INTERCEPT FOR SUBJECTS & model as a result of GP(mean distance)
#' 
#' ###########################################

#' run model
dat_scaled %>%
  brm(formula = shannon_scaled ~ gp(mean_distance_meters) + (1|subject_id),
      data = .,
      family = gaussian,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.9999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_shannon_subject

m_mean_distance_gp_shannon_subject %>% write_rds(file = "./models/gaussian/m_shannon_v_mean_distance_gp_subject.rds.gz", compress = "gz")
m_mean_distance_gp_shannon_subject <- read_rds(file = "./models/gaussian/m_shannon_v_mean_distance_gp_subject.rds.gz")

m_mean_distance_gp_shannon_subject
rstan::check_hmc_diagnostics(m_mean_distance_gp_shannon_subject$fit)
m_mean_distance_gp_shannon_subject %>% pp_check()

m_mean_distance_gp_shannon_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_shannon_subject$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         subject_id = unique(subject_id)
  ) %>%
  add_fitted_draws(m_mean_distance_gp_shannon_subject) %>%
  # rescale scaled outcome
  mutate(.value = .value * sd(dat_scaled$shannon, na.rm = TRUE) + mean(dat_scaled$shannon, na.rm = TRUE)) %>%
  identity() -> m_mean_distance_gp_shannon_subject_fitted
m_mean_distance_gp_shannon_subject_fitted


#' fitted
m_mean_distance_gp_shannon_subject_fitted %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  #scale_y_log10() +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Shannon Diversity",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        legend.direction = "horizontal",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_shannon_gp_subject
p_shannon_gp_subject


# p_shannon_gp_subject %>%
#   ggsave(filename = "./figs/p_shannon_gp_subject.pdf", height = 4, width = 6, units = "in")
# p_shannon_gp_subject %>%
#   ggsave(filename = "./figs/p_shannon_gp_subject.png", height = 4, width = 6, units = "in", dpi = 600)
# p_shannon_gp_subject %>%
#   ggsave(filename = "./figs/p_shannon_gp_subject.svg", height = 4, width = 6, units = "in")








#' ###########################################
#' 
#' Richness (Number of ASVs)
#' 
#' ###########################################

#' ###########################################
#' 
#' ADD RANDOM INTERCEPT FOR SUBJECTS & model as a result of GP(mean distance)
#' 
#' ###########################################

#' run model
dat_scaled %>%
  brm(formula = num_asvs_scaled ~ gp(mean_distance_meters) + (1|subject_id),
      data = .,
      family = gaussian,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_mean_distance_gp_num_asvs_subject

m_mean_distance_gp_num_asvs_subject %>% write_rds(file = "./models/gaussian/m_num_asvs_v_mean_distance_gp_subject.rds.gz", compress = "gz")
m_mean_distance_gp_num_asvs_subject <- read_rds(file = "./models/gaussian/m_num_asvs_v_mean_distance_gp_subject.rds.gz")

m_mean_distance_gp_num_asvs_subject
rstan::check_hmc_diagnostics(m_mean_distance_gp_num_asvs_subject$fit)
m_mean_distance_gp_num_asvs_subject %>% pp_check()

m_mean_distance_gp_num_asvs_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_mean_distance_gp_num_asvs_subject$data %>%
  as_tibble() %>%
  expand(mean_distance_meters = modelr::seq_range(mean_distance_meters, n = 20),
         subject_id = unique(subject_id)
  ) %>%
  add_fitted_draws(m_mean_distance_gp_num_asvs_subject) %>%
  # rescale scaled outcome
  mutate(.value = .value * sd(dat_scaled$num_asvs, na.rm = TRUE) + mean(dat_scaled$num_asvs, na.rm = TRUE)) %>%
  identity() -> m_mean_distance_gp_num_asvs_subject_fitted
m_mean_distance_gp_num_asvs_subject_fitted


#' fitted
m_mean_distance_gp_num_asvs_subject_fitted %>%
  ggplot(aes(x = mean_distance_meters, y = .value)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  #scale_y_log10() +
  labs(x = "Distance from Patient & Towards Wastewater Sites (meters)",
       y = "Number of ASVs",
       fill = "Posterior Credible Interval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = "top",
        legend.direction = "horizontal",
        #legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_num_asvs_gp_subject
p_num_asvs_gp_subject


# p_num_asvs_gp_subject %>%
#   ggsave(filename = "./figs/p_num_asvs_gp_subject.pdf", height = 4, width = 6, units = "in")
# p_num_asvs_gp_subject %>%
#   ggsave(filename = "./figs/p_num_asvs_gp_subject.png", height = 4, width = 6, units = "in", dpi = 600)
# p_num_asvs_gp_subject %>%
#   ggsave(filename = "./figs/p_num_asvs_gp_subject.svg", height = 4, width = 6, units = "in")









#' ###########################################
#' 
#' COMBINE FIGURES
#' 
#' ###########################################

library(patchwork)

p_qpcr_gp_subject +
  p_shannon_gp_subject + scale_y_continuous(limits = c(0,NA)) +
  p_num_asvs_gp_subject + scale_y_continuous(limits = c(0,NA)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 3, widths = c(1,1,1), guides = "collect") &
  #scale_fill_brewer(palette = "Purples") &
  theme(plot.margin = margin(5,5,5,5), plot.tag.position = c(0.05,0.85), legend.position = "top", legend.direction = "horizontal", axis.title.x = element_blank()) %>%
  identity() -> p
p

patchwork::patchworkGrob(p & theme(axis.title.x = element_blank())) -> p_gt

gridExtra::grid.arrange(p_gt, bottom = "Distance from Patient & Towards Wastewater Sites (meters)") %>%
  identity() -> p_gt_extra
p_gt_extra


# p_gt_extra %>%
#   ggsave(filename = "./figs/p_microbiome_gp_combined.pdf", height = 4, width = 8, units = "in", device = cairo_pdf)
# p_gt_extra %>%
#   ggsave(filename = "./figs/p_microbiome_gp_combined.png", height = 4, width = 8, units = "in", dpi = 600)
# p_gt_extra %>%
#   ggsave(filename = "./figs/p_microbiome_gp_combined.svg", height = 4, width = 8, units = "in")









#' ###########################################
#' 
#' qPCR, DIVERSITY, and RICHNESS vs CULTURE
#' 
#' ###########################################

qpcr %>%
  left_join(select(meta_asv, specimen_id, subject_id, specimen_site, sampleevent)) %>%
  distinct() %>%
  left_join(dat, by = c("subject_id", "specimen_site", "sampleevent")) %>%
  rename(qpcr = qpcr_copy_per_ul) %>%
  mutate_at(.vars = vars(qpcr, shannon, num_asvs), .funs = list("scaled" = ~ scale(.x)[,1], "log" = ~ log(.x))) %>%
  identity() -> dat_cx_qpcr
dat_cx_qpcr


#' run model: qPCR
dat_cx_qpcr %>%
  brm(formula = cx_positive ~ (qpcr_log + 0 | genus) + (1 | subject_id),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_cx_qpcr_mix_genus_subject

m_cx_qpcr_mix_genus_subject %>% write_rds(file = "./models/binomial/m_cx_v_qpcr_mix_genus_subject.rds.gz", compress = "gz")
m_cx_qpcr_mix_genus_subject <- read_rds(file = "./models/binomial/m_cx_v_qpcr_mix_genus_subject.rds.gz")

m_cx_qpcr_mix_genus_subject
rstan::check_hmc_diagnostics(m_cx_qpcr_mix_genus_subject$fit)
m_cx_qpcr_mix_genus_subject %>% pp_check()

m_cx_qpcr_mix_genus_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_cx_qpcr_mix_genus_subject$data %>%
  as_tibble() %>%
  expand(qpcr_log = modelr::seq_range(qpcr_log, n = 20),
         genus = unique(genus),
         subject_id = unique(subject_id)) %>%
  add_fitted_draws(m_cx_qpcr_mix_genus_subject) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  mutate(qpcr = exp(qpcr_log)) %>%
  identity() -> m_cx_qpcr_mix_genus_subject_fitted
m_cx_qpcr_mix_genus_subject_fitted


m_cx_qpcr_mix_genus_subject_fitted %>%
  ggplot(aes(x = qpcr_log, y = .value)) +
  facet_wrap(facets = ~ genus_label, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Total Bacterial Abundance by 16S rRNA Gene qPCR (log copies per μL)",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_qpcr
p_cx_qpcr




#' run model: shannon
dat_cx_qpcr %>%
  brm(formula = cx_positive ~ (shannon_scaled + 0 | genus) + (1 | subject_id),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_cx_shannon_mix_genus_subject

m_cx_shannon_mix_genus_subject %>% write_rds(file = "./models/binomial/m_cx_v_shannon_mix_genus_subject.rds.gz", compress = "gz")
m_cx_shannon_mix_genus_subject <- read_rds(file = "./models/binomial/m_cx_v_shannon_mix_genus_subject.rds.gz")

m_cx_shannon_mix_genus_subject
rstan::check_hmc_diagnostics(m_cx_shannon_mix_genus_subject$fit)
m_cx_shannon_mix_genus_subject %>% pp_check()

m_cx_shannon_mix_genus_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_cx_shannon_mix_genus_subject$data %>%
  as_tibble() %>%
  expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 20),
         genus = unique(genus),
         subject_id = unique(subject_id)) %>%
  mutate(shannon = shannon_scaled * sd(dat_cx_qpcr$shannon, na.rm = TRUE) + mean(dat_cx_qpcr$shannon, na.rm = TRUE)) %>%
  add_fitted_draws(m_cx_shannon_mix_genus_subject) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_cx_shannon_mix_genus_subject_fitted
m_cx_shannon_mix_genus_subject_fitted


m_cx_shannon_mix_genus_subject_fitted %>%
  ggplot(aes(x = shannon, y = .value)) +
  facet_wrap(facets = ~ genus_label, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Bacterial Community Shannon Diversity",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_shannon
p_cx_shannon






#' run model: num_asvs
dat_cx_qpcr %>%
  brm(formula = cx_positive ~ (num_asvs_scaled + 0 | genus) + (1 | subject_id),
      data = .,
      family = bernoulli,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 18),
      backend = "cmdstanr",
      seed = 16) -> m_cx_num_asvs_mix_genus_subject

m_cx_num_asvs_mix_genus_subject %>% write_rds(file = "./models/binomial/m_cx_v_num_asvs_mix_genus_subject.rds.gz", compress = "gz")
m_cx_num_asvs_mix_genus_subject <- read_rds(file = "./models/binomial/m_cx_v_num_asvs_mix_genus_subject.rds.gz")

m_cx_num_asvs_mix_genus_subject
rstan::check_hmc_diagnostics(m_cx_num_asvs_mix_genus_subject$fit)
m_cx_num_asvs_mix_genus_subject %>% pp_check()

m_cx_num_asvs_mix_genus_subject %>%
  posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 2:5, n_sigfig = 3)


#' fitted
m_cx_num_asvs_mix_genus_subject$data %>%
  as_tibble() %>%
  expand(num_asvs_scaled = modelr::seq_range(num_asvs_scaled, n = 20),
         genus = unique(genus),
         subject_id = unique(subject_id)) %>%
  mutate(num_asvs = num_asvs_scaled * sd(dat_cx_qpcr$num_asvs, na.rm = TRUE) + mean(dat_cx_qpcr$num_asvs, na.rm = TRUE)) %>%
  add_fitted_draws(m_cx_num_asvs_mix_genus_subject) %>%
  left_join(distinct(select(dat, genus, genus_label)), by = "genus") %>%
  identity() -> m_cx_num_asvs_mix_genus_subject_fitted
m_cx_num_asvs_mix_genus_subject_fitted


m_cx_num_asvs_mix_genus_subject_fitted %>%
  ggplot(aes(x = num_asvs, y = .value)) +
  facet_wrap(facets = ~ genus_label, scales = "free_y") +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Bacterial Community Richness (Number of ASVs)",
       y = "Probability of Positive MDRO Culture",
       fill = "Posterior\nCredible\nInterval") +
  theme_bw() +
  theme(strip.text = ggtext::element_markdown(color = "black", size = 8),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        legend.position = c(0.9, 0.12),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        strip.background = element_blank()) -> p_cx_num_asvs
p_cx_num_asvs




