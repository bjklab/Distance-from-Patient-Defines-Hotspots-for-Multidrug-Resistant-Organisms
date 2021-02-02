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



#' ###########################################
#' quantify missing data (from COVID-lost plates)
#' ###########################################

iceman_corrected %>%
  count(subject_id, specimen_site, sampleevent) %>%
  filter(n < 11) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")





#' ###########################################
#' culture data
#' ###########################################

#' media types
# haven::read_dta("./data/10_16_20 ICEMAN LabVantage File Brendan.dta") %>%
#   count(Media)


#' subject and specimen tallies
iceman_corrected %>%
  group_by(specimen_site) %>%
  summarise(subjects_sampled = n_distinct(subject_id),
            specimens_collected = length(unique(paste0(subject_id, "_", sampleevent)))) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")
  


#' subject sampling event tallies
iceman_corrected %>%
  group_by(subject_id) %>%
  summarise(sampling_events = n_distinct(sampleevent)) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")


#' subject sampling event summaries
iceman_corrected %>%
  group_by(subject_id) %>%
  summarise(sampling_events = n_distinct(sampleevent)) %>%
  ungroup() %>%
  summarise(`Median` = median(sampling_events, na.rm = TRUE),
            `IQR` = IQR(sampling_events, na.rm = TRUE),
            `Mean` = mean(sampling_events, na.rm = TRUE),
            `SD` = sd(sampling_events, na.rm = TRUE),
            Sum = sum(sampling_events, na.rm = TRUE)) %>%
  gt::gt() %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")




#' genera detection tallies
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
  group_by(genus_label) %>%
  summarise(positive_cultures = sum(cx_positive, na.rm = TRUE),
            percent_positive = positive_cultures / n(),
            ) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::cols_label(genus_label = "MDRO", positive_cultures = gt::html("Total Positive<br>Cultures"), percent_positive = "Percent Positive") %>%
  gt::fmt_percent(columns = 3) %>%
  gt::fmt_markdown(columns = 1) %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black") -> tab1

tab1

tab1 %>%
  gt::as_raw_html() %>%
  write_lines(file = "./tabs/tab1_mdro_occurrence.html")




#' ############################################
#' ASV data
#' ############################################

meta_asv %>%
  group_by(specimen_site) %>%
  summarise(n_subject = n_distinct(subject_id),
            n_specimens = length(unique(specimen_id))) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_options(table_body.hlines.color = "black",
                  column_labels.border.bottom.color = "black",
                  column_labels.border.top.color = "black",
                  table_body.border.bottom.color = "black")



