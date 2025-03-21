#####################   Meta-Analysis MetaDRRap #####################
####### load packages + source the function script #######
source("Code/missing_sd_funcs.R")
source("Code/function_orchardsummary.R")
library(tidyverse)
library(metafor)
library(orchaRd)
library(patchwork)
library(janitor)
library(gt)
library(emmeans)

sink("sessionInfo.txt")
sessionInfo()
sink()
##write_delim()####### read in data and create cv for everything if it exists #######
dat <- read_csv("Data/analysis_data.csv") %>%
  mutate(
    cv_control = na_if(sd_control / m_control, Inf),
    cv_experimental = na_if(sd_treatment / m_treatment, Inf)
  ) %>%
  janitor::clean_names() %>%
  select(-c(dose_compared_to_al, age_on_treatment, treatment_length)) %>%
  mutate(
    treatment_type = as.factor(treatment_type),
    treatment_overall = as.factor(treatment_overall),
    species = as.factor(species),
    sex = as.factor(sex),
    id = as.factor(id),
    title = as.factor(title)
  )

####### summarise overall number of titles,measures,species,sex,treatments #######
length(unique(dat$title))

dat %>%
  group_by(m_measure) %>%
  summarise(n = n(),
            n_stud = n_distinct(title))

dat %>%
  group_by(species) %>%
  summarise(n = n(),
            n_stud = n_distinct(title))

dat %>%
  group_by(sex) %>%
  summarise(n = n(),
            n_stud = n_distinct(title))

dat %>%
  group_by(treatment_overall) %>%
  summarise(n = n(),
            n_stud = n_distinct(title))

dat %>%
  group_by(treatment_type) %>%
  summarise(n = n(),
            n_stud = n_distinct(title))

####### geary function to test for outliers #######
geary <- function(mean, sd, n) {
  (1 / (sd / mean)) * ((4 * n^(3 / 2)) / (1 + 4 * n))
}

####### number of outliers pass/na #######
dat <- dat %>%
  mutate(
    geary_control = geary(m_control, sd_control, n_control),
    geary_trt = geary(m_treatment, sd_treatment, n_treatment),
    geary_test = ifelse(geary_control >= 3 & geary_trt >= 3, "pass", "fail")
  ) # %>%
# filter(is.na(geary_test) | geary_test == "pass")

####### How many fail? #######
geary_res <- dat %>%
  group_by(geary_test) %>%
  summarise(n = n()) %>%
  data.frame()


####### Calculate the average between study CV, which will replace missing values #######
# 2 = control, 1 = treatment
dat <- cv_avg(
  x = m_control, sd = sd_control,
  n = n_control, group = title, label = "2",
  data = dat
)

dat <- cv_avg(
  x = m_treatment, sd = sd_treatment,
  n = n_treatment, group = title,
  label = "1", data = dat
)

####### replace missing values #######
dat <- dat %>%
  mutate(
    cv2_cont_new = if_else(is.na(cv_control), b_CV2_2, cv_control^2),
    cv2_expt_new = if_else(is.na(cv_experimental), b_CV2_1, cv_experimental^2)
  )

##################### all_Cases + lajanueese adjustment for lnRR and variance #####################
dat <- dat %>%
  mutate(
    lnrr_laj = lnrr_laj(
      m1 = m_treatment, m2 = m_control,
      cv1_2 = cv2_expt_new, cv2_2 = cv2_cont_new,
      n1 = n_treatment, n2 = n_control
    ),
    v_lnrr_laj_1B = v_lnrr_laj(
      cv1 = b_CV2_1, n1 = n_treatment,
      cv2 = b_CV2_2, n2 = n_control
    )
  )

####### create the new variables for pub bias correction #######
dat$inv_n_tilda <- with(dat, (n_control + n_treatment) / (n_control * n_treatment))
dat$sqrt_inv_n_tilda <- with(dat, sqrt(inv_n_tilda))
dat$year.c <- as.vector(scale(dat$year, scale = F))

##################### treatment effects #####################
####### model with mean + summary #######
method_AC_meta_mean <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "mean")
)

summary(method_AC_meta_mean)
summary_mean <- orchaRd::mod_results(method_AC_meta_mean,
  mod = "treatment_overall",
  group = "title"
)

mean_treatment <- orchaRd_table(
  model = method_AC_meta_mean,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin"
  )
)

png("mean_treatment.png", res = 600, width = 8, height = 10, units = "in")
print(mean_treatment)
dev.off()

####### model with pub bias mean + summary #######
method_AC_meta_mean_pub <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "mean")
)

summary(method_AC_meta_mean_pub)

summary_mean_pub <- orchaRd::mod_results(method_AC_meta_mean_pub,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

mean_treatment_pub <- orchaRd_table(
  model = method_AC_meta_mean_pub,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin", "1/n", "Year"
  ),
  pb = c("inv_n_tilda", "year.c")
)

png("mean_treatment_pub.png", res = 600, width = 8, height = 10, units = "in")
print(mean_treatment_pub)
dev.off()
####### model with median + summary #######
method_AC_meta_median <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "median")
)

summary(method_AC_meta_median)
summary_median <- orchaRd::mod_results(method_AC_meta_median,
  mod = "treatment_overall",
  group = "title"
)

median_treatment <- orchaRd_table(
  model = method_AC_meta_median,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin"
  )
)

png("median_treatment.png", res = 600, width = 8, height = 10, units = "in")
print(median_treatment)
dev.off()

####### model with pub bias median + summary #######
method_AC_meta_median_pub <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "median")
)

summary(method_AC_meta_median_pub)

summary_median_pub <- orchaRd::mod_results(method_AC_meta_median_pub,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

median_treatment_pub <- orchaRd_table(
  model = method_AC_meta_median_pub,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin", "1/n", "Year"
  ),
  pb = c("inv_n_tilda", "year.c")
)

png("median_treatment_pub.png", res = 600, width = 8, height = 10, units = "in")
print(median_treatment_pub)
dev.off()

####### model with all + summary #######
method_AC_meta_all <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_AC_meta_all)

summary_all <- orchaRd::mod_results(method_AC_meta_all,
  mod = "treatment_overall",
  group = "title"
)

all_treatment <- orchaRd_table(
  model = method_AC_meta_all,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin"
  )
)

png("all_treatment.png", res = 600, width = 8, height = 10, units = "in")
print(all_treatment)
dev.off()

####### model with pub bias all + summary #######
method_AC_meta_all_pub <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_AC_meta_all_pub)

summary_all_pub <- orchaRd::mod_results(method_AC_meta_all_pub,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

all_treatment_pub <- orchaRd_table(
  model = method_AC_meta_all_pub,
  mod = "treatment_overall",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Metformin", "Rapamycin", "1/n", "Year"
  ),
  pb = c("inv_n_tilda", "year.c")
)

png("all_treatment_pub.png", res = 600, width = 8, height = 10, units = "in")
print(all_treatment_pub)
dev.off()


####### graph all effects #######

all_pub_model <- summary_all_pub$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures")

median_pub_model <- summary_median_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Medians")

mean_pub_model <- summary_mean_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Means")

all_model <- summary_all$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures")

median_model <- summary_median$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians")

mean_model <- summary_mean$mod_table %>%
  mutate(model = "Unadjusted", data = "Means")

all_mods <- bind_rows(
  all_pub_model,
  median_pub_model,
  mean_pub_model,
  all_model,
  median_model,
  mean_model
)

ggplot() +
  geom_point(
    data = all_mods, aes(
      x = estimate,
      y = name,
      colour = model,
      shape = data
    ),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbarh(
    data = all_mods, aes(
      x = estimate,
      y = name,
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data,
      linetype = name
    ),
    height = 0.4,
    linewidth = 1.1,
    position = position_dodge(width = 0.6),
    show.legend = F
  ) +
  theme_bw(base_size = 13) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Model Type",
    shape = "Measure Type",
    x = "Log Response Ratio",
    y = "Treatment",
  ) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") -> overall

ggsave("Overall.png", overall, dpi = 600, width = 20, height = 8)


####### i2 for full model w/ table #######
method_AC_meta_all <- rma.mv(lnrr_laj ~ 1,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_AC_meta_all)

i2.model.smd <- i2_ml(method_AC_meta_all)

table.model.smd <- data.frame(
  n = length(unique(dat$title)),
  k = nrow(dat),
  I2_obsID = round(i2.model.smd[[4]], 1),
  I2_paperID = round(i2.model.smd[[3]], 1),
  I2_species = round(i2.model.smd[[2]], 1),
  I2_total = round(i2.model.smd[[1]], 1)
)

rownames(table.model.smd) <- NULL

table.model.smd.gt <- table.model.smd %>%
  gt::gt() %>%
  cols_label(
    n = md("**n**"),
    k = md("**k**"),
    I2_obsID = md("***I*<sup>2</sup><sub>residual</sub> (%)**"),
    I2_paperID = md("***I*<sup>2</sup><sub>study</sub> (%)**"),
    I2_species = md("***I*<sup>2</sup><sub>species</sub> (%)**"),
    I2_total = md("***I*<sup>2</sup><sub>total</sub> (%)**")
  ) %>%
  cols_align(align = "center")

table.model.smd.gt

####### relevel for metformin comparison w/ or wo/ pub bias correction #######
dat$treatment_overall <- relevel(dat$treatment_overall, ref = "Metformin")

method_AC_meta_all <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_AC_meta_all)

method_AC_meta_all_pub <- rma.mv(
  lnrr_laj ~treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_AC_meta_all_pub)

####### relevel back to DR #######
dat$treatment_overall <- relevel(dat$treatment_overall, ref = "Dietary restriction")


##################### SEX effects #####################
####### Rap Mean #######
sex_AC_Rap_mean <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_mean)
summary_sex_AC_Rap_mean <- orchaRd::mod_results(sex_AC_Rap_mean,
  mod = "sex",
  group = "title"
)

mean_sex_rap <- orchaRd_table(
  model = sex_AC_Rap_mean,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("mean_sex_rap.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_rap)
dev.off()

####### Rap Pub Bias Mean #######
sex_AC_Rap_mean_pub <- rma.mv(
  lnrr_laj ~ -1+sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_mean_pub)

summary_sex_AC_Rap_mean_pub <- orchaRd::mod_results(sex_AC_Rap_mean_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

mean_sex_rap_pub <- orchaRd_table(
  model = sex_AC_Rap_mean_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("mean_sex_rap_pub.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_rap_pub)
dev.off()

####### Rap Median #######
sex_AC_Rap_median <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_median)
summary_sex_AC_Rap_median <- orchaRd::mod_results(sex_AC_Rap_median,
  mod = "sex",
  group = "title"
)

median_sex_rap <- orchaRd_table(
  model = sex_AC_Rap_median,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("median_sex_rap.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_rap)
dev.off()


####### Rap Pub Bias Median #######
sex_AC_Rap_median_pub <- rma.mv(
  lnrr_laj ~ -1+sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_median_pub)

summary_sex_AC_Rap_median_pub <- orchaRd::mod_results(sex_AC_Rap_median_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

median_sex_rap_pub <- orchaRd_table(
  model = sex_AC_Rap_median_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("median_sex_rap_pub.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_rap_pub)
dev.off()

####### Rap All #######
sex_AC_Rap_all <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_all)
summary_sex_AC_Rap_all <- orchaRd::mod_results(sex_AC_Rap_all,
  mod = "sex",
  group = "title"
)

all_sex_rap <- orchaRd_table(
  model = sex_AC_Rap_all,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("all_sex_rap.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_rap)
dev.off()

####### Rap Pub Bias All #######
sex_AC_Rap_all_pub <- rma.mv(
  lnrr_laj ~ -1+sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Rapamycin"
  )
)

summary(sex_AC_Rap_all_pub)

summary_sex_AC_Rap_all_pub <- orchaRd::mod_results(sex_AC_Rap_all_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

all_sex_rap_pub <- orchaRd_table(
  model = sex_AC_Rap_all_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("all_sex_rap_pub.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_rap_pub)
dev.off()

####### Met Mean #######
sex_AC_Met_mean <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Metformin"
  )
)

summary(sex_AC_Met_mean)
summary_sex_AC_Met_mean <- orchaRd::mod_results(sex_AC_Met_mean,
  mod = "sex",
  group = "title"
)

mean_sex_met <- orchaRd_table(
  model = sex_AC_Met_mean,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("mean_sex_met.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_met)
dev.off()

####### Met Pub Bias Mean #######
sex_AC_Met_mean_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Metformin"
  )
)

summary(sex_AC_Met_mean_pub)

summary_sex_AC_Met_mean_pub <- orchaRd::mod_results(sex_AC_Met_mean_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

mean_sex_met_pub <- orchaRd_table(
  model = sex_AC_Met_mean_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("mean_sex_met_pub.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_met_pub)
dev.off()

####### Met Median #######
sex_AC_Met_median <- rma.mv(lnrr_laj ~ sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Metformin"
  ),
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
)

summary(sex_AC_Met_median)
summary_sex_AC_Met_median <- orchaRd::mod_results(sex_AC_Met_median,
  mod = "sex",
  group = "title"
)

median_sex_met <- orchaRd_table(
  model = sex_AC_Met_median,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("median_sex_met.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_met)
dev.off()

####### Met Pub Bias Median #######
sex_AC_Met_median_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Metformin"
  )
)

summary(sex_AC_Met_median_pub)

summary_sex_AC_Met_median_pub <- orchaRd::mod_results(sex_AC_Met_median_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

median_sex_met_pub <- orchaRd_table(
  model = sex_AC_Met_median_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("median_sex_met_pub.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_met_pub)
dev.off()

####### Met All #######
sex_AC_Met_all <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Metformin"
  )
)

summary(sex_AC_Met_all)
summary_sex_AC_Met_all <- orchaRd::mod_results(sex_AC_Met_all,
  mod = "sex",
  group = "title"
)

all_sex_met <- orchaRd_table(
  model = sex_AC_Met_all,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("all_sex_met.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_met)
dev.off()

####### Met Pub Bias All #######
sex_AC_Met_all_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Metformin"
  )
)

summary(sex_AC_Met_all_pub)

summary_sex_AC_Met_all_pub <- orchaRd::mod_results(sex_AC_Met_all_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

all_sex_met_pub <- orchaRd_table(
  model = sex_AC_Met_all_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("all_sex_met_pub.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_met_pub)
dev.off()

####### DR Mean #######
sex_AC_DR_mean <- rma.mv(lnrr_laj ~ -1+sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Dietary restriction"
  )
)

summary(sex_AC_DR_mean)
summary_sex_AC_DR_mean <- orchaRd::mod_results(sex_AC_DR_mean,
  mod = "sex",
  group = "title"
)

mean_sex_dr <- orchaRd_table(
  model = sex_AC_DR_mean,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("mean_sex_dr.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_dr)
dev.off()

####### DR Pub Bias Mean #######
sex_AC_DR_mean_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Dietary restriction"
  ),
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
)

summary(sex_AC_DR_mean_pub)

summary_sex_AC_DR_mean_pub <- orchaRd::mod_results(sex_AC_DR_mean_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

mean_sex_dr_pub <- orchaRd_table(
  model = sex_AC_DR_mean_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("mean_sex_dr_pub.png", res = 600, width = 8, height = 10, units = "in")
print(mean_sex_dr_pub)
dev.off()

####### DR Median #######
sex_AC_DR_median <- rma.mv(lnrr_laj ~ sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Dietary restriction"
  )
)

summary(sex_AC_DR_median)
summary_sex_AC_DR_median <- orchaRd::mod_results(sex_AC_DR_median,
  mod = "sex",
  group = "title"
)

median_sex_dr <- orchaRd_table(
  model = sex_AC_DR_median,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("median_sex_dr.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_dr)
dev.off()

####### DR Pub Bias Median######
sex_AC_DR_median_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Dietary restriction"
  )
)

summary(sex_AC_DR_median_pub)

summary_sex_AC_DR_median_pub <- orchaRd::mod_results(sex_AC_DR_median_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

median_sex_dr_pub <- orchaRd_table(
  model = sex_AC_DR_median_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("median_sex_dr_pub.png", res = 600, width = 8, height = 10, units = "in")
print(median_sex_dr_pub)
dev.off()

####### DR All #######
sex_AC_DR_all <- rma.mv(lnrr_laj ~ sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Dietary restriction"
  )
)

summary(sex_AC_DR_all)
summary_sex_AC_DR_all <- orchaRd::mod_results(sex_AC_DR_all,
  mod = "sex",
  group = "title"
)

all_sex_dr <- orchaRd_table(
  model = sex_AC_DR_all,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed")
)

png("all_sex_dr.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_dr)
dev.off()

####### DR Pub Bias All #######
sex_AC_DR_all_pub <- rma.mv(
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    treatment_overall == "Dietary restriction"
  )
)

summary(sex_AC_DR_all_pub)

summary_sex_AC_DR_all_pub <- orchaRd::mod_results(sex_AC_DR_all_pub,
  mod = "sex",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

all_sex_dr_pub <- orchaRd_table(
  model = sex_AC_DR_all_pub,
  mod = "sex",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Male", "Mixed", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("all_sex_dr_pub.png", res = 600, width = 8, height = 10, units = "in")
print(all_sex_dr_pub)
dev.off()

####### graph all effects#######

mean_pub_model_rap <- summary_sex_AC_Rap_mean_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Means", treatment = "Rapamycin")

mean_pub_model_met <- summary_sex_AC_Met_mean_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Means", treatment = "Metformin")

mean_pub_model_dr <- summary_sex_AC_DR_mean_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Means", treatment = "Dietary restriction")

mean_model_rap <- summary_sex_AC_Rap_mean$mod_table %>%
  mutate(model = "Unadjusted", data = "Means", treatment = "Rapamycin")

mean_model_met <- summary_sex_AC_Met_mean$mod_table %>%
  mutate(model = "Unadjusted", data = "Means", treatment = "Metformin")

mean_model_dr <- summary_sex_AC_DR_mean$mod_table %>%
  mutate(model = "Unadjusted", data = "Means", treatment = "Dietary restriction")

median_pub_model_rap <- summary_sex_AC_Rap_median_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Medians", treatment = "Rapamycin")

median_pub_model_met <- summary_sex_AC_Met_median_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Medians", treatment = "Metformin")

median_pub_model_dr <- summary_sex_AC_DR_median_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Medians", treatment = "Dietary restriction")

median_model_rap <- summary_sex_AC_Rap_median$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians", treatment = "Rapamycin")

median_model_met <- summary_sex_AC_Met_median$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians", treatment = "Metformin")

median_model_dr <- summary_sex_AC_DR_median$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians", treatment = "Dietary restriction")

all_pub_model_rap <- summary_sex_AC_Rap_all_pub$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures", treatment = "Rapamycin")

all_pub_model_met <- summary_sex_AC_Met_all_pub$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures", treatment = "Metformin")

all_pub_model_dr <- summary_sex_AC_DR_all_pub$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures", treatment = "Dietary restriction")

all_model_rap <- summary_sex_AC_Rap_all$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures", treatment = "Rapamycin")

all_model_met <- summary_sex_AC_Met_all$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures", treatment = "Metformin")

all_model_dr <- summary_sex_AC_DR_all$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures", treatment = "Dietary restriction")

# bind
all_mods <- bind_rows(
  all_pub_model_rap,
  all_pub_model_met,
  all_pub_model_dr,
  all_model_rap,
  all_model_met,
  all_model_dr,
  median_pub_model_rap,
  median_pub_model_met,
  median_pub_model_dr,
  median_model_rap,
  median_model_met,
  median_model_dr,
  mean_pub_model_rap,
  mean_pub_model_met,
  mean_pub_model_dr,
  mean_model_rap,
  mean_model_met,
  mean_model_dr
)

ggplot() +
  geom_point(
    data = all_mods, aes(
      x = estimate,
      y = name,
      colour = model,
      shape = data
    ),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbarh(
    data = all_mods, aes(
      x = estimate,
      y = name,
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data,
      linetype = name
    ),
    height = 0.4,
    linewidth = 1.1,
    position = position_dodge(width = 0.6),
    show.legend = F
  ) +
  theme_bw(base_size = 13) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Model Type",
    shape = "Measure Type",
    x = "Log Response Ratio",
    y = "Sex",
  ) +
  theme(legend.position = "bottom") +
  ggthemes::scale_colour_colorblind() +
  facet_grid(. ~ treatment) -> sex_plot

ggsave("Sex.png", sex_plot, dpi = 600, width = 20, height = 8)


##################### DR method #####################
####### DR Mean #######
type_method_AC_DR_mean <- rma.mv(lnrr_laj ~ treatment_type,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Dietary restriction"
  )
)

summary(type_method_AC_DR_mean)
summary_type_method_AC_DR_mean <- orchaRd::mod_results(type_method_AC_DR_mean,
  mod = "treatment_type",
  group = "title"
)

mean_type_dr <- orchaRd_table(
  model = type_method_AC_DR_mean,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Percent Reduction")
)

png("mean_type_dr.png", res = 600, width = 9, height = 10, units = "in")
print(mean_type_dr)
dev.off()

####### DR Pub Bias Mean #######
type_method_AC_DR_mean_pub <- rma.mv(
  lnrr_laj ~ treatment_type + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    treatment_overall == "Dietary restriction"
  ),
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
)

summary(type_method_AC_DR_mean_pub)

summary_type_method_AC_DR_mean_pub <- orchaRd::mod_results(type_method_AC_DR_mean_pub,
  mod = "treatment_type",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)


mean_type_dr_pub <- orchaRd_table(
  model = type_method_AC_DR_mean_pub,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Percent Reduction", "1/n", "Year"),
  pb = c("inv_n_tilda", "year.c")
)

png("mean_type_dr_pub.png", res = 600, width = 9, height = 10, units = "in")
print(mean_type_dr_pub)
dev.off()

####### DR Median #######
type_method_AC_DR_median <- rma.mv(lnrr_laj ~ treatment_type,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Dietary restriction"
  )
)

summary(type_method_AC_DR_median)

summary_type_method_AC_DR_median <- orchaRd::mod_results(type_method_AC_DR_median,
  mod = "treatment_type",
  group = "title"
)

median_type_dr <- orchaRd_table(
  model = type_method_AC_DR_median,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Percent Reduction", "Percent Reduction + Fasted")
)

png("median_type_dr.png", res = 600, width = 9, height = 10, units = "in")
print(median_type_dr)
dev.off()

####### DR Pub Bias Median #######
type_method_AC_DR_median_pub <- rma.mv(
  lnrr_laj ~ treatment_type + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    treatment_overall == "Dietary restriction"
  )
)

summary(type_method_AC_DR_median_pub)

summary_type_method_AC_DR_median_pub <- orchaRd::mod_results(type_method_AC_DR_median_pub,
  mod = "treatment_type",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

median_type_dr_pub <- orchaRd_table(
  model = type_method_AC_DR_median_pub,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Percent Reduction", "Percent Reduction + Fasted",
    "1/n", "Year"
  ),
  pb = c("inv_n_tilda", "year.c")
)

png("median_type_dr_pub.png", res = 600, width = 9, height = 10, units = "in")
print(median_type_dr_pub)
dev.off()

####### DR All #######
type_method_AC_DR_all <- rma.mv(lnrr_laj ~ treatment_type,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(treatment_overall == "Dietary restriction")
)

summary(type_method_AC_DR_all)

summary_type_method_AC_DR_all <- orchaRd::mod_results(type_method_AC_DR_all,
  mod = "treatment_type", group = "title"
)

all_type_dr <- orchaRd_table(
  model = type_method_AC_DR_all,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c("Intercept", "Percent Reduction", "Percent Reduction + Fasted")
)

png("all_type_dr.png", res = 600, width = 9, height = 10, units = "in")
print(all_type_dr)
dev.off()

####### DR Pub Bias All #######
type_method_AC_DR_all_pub <- rma.mv(
  lnrr_laj ~ treatment_type + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(treatment_overall == "Dietary restriction")
)

summary(type_method_AC_DR_all_pub)

summary_type_method_AC_DR_all_pub <- orchaRd::mod_results(type_method_AC_DR_all_pub,
  mod = "treatment_type",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

all_type_dr_pub <- orchaRd_table(
  model = type_method_AC_DR_all_pub,
  mod = "treatment_type",
  es = "logRR",
  group = "title",
  terms = c(
    "Intercept", "Percent Reduction", "Percent Reduction + Fasted",
    "1/n", "Year"
  ),
  pb = c("inv_n_tilda", "year.c")
)

png("all_type_dr_pub.png", res = 600, width = 9, height = 10, units = "in")
print(all_type_dr_pub)
dev.off()

####### graph all effects#####

all_pub_model <- summary_type_method_AC_DR_all_pub$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures")

median_pub_model <- summary_type_method_AC_DR_median_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Medians")

mean_pub_model <- summary_type_method_AC_DR_mean_pub$mod_table %>%
  mutate(model = "Adjusted", data = "Means")

all_model <- summary_type_method_AC_DR_all$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures")

median_model <- summary_type_method_AC_DR_median$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians")

mean_model <- summary_type_method_AC_DR_mean$mod_table %>%
  mutate(model = "Unadjusted", data = "Means")

all_mods <- bind_rows(
  all_pub_model,
  median_pub_model,
  mean_pub_model,
  all_model,
  median_model,
  mean_model
)

ggplot() +
  geom_point(
    data = all_mods, aes(
      x = estimate,
      y = name,
      colour = model,
      shape = data
    ),
    size = 4,
    position = position_dodge(width = 0.6)
  ) +
  geom_errorbarh(
    data = all_mods, aes(
      x = estimate,
      y = name,
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data,
      linetype = name
    ),
    height = 0.4,
    linewidth = 1.1,
    position = position_dodge(width = 0.6),
    show.legend = F
  ) +
  theme_bw(base_size = 13) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Model Type",
    shape = "Measure Type",
    x = "Log Response Ratio",
    y = "Restirction Method",
  ) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set2") -> DR_type

ggsave("DR_type.png", DR_type, dpi = 600, width = 20, height = 8)

##################### missing cases+ lajanueese adjustment for lnRR and variance #####################
####### create new data with different method #######
dat <- dat %>%
  mutate(
    lnrr_laj = lnrr_laj(
      m1 = m_treatment, m2 = m_control,
      cv1_2 = cv2_expt_new, cv2_2 = cv2_cont_new,
      n1 = n_treatment, n2 = n_control
    ),
    v_lnrr_laj = v_lnrr_laj(
      cv1_2 = cv2_expt_new, n1 = n_treatment,
      cv2_2 = cv2_cont_new, n2 = n_control
    )
  ) %>%
  mutate(diff = v_lnrr_laj - v_lnrr_laj_1B)


####### model with mean + summary + plot #######
method_MC_meta_mean <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "mean")
)

summary(method_MC_meta_mean)
summary_mean_2 <- orchaRd::mod_results(method_MC_meta_mean, mod = "treatment_overall", group = "title")
plot_mean_2 <- orchaRd::orchard_plot(summary_mean_2,
  xlab = "LogRR") + ggtitle("Mean")

plot_mean_2

####### model with median + summary + plot #######
method_MC_meta_median <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat %>% filter(m_measure == "median")
)

summary(method_MC_meta_median)
summary_median_2 <- orchaRd::mod_results(method_MC_meta_median, mod = "treatment_overall", group = "title")
plot_median_2 <- orchaRd::orchard_plot(summary_median_2,
  xlab = "LogRR") + ggtitle("Median")

plot_median_2

####### model with all + summary + plot #######
method_MC_meta_all <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj,
  test = "t",
  method = "REML",
  random = list(~ 1 | species, ~ 1 | title, ~ 1 | id),
  data = dat
)

summary(method_MC_meta_all)
summary_all_2 <- orchaRd::mod_results(method_MC_meta_all, mod = "treatment_overall", group = "title")
plot_all_2 <- orchaRd::orchard_plot(summary_all_2,
  xlab = "LogRR") + ggtitle("All")

plot_all_2

####### put all plots together #######
all_plots <- (plot_mean_2 | plot_median_2 | plot_all_2) + plot_annotation(tag_levels = "A")

ggsave("MC_plots.png", all_plots, dpi = 600, width = 20, height = 8)
