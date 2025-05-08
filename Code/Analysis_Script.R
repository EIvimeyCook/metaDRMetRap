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
library(showtext)
library(ggbeeswarm)
library(ggnewscale)
library(cowplot)
library(rphylopic)
library(ggthemes)
font_add_google("Montserrat")
showtext_auto()
showtext_opts(dpi = 600)


sink("sessionInfo.txt")
sessionInfo()
sink()
####### read in data and create cv for everything if it exists #######
dat <- read_csv("Data/analysis_data.csv") %>%
  mutate(
    cv_control = na_if(sd_control / m_control, Inf),
    cv_experimental = na_if(sd_treatment / m_treatment, Inf)
  ) %>%
  janitor::clean_names() %>%
  mutate(
    treatment_type = as.factor(treatment_type),
    treatment_overall = as.factor(treatment_overall),
    species = as.factor(species),
    sex = as.factor(sex),
    id = as.factor(id),
    m_measure = as.factor(m_measure),
    title = as.factor(title),
    mice900_keep = as.factor(mice900_keep)
  ) %>%
  mutate(measure_plot = case_when(
    m_measure == "mean" ~ "Means",
    m_measure == "median" ~ "Medians"
  )) %>%
  mutate(sex_plot = case_when(
    sex == "male" ~ "Male",
    sex == "female" ~ "Female",
    sex == "mixed" ~ "Mixed"
  )) %>%
  mutate(species = fct_relevel(species, c("dogs", 
                          "mice",
                          "mouse lemur",
                          "rats",
                          "rhesus monkeys",
                          "stickleback",
                         "redtail killifish", 
                         "turqoise killifish")))

####### summarise overall number of titles,measures,species,sex,treatments #######
length(unique(dat$title))

dat %>%
  group_by(m_measure) %>%
  summarise(
    n = n(),
    n_stud = n_distinct(title)
  )

dat %>%
  group_by(species) %>%
  summarise(
    n = n(),
    n_stud = n_distinct(title)
  )

dat %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    n_stud = n_distinct(title)
  )

dat %>%
  group_by(treatment_overall) %>%
  summarise(
    n = n(),
    n_stud = n_distinct(title)
  )

dat %>%
  group_by(treatment_type) %>%
  summarise(
    n = n(),
    n_stud = n_distinct(title)
  )

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
  )

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
) %>%
  mutate(label = paste0(
    formatC(estimate, format = "f", digits = 3), " [",
    formatC(lowerCL, format = "f", digits = 3), ", ",
    formatC(upperCL, format = "f", digits = 3), "]"
  ))

# rphylopic images
mice_image <- pick_phylopic(name = "Mus musculus", n = 1, view = 1)
rat_image <- pick_phylopic(name = "Rattus norvegicus", n = 1, view = 1)
dog_image <- pick_phylopic(name = "Canis lupus familiaris", n = 1, view = 1)
killifish1_image <- pick_phylopic(name = "Nothobranchius", n = 1, view = 1)
killifish2_image <- pick_phylopic(name = "Nothobranchius furzeri", n = 1, view = 1)
rhesus_image <- pick_phylopic(name = "Macaca mulatta", n = 1, view = 1)
lemur_image <- pick_phylopic(name = "Microcebus murinus", n = 1, view = 1)
sb_image <- pick_phylopic(name = "Gasterosteus aculeatus", n = 1, view = 1)

get_attribution(uuid = "128043bd-2d06-47bc-a53b-b6f0bcf214e4")



ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat,
    aes(
      x = lnrr_laj,
      y = fct_rev(treatment_overall),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B)
    ),
    alpha = 0.4, show.legend = F
  ) +
  labs(colour = "Species") +
  ggthemes::scale_colour_colorblind() +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods,
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.5,
    linewidth = 1.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods, aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  ) +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 25, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 21),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "inside", 
    legend.position.inside =  c(0.8, 0.5),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 150, 20)
  ) +
  guides(size = "none") +
  #ggrepel::geom_text_repel(
   # data = all_mods,
    #aes(
    #  y = fct_rev(name),
     # x = upperCL,
     # label = label,
      #colour = model,
      #shape = data
    #),
    #position = position_dodge(width = 0.6),
    #show.legend = F,
    #size = ,
    ##hjust = -1.2,
    #point.padding = 0.3,
    #box.padding = 0.3,
    #segment.curvature = -0.3,
    #segment.ncp = 3,
    #segment.angle = 20
  #) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top"),
    shape = guide_legend(reverse = TRUE, title.position = "top")
  ) +
  scale_x_continuous(limits = c(-2, 2)) -> overall

legend_species <- cowplot::get_legend(
  ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = dat,
      aes(
        x = lnrr_laj,
        y = fct_rev(treatment_overall),
        colour = species
      ),
      alpha = 0.4
    ) +
    theme_bw() +
    theme(
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.spacing.x = unit(10, "cm"),
      legend.position = "right",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.width = unit(8, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 20),
      text = element_text(size = 20, family = "Montserrat"),
    ) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggthemes::scale_colour_colorblind() +
    labs(colour = "Species")
)

plot_with_species_inside <- cowplot::ggdraw(overall) +
  cowplot::draw_grob(legend_species, x = 0.5, y = -0.43, hjust = 0.5)


ggsave("Overall.pdf", plot_with_species_inside, dpi = 600, width = 22, height = 12.5)


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
  lnrr_laj ~ treatment_overall + inv_n_tilda +
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
sex_AC_Rap_mean <- rma.mv(lnrr_laj ~ sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
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
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
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
sex_AC_Rap_median <- rma.mv(lnrr_laj ~sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
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
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list( ~ 1 | title, ~ 1 | id),
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
sex_AC_Rap_all <- rma.mv(lnrr_laj ~ sex,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
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
  lnrr_laj ~ sex + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
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
sex_AC_Met_mean <- rma.mv(lnrr_laj ~ sex,
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
sex_AC_Met_all <- rma.mv(lnrr_laj ~ sex,
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
sex_AC_DR_mean <- rma.mv(lnrr_laj ~ sex,
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
  ),
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
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
) %>%
  mutate(label = paste0(
    formatC(estimate, format = "f", digits = 3), " [",
    formatC(lowerCL, format = "f", digits = 3), ", ",
    formatC(upperCL, format = "f", digits = 3), "]"
  ))

ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(treatment_overall == "Dietary restriction"),
    aes(
      x = lnrr_laj,
      y = fct_rev(sex_plot),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B),
      group = treatment_overall
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  ggthemes::scale_colour_colorblind() +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods %>% filter(treatment == "Dietary restriction"),
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.2,
    linewidth = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods %>% filter(treatment == "Dietary restriction"),
    aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 1.4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  ) +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 25, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 21),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none", 
    legend.position.inside =  c(0.8, 0.5),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(size = "none") +
  scale_x_continuous(limits = c(-2, 2)) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5)),
    shape = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5))
  )   +
  ggtitle("Dietary Restriction") -> sex_effect_plot_dr



ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(treatment_overall == "Metformin"),
    aes(
      x = lnrr_laj,
      y = fct_rev(sex_plot),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B),
      group = treatment_overall
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442")) +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods %>% filter(treatment == "Metformin"),
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.2,
    linewidth = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods %>% filter(treatment == "Metformin"),
    aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 1.4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  ) +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 25, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 21),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "inside", 
    legend.position.inside =  c(0.8, 0.5),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(size = "none") +
  scale_x_continuous(limits = c(-2, 2)) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5)),
    shape = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5))
  ) +
  ggtitle("Metformin") -> sex_effect_plot_met


ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(treatment_overall == "Rapamycin"),
    aes(
      x = lnrr_laj,
      y = fct_rev(sex_plot),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B),
      group = treatment_overall
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  scale_colour_manual(values = c("#E69F00")) +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods %>% filter(treatment == "Rapamycin"),
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.2,
    linewidth = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods %>% filter(treatment == "Rapamycin"),
    aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 1.4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  )  +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 25, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 21),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none", 
    legend.position.inside =  c(0.8, 0.5),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 150, 20)
  )  +
  guides(size = "none") +
  scale_x_continuous(limits = c(-2, 2)) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5)),
    shape = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5))
  ) +
  ggtitle("Rapamycin") -> sex_effect_plot_rap

all_sex_plot <- sex_effect_plot_dr/sex_effect_plot_met/sex_effect_plot_rap

legend_species <- cowplot::get_legend(
  ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = dat,
      aes(
        x = lnrr_laj,
        y = fct_rev(sex),
        colour = species
      ),
      alpha = 0.4
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.spacing.x = unit(4, "cm"),
      legend.spacing.y = unit(2, "cm"),
      legend.key.width = unit(3, "cm"),
      legend.key.height = unit(3, "cm"),
      legend.text = element_text(size = 14),
      text = element_text(size = 16, family = "Montserrat"),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.box.background = element_rect(fill = "transparent", colour = "transparent")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggthemes::scale_colour_colorblind() +
    labs(colour = "Species")
)

plot_with_species_inside <- cowplot::ggdraw(all_sex_plot) +
  cowplot::draw_grob(legend_species, x = 0.5, y = -0.46, hjust = 0.5)

ggsave("Sex.pdf", plot_with_species_inside, dpi = 600, width = 12.5, height = 23)


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
) %>%
  mutate(label = paste0(
    formatC(estimate, format = "f", digits = 3), " [",
    formatC(lowerCL, format = "f", digits = 3), ", ",
    formatC(upperCL, format = "f", digits = 3), "]"
  ))

ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(treatment_overall == "Dietary restriction"),
    aes(
      x = lnrr_laj,
      y = fct_rev(treatment_type),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B)
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  ggthemes::scale_colour_colorblind() +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods,
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.5,
    linewidth = 1.3,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods, aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  )+
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 25, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 21),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "inside", 
    legend.position.inside =  c(0.8, 0.5),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 150, 50)
  ) +
  guides(size = "none")  +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top"),
    shape = guide_legend(reverse = TRUE, title.position = "top")
  ) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 10))-> DR_Type

legend_species <- cowplot::get_legend(
  ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = dat,
      aes(
        x = lnrr_laj,
        y = fct_rev(treatment_type),
        colour = species
      ),
      alpha = 0.4
    ) +
    theme_bw() +
    theme(
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.spacing.x = unit(10, "cm"),
      legend.position = "right",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.key.width = unit(8, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 20),
      text = element_text(size = 20, family = "Montserrat"),
    ) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggthemes::scale_colour_colorblind() +
    labs(colour = "Species")
)

plot_with_species_inside <- cowplot::ggdraw(DR_Type) +
  cowplot::draw_grob(legend_species, x = 0.5, y = -0.43, hjust = 0.5)

ggsave("DR_type.pdf", plot_with_species_inside, dpi = 600, width = 24, height = 12.5)


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
  xlab = "LogRR"
) + ggtitle("Mean")

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
  xlab = "LogRR"
) + ggtitle("Median")

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
  xlab = "LogRR"
) + ggtitle("All")

plot_all_2

####### put all plots together #######
all_plots <- (plot_mean_2 | plot_median_2 | plot_all_2) + plot_annotation(tag_levels = "A")

ggsave("MC_plots.png", all_plots, dpi = 600, width = 20, height = 8)

##################### SM mice only DR #####################
####### model with mean + summary #######
method_AC_meta_mean_mice <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    species == "mice"
  )
)

summary(method_AC_meta_mean_mice)
summary_mean_mice <- orchaRd::mod_results(method_AC_meta_mean_mice,
  mod = "treatment_overall",
  group = "title"
)


####### model with pub bias mean + summary #######
method_AC_meta_mean_pub_mice <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    species == "mice"
  )
)

summary(method_AC_meta_mean_pub_mice)

summary_mean_pub_mice <- orchaRd::mod_results(method_AC_meta_mean_pub_mice,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

####### model with median + summary #######
method_AC_meta_median_mice <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% dplyr::filter(
    m_measure == "median",
    species == "mice"
  )
)


summary(method_AC_meta_median_mice)
summary_median_mice <- orchaRd::mod_results(method_AC_meta_median_mice,
  mod = "treatment_overall",
  group = "title"
)

####### model with pub bias median + summary #######
method_AC_meta_median_pub_mice <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    species == "mice"
  )
)

summary(method_AC_meta_median_pub_mice)

summary_median_pub_mice <- orchaRd::mod_results(method_AC_meta_median_pub_mice,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

####### model with all + summary #######
method_AC_meta_all_mice <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    species == "mice"
  )
)

summary(method_AC_meta_all_mice)

summary_all_mice <- orchaRd::mod_results(method_AC_meta_all_mice,
  mod = "treatment_overall",
  group = "title"
)

####### model with pub bias all + summary #######
method_AC_meta_all_pub_mice <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    species == "mice"
  )
)

summary(method_AC_meta_all_pub_mice)

summary_all_pub_mice <- orchaRd::mod_results(method_AC_meta_all_pub_mice,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)


####### model with mean + summary >keep900 #######
method_AC_meta_mean_mice_900 <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    species == "mice",
    mice900_keep == "Y"
  )
)

summary(method_AC_meta_mean_mice_900)
summary_mean_mice_900 <- orchaRd::mod_results(method_AC_meta_mean_mice_900,
  mod = "treatment_overall",
  group = "title"
)


####### model with pub bias mean + summary >keep900 #######
method_AC_meta_mean_pub_mice_900 <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "mean",
    species == "mice",
    mice900_keep == "Y"
  )
)

summary(method_AC_meta_mean_pub_mice_900)

summary_mean_pub_mice_900 <- orchaRd::mod_results(method_AC_meta_mean_pub_mice_900,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

####### model with median + summary >keep900 #######
method_AC_meta_median_mice_900 <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% dplyr::filter(
    m_measure == "median",
    species == "mice",
    mice900_keep == "Y"
  )
)


summary(method_AC_meta_median_mice_900)
summary_median_mice_900 <- orchaRd::mod_results(method_AC_meta_median_mice_900,
  mod = "treatment_overall",
  group = "title"
)

####### model with pub bias median + summary >keep900 #######
method_AC_meta_median_pub_mice_900 <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    m_measure == "median",
    species == "mice",
    mice900_keep == "Y"
  )
)

summary(method_AC_meta_median_pub_mice_900)

summary_median_pub_mice_900 <- orchaRd::mod_results(method_AC_meta_median_pub_mice_900,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)

####### model with all + summary >keep900 #######
method_AC_meta_all_mice_900 <- rma.mv(lnrr_laj ~ treatment_overall,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    species == "mice",
    mice900_keep == "Y"
  )
)

summary(method_AC_meta_all_mice_900)

summary_all_mice_900 <- orchaRd::mod_results(method_AC_meta_all_mice_900,
  mod = "treatment_overall",
  group = "title"
)

####### model with pub bias all + summary >keep900 #######
method_AC_meta_all_pub_mice_900 <- rma.mv(
  lnrr_laj ~ treatment_overall + inv_n_tilda +
    year.c,
  V = v_lnrr_laj_1B,
  test = "t",
  method = "REML",
  random = list(~ 1 | title, ~ 1 | id),
  data = dat %>% filter(
    species == "mice",
    mice900_keep == "Y"
  )
)

summary(method_AC_meta_all_pub_mice_900)

summary_all_pub_mice_900 <- orchaRd::mod_results(method_AC_meta_all_pub_mice_900,
  mod = "treatment_overall",
  group = "title",
  at = list(inv_n_tilda = 0, year.c = 0)
)


####### graph all effects #######

all_pub_model_mice <- summary_all_pub_mice$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures", studies = "No Studies Removed")

median_pub_model_mice <- summary_median_pub_mice$mod_table %>%
  mutate(model = "Adjusted", data = "Medians", studies = "No Studies Removed")

mean_pub_model_mice <- summary_mean_pub_mice$mod_table %>%
  mutate(model = "Adjusted", data = "Means", studies = "No Studies Removed")

all_model_mice <- summary_all_mice$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures", studies = "No Studies Removed")

median_model_mice <- summary_median_mice$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians", studies = "No Studies Removed")

mean_model_mice <- summary_mean_mice$mod_table %>%
  mutate(model = "Unadjusted", data = "Means", studies = "No Studies Removed")

all_pub_model_mice_900 <- summary_all_pub_mice_900$mod_table %>%
  mutate(model = "Adjusted", data = "All Measures", studies = "900-Day Only")

median_pub_model_mice_900 <- summary_median_pub_mice_900$mod_table %>%
  mutate(model = "Adjusted", data = "Medians", studies = "900-Day Only")

mean_pub_model_mice_900 <- summary_mean_pub_mice_900$mod_table %>%
  mutate(model = "Adjusted", data = "Means", studies  = "900-Day Only")

all_model_mice_900 <- summary_all_mice_900$mod_table %>%
  mutate(model = "Unadjusted", data = "All Measures", studies = "900-Day Only")

median_model_mice_900 <- summary_median_mice_900$mod_table %>%
  mutate(model = "Unadjusted", data = "Medians", studies = "900-Day Only")

mean_model_mice_900 <- summary_mean_mice_900$mod_table %>%
  mutate(model = "Unadjusted", data = "Means", studies = "900-Day Only")


all_mods_mice <- bind_rows(
  all_pub_model_mice,
  median_pub_model_mice,
  mean_pub_model_mice,
  all_model_mice,
  median_model_mice,
  mean_model_mice,
  all_pub_model_mice_900,
  median_pub_model_mice_900,
  mean_pub_model_mice_900,
  all_model_mice_900,
  median_model_mice_900,
  mean_model_mice_900
) %>%
  mutate(label = paste0(
    formatC(estimate, format = "f", digits = 3), " [",
    formatC(lowerCL, format = "f", digits = 3), ", ",
    formatC(upperCL, format = "f", digits = 3), "]"
  ))

ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(species == "mice"),
    aes(
      x = lnrr_laj,
      y = fct_rev(treatment_overall),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B),
      group = treatment_overall
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  scale_colour_manual(values = c("#E69F00")) +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods_mice %>% filter(studies == "No Studies Removed"),
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.2,
    linewidth = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods_mice %>% filter(studies == "No Studies Removed"),
    aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 1.4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  ) +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 18, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 18, face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(20, 20, 80, 20)
  ) +
  guides(size = "none") +
  ggrepel::geom_text_repel(
    data = all_mods_mice %>% filter(studies == "No Studies Removed"),
    aes(
      x = upperCL,
      y = fct_rev(name),
      label = label,
      colour = model,
      shape = data
    ),
    position = position_dodge(width = 0.6),
    show.legend = F,
    hjust = -1.2,
    point.padding = 0.6,
    box.padding = 0.5,
    segment.curvature = -0.3,
    segment.ncp = 3,
    segment.angle = 20,
    size = 3.5,
    fontface = "italic"
  ) +
  scale_x_continuous(limits = c(-2, 2)) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5)),
    shape = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5))
  ) +
  ggtitle("All Effects") -> mice_effect_all

ggplot() +
  ggbeeswarm::geom_quasirandom(
    data = dat %>% filter(species == "mice",mice900_keep == "Y"),
    aes(
      x = lnrr_laj,
      y = fct_rev(treatment_overall),
      colour = species,
      shape = measure_plot,
      size = 1 / sqrt(v_lnrr_laj_1B),
      group = treatment_overall
    ),
    alpha = 0.4,
    show.legend = F
  ) +
  labs(colour = "Species") +
  scale_colour_manual(values = c("#E69F00")) +
  ggnewscale::new_scale_color() +
  geom_errorbarh(
    data = all_mods_mice %>% filter(studies == "900-Day Only"),
    aes(
      x = estimate,
      y = fct_rev(name),
      xmin = lowerCL,
      xmax = upperCL,
      colour = model,
      shape = data
    ),
    height = 0.2,
    linewidth = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    data = all_mods_mice %>% filter(studies == "900-Day Only"),
    aes(
      x = estimate,
      y = fct_rev(name),
      colour = model,
      shape = data
    ),
    size = 1.4,
    stroke = 1.7,
    position = position_dodge(width = 0.6)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    colour = "Publication Bias",
    shape = "Measure",
    x = "Log Response Ratio",
    y = NULL,
    size = NULL
  ) +
  theme_bw() +
  scale_colour_manual(values = c("purple", "black")) +
  scale_shape_manual(values = c("circle", "triangle", "square")) +
  theme(
    axis.title.y = element_blank(),
    text = element_text(size = 18, family = "Montserrat"),
    legend.key.width = unit(1.6, "cm"),
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    axis.text.y = element_text(size = 18, face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = margin(20, 20, 80, 20)
  ) +
  guides(size = "none") +
  ggrepel::geom_text_repel(
    data = all_mods_mice %>% filter(studies == "900-Day Only"),
    aes(
      x = upperCL,
      y = fct_rev(name),
      label = label,
      colour = model,
      shape = data
    ),
    position = position_dodge(width = 0.6),
    show.legend = F,
    hjust = -1.2,
    point.padding = 0.6,
    box.padding = 0.5,
    segment.curvature = -0.3,
    segment.ncp = 3,
    segment.angle = 20,
    size = 3.5,
    fontface = "italic"
  ) +
  scale_x_continuous(limits = c(-2, 2)) +
  guides(
    colour = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5)),
    shape = guide_legend(reverse = TRUE, title.position = "top", override.aes = list(size = 5))
  ) +
  ggtitle("900 day rule only") -> mice_effect_900

overall_mice <- mice_effect_all + mice_effect_900 +
  plot_layout(guides = "collect")


ggsave("Overall_mice.png", overall_mice, dpi = 600, width = 23, height = 12.5)
