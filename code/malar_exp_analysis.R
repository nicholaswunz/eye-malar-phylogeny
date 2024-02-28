# Load library
library(Matrix) # tidyverse not loading without it
library(tidyverse)
library(cowplot) # organise figures
library(colorspace)
library(brms)
library(rstan)

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(0.2, 0.2, 0.2, 0.2), units = , "cm"), 
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

# Set directory
setwd('/Users/nicholaswu/Dropbox/Eye malar study')

# Load data
clean_dat <- read.csv("exp_raw_data.csv") %>%
  dplyr::mutate(light = factor(light, levels = c("low", "medium", "high")),
                treatment = factor(treatment, levels = c("control", "painted", "clear")),
                body_bitten = factor(body_bitten, levels = c("front", "back")),
                time_strike_prop = time_strike_s / 300)
str(clean_dat)

clean_dat %>%
  ggplot(aes(x = svl_cm, y = mass_g, colour = sex)) +
  geom_point() +
  mytheme()

clean_dat %>%
  dplyr::group_by(light, treatment) %>%
  summarise(median = median(time_strike_s))


## STRIKE ATTEMPT ANALYSIS ## ------------------------------------------------------------------------------------------
# General STAN specs
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores avaliable to use
set.seed(10)

# Account for variation in individual prey choice (1|juvenile_ID) and individual variation in behaviour of the test subjects (1|adult_ID)
# No covariates
strike_model <- brms::brm(strike_attempt ~ light + treatment + sex + (1 | ID), 
                          data = clean_dat,
                          family = 'bernoulli', 
                          prior = set_prior('normal(0, 3)'), 
                          iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, 
                          control = list(adapt_delta = .99, max_treedepth = 15))

summary(strike_model)
brms::pp_check(strike_model)
plot(conditional_effects(strike_model), points = TRUE)

## PLOT FIG 2a ##
strike_marg_eff <- as.data.frame(ggeffects::ggpredict(strike_model, terms = c("light", "treatment"))) 

strike_plot <- strike_marg_eff %>%
  dplyr::rename(light = x, 
                treatment = group) %>%
  mutate(predicted = predicted * 100,
         conf.low = conf.low * 100,
         conf.high = conf.high * 100) %>%
  ggplot(aes(x = light, y = predicted, colour = light, shape = light)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 0.8, width = 0.1, show.legend = FALSE) + 
  geom_point(size = 3, show.legend = FALSE) +
  xlab(NULL) + ylab("Probability of successful strike (%)") +
  ylim(50,100) +
  colorspace::scale_colour_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  facet_wrap(~ treatment) +
  mytheme()


 ## TIME TO STRIKE ANALYSIS ## ----------------------------------------------------------------------------------------------------
# brms model
time_model <- brms::brm(time_strike_prop ~ light + treatment + sex + (1 | ID), 
                       data = clean_dat, 
                       family = zero_one_inflated_beta(), 
                       iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, 
                       control = list(adapt_delta = .99, max_treedepth = 15)) 

summary(time_model)
brms::pp_check(time_model)

as.data.frame(ggeffects::ggpredict(time_model, terms = c("light", "treatment"))) 

## PLOT FIG 2b ##
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

time_plot <- clean_dat %>%
  ggplot(aes(x = light, y = time_strike_s, fill = light)) +
  geom_jitter(aes(colour = light), position = position_jitter(0.05), size = 2, show.legend = FALSE) +
  geom_flat_violin(alpha = 0.5, position = position_nudge(x = 0.1), show.legend = FALSE) + 
  geom_hline(yintercept = 300, linetype = "dashed") + 
  colorspace::scale_fill_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  xlab(NULL) + ylab("Time to first strike (s)") +
  facet_wrap(~ treatment) +
  scale_y_reverse() +
  mytheme()

## NUMBER OF STRIKE ATTEMPTS ## -------------------------------------------
attempt_model <- brms::brm(strike_n ~ light + treatment + sex + (1 | ID), 
                        data = clean_dat, 
                        family = 'hurdle_poisson', 
                        iter = 5e3, warmup = 2.5e3, chains = 4, cores = 4, 
                        control = list(adapt_delta = .99, max_treedepth = 15)) 

summary(attempt_model)
brms::pp_check(attempt_model)

attempt_marg_eff <- as.data.frame(ggeffects::ggpredict(attempt_model, terms = c("light", "treatment"))) 

attempt_plot <- attempt_marg_eff %>%
  dplyr::rename(light = x, 
                treatment = group) %>%
  ggplot(aes(x = light, y = predicted, colour = light, shape = light)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 0.8, width = 0.1, show.legend = FALSE) + 
  geom_point(size = 3, show.legend = FALSE) +
  xlab(NULL) + ylab("Number of strike attempts") +
  colorspace::scale_colour_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  facet_wrap(~ treatment) +
  mytheme()

attempt_plot <- clean_dat %>%
  ggplot(aes(x = light, y = strike_n, fill = light)) +
  geom_point(aes(colour = light), size = 2, show.legend = FALSE) +
  geom_flat_violin(alpha = 0.5, position = position_nudge(x = 0.1), show.legend = FALSE) + 
  colorspace::scale_fill_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "SunsetDark", rev = FALSE) +
  xlab(NULL) + ylab("Number of strike attempts") +
  facet_wrap(~ treatment) +
  mytheme()

# CREATE FIG 2
plot_grid(time_plot, attempt_plot, 
          align = "h", axis = "bt", nrow = 2, labels = "auto")

clean_dat %>%
  group_by(light, treatment, strike_n) %>%
  count() %>%
  data.frame()
