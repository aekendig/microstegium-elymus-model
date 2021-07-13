#### set-up ####

# goal: litter effects on establishment

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)

# import data
MvEstDat <- read_csv("../microstegium-litter-reu/intermediate-data/mv_establishment_data.csv")
EvEstDat <- read_csv("../microstegium-litter-reu/intermediate-data/ev_establishment_data.csv")


#### edit data ####

# remove competition treatments
MvEstDat2 <- MvEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstMvDenCor)

EvEstDat2 <- EvEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstEv)

# initial visualizations
ggplot(MvEstDat2, aes(x = Litter.g, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()

ggplot(EvEstDat2, aes(x = Litter.g, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()


#### fit models ####

# Mv
mv_bh_mod <- brm(data = MvEstDat2, family = gaussian,
                 bf(PropEst ~ e0/(1 + beta * Litter.g), 
                    e0 ~ 1, 
                    beta ~ 1, 
                    nl = T),
                 prior <- c(prior(normal(0.85, 1), nlpar = "e0", lb = 0),
                            prior(exponential(0.5), nlpar = "beta", lb = 0),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

prior_summary(mv_bh_mod)
summary(mv_bh_mod)
mv_bh_mod <- update(mv_bh_mod, chains = 3)
summary(mv_bh_mod)
plot(mv_bh_mod)

# Ev
ev_bh_mod <- update(mv_bh_mod, newdata = EvEstDat2)

summary(ev_bh_mod)
plot(ev_bh_mod)

# save
save(mv_bh_mod, file = "output/litter_reu_mv_establishment_model.rda")
save(ev_bh_mod, file = "output/litter_reu_ev_establishment_model.rda")


#### figure ####

# simulate data
bh_sim_dat <- tibble(Litter.g = seq(min(MvEstDat2$Litter.g), max(MvEstDat2$Litter.g), length.out = 200)) %>%
  mutate(PropEst = fitted(mv_bh_mod, newdata = .)[, "Estimate"],
         lower = fitted(mv_bh_mod, newdata = .)[, "Q2.5"],
         upper = fitted(mv_bh_mod, newdata = .)[, "Q97.5"],
         species = "M. vimineum") %>%
  full_join(tibble(Litter.g = seq(min(MvEstDat2$Litter.g), max(MvEstDat2$Litter.g), length.out = 200)) %>%
              mutate(PropEst = fitted(ev_bh_mod, newdata = .)[, "Estimate"],
                     lower = fitted(ev_bh_mod, newdata = .)[, "Q2.5"],
                     upper = fitted(ev_bh_mod, newdata = .)[, "Q97.5"],
                     species = "E. virginicus"))

# combine raw data
raw_dat <- MvEstDat2 %>%
  select(Litter.g, PropEst) %>%
  mutate(species = "M. vimineum") %>%
  full_join(EvEstDat2 %>%
              select(Litter.g, PropEst) %>%
              mutate(species = "E. virginicus"))

# plot
pdf("output/litter_reu_establishment_analysis.pdf", width = 6, height = 3)
ggplot(bh_sim_dat, aes(x = Litter.g, y = PropEst)) +
  geom_ribbon(alpha = 0.5, aes(ymin = lower, ymax = upper)) +
  geom_line() +
  stat_summary(data = raw_dat, fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(data = raw_dat, fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic")) +
  facet_wrap(~ species) +
  labs(x = expression(paste(italic("M. vimineum"), " litter (g)", sep = "")), y = "Proportion seeds established")
dev.off()
