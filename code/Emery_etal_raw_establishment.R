#### set-up ####

# goal: estimate survival from seedling to adult stages, see how this varies with seedling density

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)

# import data
dat <- read_csv("data/Emery_etal_2013_raw.csv")


#### edit data ####

# check treatment spelling
unique(dat$trt)

# subset for reference plots (no treatment)
# make long by year
# sects = seedling counts, stems = proxy for adults
# calculate proportion survival
dat2 <- dat %>%
  filter(trt == "REF") %>%
  select(site, plot, trt, S09sects, F09stems) %>%
  rename("sects" = "S09sects", "stems" = "F09stems") %>%
  mutate(year = 2009) %>%
  full_join(dat %>%
              filter(trt == "REF") %>%
              select(site, plot, trt, S10sects, F10stems) %>%
              rename("sects" = "S10sects", "stems" = "F10stems") %>%
              mutate(year = 2010)) %>%
  mutate(survival = stems/sects,
         survival2 = ifelse(survival > 1, 1, survival),
         stems2 = ifelse(stems > sects, sects, stems),
         log_surv2 = log(survival2),
         no_stems2 = sects - stems2) %>%
  filter(!is.na(survival))

# separate by year
dat2_y1 <- dat2 %>%
  filter(year == 2009)

dat2_y2 <- dat2 %>%
  filter(year == 2010)


#### summarize data ####

# overall year mean
dat2 %>%
  group_by(year) %>%
  summarise(mean = mean(survival),
            sd = sd(survival),
            n = length(survival),
            mean2 = mean(survival2),
            sd2 = sd(survival2),
            mean_stems2 = mean(stems2))
# 0.38 and 0.21
# 0.36 and 0.21

# overall mean
dat2 %>%
  summarise(mean = mean(survival),
            sd = sd(survival),
            n = length(survival),
            mean2 = mean(survival2),
            sd2 = sd(survival2))
# 0.30
# 0.29


#### visualize ####

# survival by density
ggplot(dat2, aes(sects, survival)) +
  geom_point() +
  facet_wrap(~ year)

ggplot(dat2, aes(sects, survival2)) +
  geom_point() +
  facet_wrap(~ year)

# log-transformed
ggplot(dat2, aes(sects, log_surv2)) +
  geom_point() +
  facet_wrap(~ year, scales = "free")


#### model ####

# initial fit
emeryMod1 <- brm(stems2 | trials(sects) ~ sects + (1|year),
                  data = dat2, family = binomial,
                  prior = c(prior(normal(0, 10), class = Intercept),
                            prior(normal(0, 10), class = b),
                            prior(cauchy(0, 1), class = sd)),
                  iter = 6000, warmup = 1000, chains = 1)
# 60 divergent transitions
summary(emeryMod1)

# update model
emeryMod2 <- update(emeryMod1,
                    chains = 3,
                    control = list(adapt_delta = 0.999))
# more errors
summary(emeryMod2)
# coefficient is still zero

# normal fit
emeryMod3 <- brm(log_surv2 ~ sects + (1|year),
                 data = dat2, family = gaussian,
                 prior = c(prior(normal(0, 10), class = Intercept),
                           prior(normal(0, 10), class = b),
                           prior(cauchy(0, 1), class = sd)),
                 iter = 6000, warmup = 1000, chains = 1)
# 54 divergent transitions
summary(emeryMod3)
# still a near zero coefficient

# first year only
emeryMod4 <- brm(stems2 | trials(sects) ~ sects,
                 data = dat2_y1, family = binomial,
                 prior = c(prior(normal(110, 10), class = Intercept),
                           prior(normal(0, 10), class = b)),
                 iter = 6000, warmup = 1000, chains = 1)

summary(emeryMod4)
# still a near zero coefficient

# update model
emeryMod5 <- update(emeryMod4,
                    chains = 3)
summary(emeryMod5)
pp_check(emeryMod5, nsamples = 50)
# estimates are way off

# fit with glm
emeryMod6 <- glm(cbind(stems2, no_stems2) ~ sects,
                 data = dat2_y1,
                 family = binomial)
summary(emeryMod6)
# it is a small effect, but significant