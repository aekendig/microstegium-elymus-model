#### establishment statistics for Amy's model ####

# select data with Ev by itself
ev_bh_dat <- filter(EvEstDat1, SpPresent == "Ev")

# initial value
filter(ev_bh_dat, Treatment == "None") %>%
  summarise(mean_germ = mean(NewGermEv))

# fit model
# model
ev_bh_mod <- brm(data = ev_bh_dat, family = gaussian,
                 bf(NewGermEv ~ g0/(1 + alpha * Litter.g), 
                    g0 ~ 1, 
                    alpha ~ 1, 
                    nl = T),
                 prior <- c(prior(normal(42, 10), nlpar = "g0", lb = 0),
                            prior(exponential(0.5), nlpar = "alpha", lb = 0),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

prior_summary(ev_bh_mod)
summary(ev_bh_mod)
ev_bh_mod <- update(ev_bh_mod, chains = 3)
plot(ev_bh_mod)

# simulate data
ev_bh_sim_dat <- tibble(Litter.g = seq(min(ev_bh_dat$Litter.g), max(ev_bh_dat$Litter.g), length.out = 200)) %>%
  mutate(NewGermEv = fitted(ev_bh_mod, newdata = .)[, "Estimate"],
         lower = fitted(ev_bh_mod, newdata = .)[, "Q2.5"],
         upper = fitted(ev_bh_mod, newdata = .)[, "Q97.5"])

# plot raw
ggplot(ev_bh_dat, aes(x = Litter.g, y = NewGermEv)) +
  geom_point() +
  geom_ribbon(data = ev_bh_sim_dat, alpha = 0.5, aes(ymin = lower, ymax = upper)) +
  geom_line(data = ev_bh_sim_dat)

ggplot(ev_bh_dat, aes(x = Litter.g, y = NewGermEv/50)) +
  geom_point() +
  geom_ribbon(data = ev_bh_sim_dat, alpha = 0.5, aes(ymin = lower/50, ymax = upper/50)) +
  geom_line(data = ev_bh_sim_dat)
