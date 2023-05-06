# Library
library(ggplot2)
library(readr)
library(tidytable)
library(latex2exp)

# create a data set
alpha = read_csv("Data/alpha.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "alpha")


gamma = read_csv("Data/gamma.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "gamma")

p_gamma = read_csv("Data/pgamma.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "pgamma")


omega = read_csv("Data/omega.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "omega")


eta = read_csv("Data/eta.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "eta")


tau = read_csv("Data/tau.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "tau")

theta = read_csv("Data/theta.csv", col_names = F) %>% 
  rename("All Groups" = X1) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "theta")

ptheta = read_csv("Data/ptheta.csv", col_names = F) %>% 
  rename("All Groups" = X1) %>% 
  slice(600000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "ptheta")

P_Table = bind_rows(alpha, gamma, p_gamma, omega)

# Most basic violin chart
p = ggplot(P_Table, aes(x = Group, y = Value, fill = Group)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() +
  theme_bw() + 
  facet_wrap(vars(Class), scales = "free_y")
p

ggsave("violin.png", p)
# Most basic ridge chart
library(ggridges)


P2_Table = bind_rows(eta, tau, theta, ptheta) 

p2 = ggplot(P2_Table, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") + 
  facet_wrap(vars(Class), scales = "free") + 
  theme_bw() 
  
p2
ggsave("ridge.png", p2)








