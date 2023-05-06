# Library
library(ggplot2)
library(readr)
library(tidytable)
library(latex2exp)
library(cowplot)
library(ggridges)


alpha = read_csv("Data/alpha.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "alpha")

eta = read_csv("Data/eta.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "eta")


tau = read_csv("Data/tau.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "tau")

theta = read_csv("Data/theta.csv", col_names = F) %>% 
  rename("All Groups" = X1) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "theta")


gamma = read_csv("Data/gamma.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "gamma")

p_gamma = read_csv("Data/pgamma.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "pgamma")


omega = read_csv("Data/omega.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "omega")

beta = read_csv("Data/beta.csv", col_names = F) %>% 
  rename("Not Vaccinated" = X1, "Fully Vaccinated" = X2, "Boosted" = X3) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "beta")

ptheta = read_csv("Data/ptheta.csv", col_names = F) %>% 
  rename("All Groups" = X1) %>% 
  slice(900000:1000000) %>%
  pivot_longer(names_to = "Group", values_to = "Value") %>%
  mutate(Class = "ptheta")


beta_plot = ggplot(beta, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\tilde{\\beta}$")) +
  ylab("") +
  xlab("") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = -0.2))


alpha_plot = ggplot(alpha, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\alpha$")) +
  ylab("") +
  xlab("") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = -0.2))

eta_plot = ggplot(eta, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\eta$")) +
  ylab("") + 
  xlab("") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = -0.2))

tau_plot = ggplot(tau, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\tau$")) +
  ylab("") +
  xlab("") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = -0.2))

pgamma_plot = ggplot(p_gamma, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$p_\\gamma$")) +
  ylab("") +
  xlab("") + 
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -0.2))

gamma_plot = ggplot(gamma, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\gamma$")) +
  ylab("") +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -0.2))



ridge_all = plot_grid(beta_plot, alpha_plot, eta_plot, tau_plot, pgamma_plot, gamma_plot, align = "v", ncol = 1)
ridge_all
ggsave("ridge_all.png", ridge_all, height = 10)

theta_plot = ggplot(theta, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$\\theta$")) +
  ylab("") +
  xlab("") +
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -0.2))

ptheta_plot = ggplot(ptheta, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme_bw() + 
  labs(title = TeX("$p_\\theta$")) +
  ylab("") +
  xlab("value") +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5, vjust = -0.2))

ridge_another = plot_grid(theta_plot, ptheta_plot, align = "v", ncol = 1)
ridge_another
ggsave("ridge_another.png", ridge_another)