library(dtw)
library(readr)
library(coda)
library(tidytable)
library(purrr)

# Data Read in ------------------------------------------------------------
beta = read_csv("Data/beta.csv", col_names = FALSE) %>% 
  rename_with( ~paste0("beta", seq_along(.)), starts_with("X")) %>% 
  slice(900000:1000000)  %>% 
  map_dbl(~mean(.x))
Gaus = read_csv("Data/TimeTable.csv", col_names = FALSE) %>% 
  filter(X2 == 1 & !is.na(X4)) %>% 
  select(X1, X3, X4) %>% 
  rename(Time = X1, Group = X3, Gaus = X4) %>% 
  mutate(Group = paste0("beta", Group + 1)) %>%
  mutate(Gaus = (1 + exp(-Gaus))^(-1)) %>%
  group_split(Group, .named = T, .keep = F)


# Modify the data ---------------------------------------------------------
betat = Gaus[names(beta)] %>% 
  map2(beta, ~mutate(.x, beta_t = Gaus * .y) %>%
       select(Time, beta_t))

# construct comp data table ---------------------------------------------------------
alignment = list()


alignment$beta2 = betat$beta1 %>% 
  filter(Time >= min(betat$beta2$Time)) %>% 
  pull(beta_t) %>%
  dtw(betat$beta2$beta_t, step.pattern = asymmetric)

alignment$beta3 = betat$beta1 %>% 
  filter(Time >= min(betat$beta3$Time)) %>% 
  pull(beta_t) %>%
  dtw(betat$beta3$beta_t, step.pattern = asymmetric)

index = alignment %>% 
  map(~list(index1 = (.x$index1) %>% as.integer(), index2 = (.x$index2) %>% as.integer()))

comparison = function(betat, index, base_name, name){
  Table1 = betat[[base_name]] %>% filter(Time >= (betat[[name]] %>% pull(Time) %>% min())) %>% slice(index[[name]]$index1) %>% rename(base = beta_t)
  Table2 = betat[[name]] %>% slice(index[[name]]$index2) %>% rename(comp = beta_t) %>% select(comp)
  Table3 = Table1 %>% bind_cols(Table2) %>% mutate(Curve = (base - comp) / base) %>% select(Time, Curve) %>% mutate(Group = name)
  return(Table3)
}


curve1 = comparison(betat, index, "beta1" ,"beta2") 
curve2 = comparison(betat, index, "beta1" ,"beta3") 

curve = curve1 %>% 
  bind_rows(curve2) %>% 
  filter(Curve >= 0) %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Group = case_when(Group == "beta2" ~ "Fully Vaccinated", 
                           Group == "beta3" ~ "Boosted")) 

# plot --------------------------------------------------------------------

library(ggplot2)
VE = ggplot(curve, aes(x = Time, y = Curve, color = Group)) +
  geom_smooth(method = "gam", se = TRUE) +
  ggtitle("Vaccine Effectiveness Over Time") +
  xlab("Days") +
  ylab("VE(t)") +
  theme_bw() 

VE
ggsave("VE.png", VE)