library(dtw)
library(readr)
library(coda)
library(tidytable)
library(purrr)

# Data Read in ------------------------------------------------------------
beta = read_csv("Data/beta.csv", col_names = FALSE) %>% rename_with( ~paste0("beta", seq_along(.)), starts_with("X"))  %>% map_dbl(~mean(.x))
Gaus = read_csv("Data/TimeTable.csv", col_names = FALSE) %>% 
  filter(X2 == 1 & !is.na(X4)) %>% 
  select(X1, X3, X4) %>% 
  rename(Time = X1, Group = X3, Gaus = X4) %>% 
  mutate(Group = paste0("beta", Group + 1)) %>%
  mutate(Gaus = (1 + exp(-Gaus))^(-1)) %>%
  group_split(Group, .named = T, .keep = F)


# Modify the data ---------------------------------------------------------
betat = Gaus %>% .[names(beta)] %>% map2(beta, ~mutate(.x, betat = Gaus * .y) %>% select(Time, betat))



# DTW ---------------------------------------------------------------------
alignment = betat[seq(1, 24)[-seq(1, 22, by = 3)]] %>% map2(betat[seq(1, 22, by = 3) %>% rep(2) %>% sort()], ~dtw(.y$betat, .x$betat))
index = alignment %>% map(~list(index1 = (.x$index1) %>% as.integer(), index2 = (.x$index2) %>% as.integer()))


# Comp --------------------------------------------------------------------

comparison = function(betat, index, base_name, name){
  Table1 = betat[[base_name]] %>% slice(index[[name]]$index1) %>% rename(base = betat)
  Table2 = betat[[name]] %>% slice(index[[name]]$index2) %>% rename(comp = betat) %>% select(comp)
  Table3 = Table1 %>% bind_cols(Table2) %>% mutate(Curve = (base - comp) / base) %>% select(Time, Curve) %>% mutate(Group = name)
  return(Table3)
}


curve1 = comparison(betat, index, "beta22" ,"beta23")
curve2 = comparison(betat, index, "beta22" ,"beta24")
curve = curve1 %>% bind_rows(curve2) %>% mutate(Group = as.factor(Group))

# plot --------------------------------------------------------------------

library(ggplot2)
ggplot(curve, aes(x = Time, y = Curve, color = Group)) +
  geom_smooth(method = "gam") +
  geom_point() +
  ggtitle("Vaccine Effectiveness over time") +
  xlab("Days") +
  ylab("VE(t)")


