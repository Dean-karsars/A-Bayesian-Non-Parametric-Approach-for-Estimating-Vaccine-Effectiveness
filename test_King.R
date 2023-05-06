library(dtw)
library(readr)
library(coda)
library(tidytable)
library(purrr)

# Data Read in ------------------------------------------------------------
beta = read_csv("Data/beta.csv", col_names = FALSE) %>% 
  rename_with( ~paste0("beta", seq_along(.)), starts_with("X")) %>% 
  slice(600000:1000000)  %>% 
  map_dbl(~mean(.x))
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
alignment = betat[c(2, 3)] %>% map2(betat[c(1)], ~dtw(.y$betat, .x$betat, .keep = T, step.pattern = asymmetric))

index = alignment %>% map(~list(index1 = (.x$index1) %>% as.integer(), index2 = (.x$index2) %>% as.integer()))


# Comp --------------------------------------------------------------------

comparison = function(betat, index, base_name, name){
  Table1 = betat[[base_name]] %>% slice(index[[name]]$index1) %>% rename(base = betat)
  Table2 = betat[[name]] %>% slice(index[[name]]$index2) %>% rename(comp = betat) %>% select(comp)
  Table3 = Table1 %>% bind_cols(Table2) %>% mutate(Curve = (base - comp) / base) %>% select(Time, Curve) %>% mutate(Group = name)
  return(Table3)
}


curve1 = comparison(betat, index, "beta1" ,"beta2") 
curve2 = comparison(betat, index, "beta1" ,"beta3") 
curve = curve1 %>% 
  bind_rows(curve2) %>% 
  mutate(Group = as.factor(Group)) %>%
  mutate(Group = case_when(Group == "beta2" ~ "Fully Vaccinated", 
                           Group == "beta3" ~ "Boosted"))  %>%
  filter(Curve >= 0)

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


# beta-t curve --------------------------------------------------------------------

A = betat$beta1 %>% mutate(Group = "Not Vaccinated") %>% slice(seq(1, nrow(.), by = 50))
B = betat$beta2 %>% mutate(Group = "Fully Vaccinated") %>% slice(seq(1, nrow(.), by = 50))
C = betat$beta3 %>% mutate(Group = "Boosted")  %>% slice(seq(1, nrow(.), by = 50))
Beta_t_Curve = bind_rows(A, B, C) 
remove(A, B, C)
BT = ggplot(Beta_t_Curve, aes(x = Time, y = betat, color = Group)) +
  geom_line() + 
  geom_point() +
  ggtitle("Transmission Risk Over Time") +
  theme_bw()
BT
ggsave("BT.png", BT)



