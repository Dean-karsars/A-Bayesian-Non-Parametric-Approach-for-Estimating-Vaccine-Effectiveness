library(readr)
library(readxl)
library(tidytable)
library(lubridate)
library(purrr)
library(magrittr)
outcomes_daily <- read_excel("outcomes_daily.xlsx", 
                             col_types = c("text", "numeric", "numeric", 
                             "numeric", "numeric", "skip", "skip", 
                             "numeric", "skip", "numeric", "skip", 
                             "numeric", "skip", "numeric", "skip", 
                             "numeric", "numeric", "numeric", 
                             "numeric", "numeric", "skip", "numeric", 
                             "skip", "numeric", "skip", "numeric", 
                             "skip", "numeric", "skip")) %>% 
  as_tidytable() %>%
  mutate(Date = ymd(Date)) %>%
  select(Date, Pop_Not_Fully_Vax, Pop_Fully_Vax, Pop_Boosted, 
         Positive_Not_Fully_Vax, Positive_Fully_Vax, Positive_Boosted, 
         Hospitalization_Not_Fully_Vax, Hospitalization_Fully_Vax, Hospitalization_Boosted,
         Death_Not_Fully_Vax, Death_Fully_Vax, Death_Boosted) %>% 
  mutate(Date = Date - min(Date)) %>%
  mutate(Date = as.numeric(Date)) 

Time_Table = function(Table, Group, scale_level){
  
  D = Table %>% pull(3) %>% imap_dfr(~tidytable(Time = runif(ceiling(.x/scale_level), .y - 1, .y - 0.5), Outcome = 2, Group = Group))
  H = Table %>% pull(4) %>% imap_dfr(~tidytable(Time = runif(ceiling(.x/scale_level), .y - 0.5, .y), Outcome = 4, Group = Group))
  Death = Table %>% pull(5) %>% imap_dfr(~tidytable(Time = runif(ceiling(.x/scale_level), .y - 1, .y), Outcome = 9, Group = Group))
  
  m = Table %>% 
    pull(3) %>% 
    sum() %>% 
    divide_by(scale_level * 0.5) %>% 
    ceiling()
  
  I = D %>% mutate(Time = Time - 1, Outcome = 1, Group = Group)
  Th = D %>% mutate(Time = Time + 1, Outcome = -1, Group = Group)
  
  Final = bind_rows(D, H, Death, I, Th) %>% arrange(Time)
  return(Final)
}

scale_level = 100
# Un-vaccinated ------------------------------------------------------------

U = outcomes_daily %>% select(Date, Pop_Not_Fully_Vax, Positive_Not_Fully_Vax, Hospitalization_Not_Fully_Vax, Death_Not_Fully_Vax)
U_Table = U %>% Time_Table(0, scale_level)
UV = U %>% pull(2) %>% diff() %>% abs() %>% imap_dfr(~tidytable(Time = runif(ceiling(.x/scale_level), .y - 1, .y), Outcome = 10, Group = 0))

# vaccinated ------------------------------------------------------------
V = outcomes_daily %>% select(Date, Pop_Fully_Vax, Positive_Fully_Vax, Hospitalization_Fully_Vax, Death_Fully_Vax)
V_Table = V %>% Time_Table(1, scale_level)

# Boosted ------------------------------------------------------------
B = outcomes_daily %>% select(Date, Pop_Boosted, Positive_Boosted, Hospitalization_Boosted, Death_Boosted)
B_Table = B %>% Time_Table(2, scale_level)
VV = B %>% pull(2) %>% diff() %>% abs() %>% imap_dfr(~tidytable(Time = runif(ceiling(.x/scale_level), .y - 1, .y), Outcome = 10, Group = 1))

# Final ------------------------------------------------------------
head_tail = tidytable(Time = c(-10, max(outcomes_daily$Date) + 3), Outcome = c(1, 11), Group = 0)

Final = bind_rows(U_Table, V_Table, B_Table, UV, VV, head_tail) %>% mutate(Time = jitter(Time)) %>% arrange(Time)




N = outcomes_daily %>% select(Pop_Not_Fully_Vax, Pop_Fully_Vax, Pop_Boosted) %>% slice_head(1) %>% divide_by(scale_level) %>% ceiling()
library(readr)
write_csv(Final, "KingsData.csv", col_names = F)


