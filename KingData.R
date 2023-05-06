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
  filter(Date > ymd(20230201)) %>%
  mutate(Date = Date - min(Date)) %>%
  mutate(Date = as.numeric(Date)) 

#    
D1 = outcomes_daily %>% pull(Positive_Not_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
D2 = outcomes_daily %>% pull(Positive_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
D3 = outcomes_daily %>% pull(Positive_Boosted) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
D = list(Not_Fully_Vax = D1, Fully_Vax = D2, Boosted = D3) 
D %>% unlist() %>% duplicated() %>% sum()

H1 = outcomes_daily %>% pull(Hospitalization_Not_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
H2 = outcomes_daily %>% pull(Hospitalization_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F) %>% .[c(-1, -2)]
H3 = outcomes_daily %>% pull(Hospitalization_Boosted) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
H = list(Not_Fully_Vax = H1, Fully_Vax = H2, Boosted = H3)
H %>% unlist() %>% duplicated() %>% sum()

Death1 = outcomes_daily %>% pull(Death_Not_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
Death2 = outcomes_daily %>% pull(Death_Fully_Vax) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
Death3 = outcomes_daily %>% pull(Death_Boosted) %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
Death = list(Not_Fully_Vax = Death1, Fully_Vax = Death2, Boosted = Death3)
Death %>% unlist() %>% duplicated() %>% sum()

V1 = outcomes_daily %>% pull(Pop_Not_Fully_Vax) %>% diff() %>% abs() %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
V2 = outcomes_daily %>% pull(Pop_Boosted) %>% diff() %>% abs() %>% imap(~runif(.x, .y - 1, .y)) %>% unlist(use.names = F)
Vaccination = list(Not_Fully_Vax = V1, Fully_Vax = V2) %>% map(~jitter(.x))
Vaccination %>% unlist() %>% duplicated() %>% sum()


N_i = outcomes_daily %>% slice_head(1) %>% select(Pop_Not_Fully_Vax, Pop_Fully_Vax, Pop_Boosted) %>% as.list()
Time_Last = outcomes_daily %>% pull(Date) %>% max() + 5

Time = list(D = D, H = H, Death = Death, Vaccination = Vaccination)
Time %>% unlist() %>% duplicated() %>% sum()



Knwon = list(D = D, H = H, Death = Death, Vaccination = Vaccination,  N_i = N_i, Time_Last = Time_Last)
saveRDS(Knwon, "Kings_Knwon")


Knwon = readRDS("Kings_Knwon")
I1 = Knwon$D$Not_Fully_Vax %>% map(~.x + runif(1, -10, 0)) %>% unlist() %>% sort()
I1_1 = I1 %>% sample(0.5 * length(I1), FALSE) %>% jitter() %>% sort()
I1 = c(I1, I1_1) %>% jitter() %>% sort()
I2 = Knwon$D$Fully_Vax %>% map(~.x + runif(2, -9, 0)) %>% unlist() %>% sort()
I3 = Knwon$D$Boosted %>% map(~.x + runif(2, -8, 0)) %>% unlist() %>% sort()
I = list(Not_Fully_Vax = I1, Fully_Vax = I2, Boosted = I3)
I %>% unlist() %>% duplicated() %>% sum()

U = map(I, ~sample(.x, 0.01 * length(.x), F) %>% + runif(1, 0, 1))
U %>% unlist() %>% duplicated() %>% sum()

Rd = map(Knwon$D, ~sample(.x, 0.01 * length(.x), F) + runif(1, 1, 10)) %>% map(~subset(.x, .x < Knwon$Time_Last))
Rd %>% unlist() %>% duplicated() %>% sum()

Ru = map(U, ~sample(.x, 0.1 * length(.x), F) + runif(1, 0, 10)) %>% map(~subset(.x, .x < Knwon$Time_Last))
Ru %>% unlist() %>% duplicated() %>% sum()


Time = list(I, U, Knwon$D, Knwon$H, Knwon$Death, Knwon$Vaccination, Rd, Ru)
Time %>% unlist() %>% duplicated() %>% sum()

Initial_Data = list(I = I, U = U, 
                    D = Knwon$D, H = Knwon$H, 
                    Rd = Rd, Ru = Ru, 
                    Vaccination = Knwon$Vaccination, Time_Last = Knwon$Time_Last,
                    Death = Knwon$Death,
                    N_i = Knwon$N_i)

saveRDS(Initial_Data, "Kings_Initial_Data")


Initial_Data = readRDS("Kings_Initial_Data")

Time = Initial_Data %$% list(I ,U, D, H, Death, Vaccination, Rd, Ru)
Time %>% unlist() %>% duplicated() %>% sum()
