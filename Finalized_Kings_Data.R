library(tidytable)
library(purrr)
Data = readRDS("Kings_Initial_Data")

Event_Table = function(x, outcome, group) {
 df = data.table(Time = x, outcome = outcome, group = group)
}

Data$Thinned$Not_Fully_Vax = runif(100, 10, 40) %>% jitter()
Data$Thinned$Fully_Vax = runif(100, 10, 40) %>% jitter()
Data$Thinned$Boosted = runif(100, 10, 40) %>% jitter()

Thinned = map2_dfr(Data$Thinned, c(1, 2, 3) - 1 , ~Event_Table(.x, -1, .y)) %>% arrange(Time)
I =  map2_dfr(Data$I, c(1, 2, 3) - 1, ~Event_Table(.x, 1, .y)) %>% arrange(Time)
D = map2_dfr(Data$D, c(1, 2, 3) - 1, ~Event_Table(.x, 2, .y))
U = map2_dfr(Data$U, c(1, 2, 3) - 1, ~Event_Table(.x, 3, .y))
H = map2_dfr(Data$H, c(1, 2, 3) - 1, ~Event_Table(.x, 4, .y))
DR = map2_dfr(Data$Rd, c(1, 2, 3) - 1, ~Event_Table(.x, 5, .y))
UR = map2_dfr(Data$Ru, c(1, 2, 3) - 1, ~Event_Table(.x, 6, .y))
Death = map2_dfr(Data$Death, c(1, 2, 3) - 1 , ~Event_Table(.x, 9, .y))
V = map2_dfr(Data$Vaccination[c(1, 2)], c(1, 2) - 1, ~Event_Table(.x, 10, .y))
Time_Last = data.table(Time = Data$Time_Last + 5, outcome = 11, group = 0)

Final = Thinned %>% bind_rows(I, D, U, H, DR, UR, Death, Time_Last) %>% arrange(Time)


library(readr)
write_csv(Final, "KingsData.csv", col_names = F)
