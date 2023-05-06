library("RSocrata")
library(lubridate)
library(tidytable)
library(purrr)

df <- read.socrata(
  "https://data.cdc.gov/resource/54ys-qyzm.json",
  app_token = "urLMJ5lvLUP65hAnK2ttm04FN",
  email     = "275261464@qq.com",
  password  = "rAzRxX:hdb-!E4Q"
)

dfh <- read.socrata(
  "https://data.cdc.gov/resource/k3na-u7xf.json",
  app_token = "urLMJ5lvLUP65hAnK2ttm04FN",
  email     = "275261464@qq.com",
  password  = "rAzRxX:hdb-!E4Q"
) %>% filter(month == "December 2022")

dfh1 = dfh %>% 
  mutate(rate = as.numeric(rate)) %>% 
  mutate(rate = 3/4 * rate)


cdc = df %>% 
  select(outcome, mmwr_week, age_group, vaccination_status, 
    vaccinated_with_outcome, vaccinated_population, unvaccinated_with_outcome, unvaccinated_population) %>%
  mutate(vaccinated_population = vaccinated_population %>% as.numeric(), 
    vaccinated_with_outcome = vaccinated_with_outcome %>% as.numeric(),
    unvaccinated_population = unvaccinated_population %>% as.numeric(),
    unvaccinated_with_outcome = unvaccinated_with_outcome %>% as.numeric(),
    mmwr_week = mmwr_week %>% as.numeric()) %>%
  filter(mmwr_week >= 202301 & mmwr_week <= 202303, age_group != "all_ages") %>% 
  mutate(mmwr_week = mmwr_week - min(mmwr_week)) %>%
  mutate(Group = case_when(age_group == "0.5-4" & vaccination_status == "vaccinated" ~ 1, 
                           age_group == "0.5-4" & vaccination_status == "vax with updated booster" ~ 2,
                           age_group == "5-11" & vaccination_status == "vaccinated" ~ 4,
                           age_group == "5-11" & vaccination_status == "vax with updated booster" ~ 5,
                           age_group == "12-17" & vaccination_status == "vaccinated" ~ 7,
                           age_group == "12-17" & vaccination_status == "vax with updated booster" ~ 8,
                           age_group == "18-29" & vaccination_status == "vaccinated" ~ 10,
                           age_group == "18-29" & vaccination_status == "vax with updated booster" ~ 11, 
                           age_group == "30-49" & vaccination_status == "vaccinated" ~ 13,
                           age_group == "30-49" & vaccination_status == "vax with updated booster" ~ 14,
                           age_group == "50-64" & vaccination_status == "vaccinated" ~ 16,
                           age_group == "50-64" & vaccination_status == "vax with updated booster" ~ 17,
                           age_group == "65-79" & vaccination_status == "vaccinated" ~ 19,
                           age_group == "65-79" & vaccination_status == "vax with updated booster" ~ 20,
                           age_group == "80+" & vaccination_status == "vaccinated" ~ 22,
                           age_group == "80+" & vaccination_status == "vax with updated booster" ~ 23)) %>%
  mutate(Outcome = case_when(outcome == "case" ~ 2, 
                                        "death" ~ 9)) %>%
  arrange(mmwr_week, Group)

cdc1 = cdc %>% 
  select(Outcome, mmwr_week, vaccinated_with_outcome, vaccinated_population, Group)


rowwise_function = function(time, outcome, group, n, scale_level){
  Outcome = outcome
  Group = group
  Time = runif(ceiling(n/scale_level), time, time + 1)
  return(tidytable(Time = Time, Outcome = Outcome, Group = Group))
}
vaccination_function = function(table, scale_level){
  n = table %>% pull(2) %>% diff() %>% abs()
  if (identical(n, numeric(0))) {
    Time = numeric(0)
    Outcome = numeric(0)
    Group = numeric(0)
    return(tidytable(Time = Time, Outcome = Outcome, Group = Group))
  } else{
    Group = table %>% distinct(Group) %>% pull()
    Time = n %>% imap(~(runif(ceiling(.x/scale_level), .y - 1, .y))) %>% unlist()
    return(tidytable(Time = Time, Outcome = 10, Group = Group))
  }
}


scale_level = 10

# vaccinated groups data frame --------------------------------------------

vd = cdc1 %>% 
  pmap_dfr(~rowwise_function(..2, ..1, ..5, ..3, scale_level)) %>% 
  arrange(Time)

vi = vd %>% group_by(Group) %>% summarise(total = ceiling(n(Time) * 1.5)) %>%
  pmap_dfr(~tidytable(Time = runif(..2, -10, 0), Outcome = 1, Group = ..1)) %>% 
  arrange(Time)

vt = vd %>% group_by(Group) %>% summarise(total = ceiling(n(Time) * 1.5)) %>%
  pmap_dfr(~tidytable(Time = runif(..2, -5, 2), Outcome = -1, Group = ..1)) %>% 
  arrange(Time)

vv = cdc1 %>% 
  select(mmwr_week, vaccinated_population, Group) %>% 
  filter(Group %in% seq(2, 23, by = 3)) %>%
  mutate(Group = Group - 1) %>%
  group_split(Group, .keep = T, .named = T) %>% 
  map_dfr(~vaccination_function(.x, scale_level)) %>%
  arrange(Time)

vh = cdc1 %>% 
  filter(mmwr_week == 0) %>%
  mutate(Hosp = case_when(Group %in% c(1, 2, 5, 8) ~ 0, 
                          Group == 4  ~ 1.58, 
                          Group == 7  ~ 2.03, 
                          Group == 10 ~ 9.68, 
                          Group == 11 ~ 3.38, 
                          Group == 13 ~ 9.68, 
                          Group == 14 ~ 3.38, 
                          Group == 16 ~ 20.6, 
                          Group == 17 ~ 6.22, 
                          Group == 19 ~ 90.8, 
                          Group == 20 ~ 37, 
                          Group == 22 ~ 90.8, 
                          Group == 23 ~ 37)) %>%
  mutate(Hosp = round(Hosp * vaccinated_population/(100000 * scale_level))) %>%
  select(Hosp, Group) %>%
  pmap_dfr(~tidytable(Time = runif(..1, 1, 3), Outcome = 4, Group = ..2)) %>%
  arrange(Time)
# un-vaccinated groups data frame --------------------------------------------
ugroup = cdc %>% 
  distinct(mmwr_week, age_group, .keep_all = T) %>% 
    mutate(Group = case_when(age_group == "0.5-4" ~ 0, 
                             age_group == "5-11"  ~ 3,
                             age_group == "12-17" ~ 6,
                             age_group == "18-29" ~ 9,
                             age_group == "30-49" ~ 12,
                             age_group == "50-64" ~ 15,
                             age_group == "65-79" ~ 18,
                             age_group == "80+"   ~ 21)) %>%
  select(Outcome, mmwr_week, unvaccinated_with_outcome, unvaccinated_population, Group)

ud = ugroup %>%
  pmap_dfr(~rowwise_function(..2, ..1, ..5, ..3, scale_level)) %>% 
  arrange(Time) %>% 
  mutate(Time = Time)

ui = ud %>% group_by(Group) %>% summarise(total = ceiling(n(Time) * 1.5)) %>%
  pmap_dfr(~tidytable(Time = runif(..2, -10, 0), Outcome = 1, Group = ..1)) %>% 
  arrange(Time)

ut = ud %>% group_by(Group) %>% summarise(total = ceiling(n(Time) * 1.5)) %>%
  pmap_dfr(~tidytable(Time = runif(..2, -5, 2), Outcome = -1, Group = ..1)) %>% 
  arrange(Time)



uv = ugroup %>% 
  select(mmwr_week, unvaccinated_population , Group) %>% 
  filter(Group %in% seq(0, 21, by = 3)) %>%
  group_split(Group, .keep = T, .named = T) %>% 
  map(~arrange(.x, mmwr_week)) %>%
  map_dfr(~vaccination_function(.x, scale_level)) %>%
  arrange(Time)



uh = ugroup %>% 
  filter(mmwr_week == 0) %>%
  mutate(Hosp = case_when(Group == 0  ~ 0, 
                          Group == 3  ~ 7.35, 
                          Group == 6  ~ 28, 
                          Group == 9  ~ 73.1, 
                          Group == 12 ~ 73.1, 
                          Group == 15 ~ 149, 
                          Group == 18 ~ 474, 
                          Group == 21 ~ 474)) %>%
  mutate(Hosp = round(Hosp * unvaccinated_population / (100000 * scale_level))) %>%
  select(Hosp, Group) %>%
  pmap_dfr(~tidytable(Time = runif(..1, 1, 3), Outcome = 4, Group = ..2)) %>%
  arrange(Time)

head_tail = tidytable(Time = c(-11, 4), Outcome = c(1, 11), Group = c(15, 0))

Final = ui %>% bind_rows(ut, uv, uh, ud, vt, vv, vh, vd, vi, head_tail) %>% arrange(Time) %>% mutate(Time = jitter(Time))
VN = cdc %>% filter(mmwr_week == 0) %>% pull(vaccinated_population, name = Group)
UN = ugroup %>% filter(mmwr_week == 0) %>% pull(unvaccinated_population, name = Group)
subset = seq(0, 23) %>% as.character()
N = (c(UN, VN)[subset] / scale_level) %>% ceiling()

library(readr)
write_csv(Final, "CDC.csv", col_names = F)



