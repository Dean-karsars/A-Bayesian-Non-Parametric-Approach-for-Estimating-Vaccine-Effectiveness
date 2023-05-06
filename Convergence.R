library(readr)
library(coda)
library(tidytable)

myplot = function(path, range, thinned){
  subset = seq(1, range, by = thinned)
  var = read_csv(path, col_names = FALSE)
  var %>% map( ~as.mcmc(.x[subset]) %>% plot())
  var %>% map( ~as.mcmc(.x) %>% mean())
  
}

range = 1000000
thinned = 1000

myplot("Data/beta.csv", range, thinned)
myplot("Data/gamma.csv", range, thinned)
myplot("Data/alpha.csv", range, thinned)
myplot("Data/theta.csv", range, thinned)
myplot("Data/ptheta.csv", range, thinned)
myplot("Data/pgamma.csv", range, thinned)
myplot("Data/eta.csv", range, thinned)
myplot("Data/tau.csv", range, thinned)
myplot("Data/omega.csv", range, thinned)



read_csv("Data/TimeTable.csv", col_names = FALSE) %>% group_by(2) %>% summarise(Total = n(1))



