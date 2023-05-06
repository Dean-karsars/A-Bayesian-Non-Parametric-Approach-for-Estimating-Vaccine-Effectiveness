
# function def ----------------------------------------------------------------

library(tidytable)
parms_table = function(path, name, fn, ...){
  
  Table = read_csv(path, col_names = FALSE)
  name_vec = paste(name, seq(1, Table %>% ncol()), sep = "")
  names(Table) = name_vec
  fn_name = paste(substitute(fn) %>% deparse(), as.character(c(...)), sep = "")
  
  Inference_Table = Table %>% 
    slice_tail(0.1 * nrow(Table)) %>% 
    summarize(across(everything(), ~ fn(.x, ...))) %>%
    mutate(type = fn_name) %>%
    relocate(type, .before = c(1))
  
  Inference_Table %>% pivot_longer(cols = -type, values_to = fn_name, names_to = "parms") %>% select(-c("type"))
  
}

# beta table ----------------------------------------------------------------

beta_table = 
  full_join(
  parms_table("Data/beta.csv", "beta", mean), 
  parms_table("Data/beta.csv", "beta", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/beta.csv", "beta", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/beta.csv", "beta", quantile, 0.95), 
    by = "parms"
  )
  
# theta table ----------------------------------------------------------------

theta_table =
  full_join(
  parms_table("Data/theta.csv", "theta", mean), 
  parms_table("Data/theta.csv", "theta", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/theta.csv", "theta", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/theta.csv", "theta", quantile, 0.95), 
    by = "parms"
  )


# ptheta table ----------------------------------------------------------------
ptheta_table =
  full_join(
  parms_table("Data/ptheta.csv", "ptheta", mean), 
  parms_table("Data/ptheta.csv", "ptheta", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/ptheta.csv", "ptheta", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/ptheta.csv", "ptheta", quantile, 0.95), 
    by = "parms"
  )

# alpha table ----------------------------------------------------------------
alpha_table =
  full_join(
  parms_table("Data/alpha.csv", "alpha", mean), 
  parms_table("Data/alpha.csv", "alpha", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/alpha.csv", "alpha", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/alpha.csv", "alpha", quantile, 0.95), 
    by = "parms"
  )
# gamma table ----------------------------------------------------------------
gamma_table =
  full_join(
  parms_table("Data/gamma.csv", "gamma", mean), 
  parms_table("Data/gamma.csv", "gamma", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/gamma.csv", "gamma", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/gamma.csv", "gamma", quantile, 0.95), 
    by = "parms"
  )

# pgamma table ----------------------------------------------------------------
pgamma_table =
  full_join(
  parms_table("Data/pgamma.csv", "pgamma", mean), 
  parms_table("Data/pgamma.csv", "pgamma", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/pgamma.csv", "pgamma", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/pgamma.csv", "pgamma", quantile, 0.95), 
    by = "parms"
  )
# omega table ----------------------------------------------------------------
omega_table =
  full_join(
  parms_table("Data/omega.csv", "omega", mean), 
  parms_table("Data/omega.csv", "omega", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/omega.csv", "omega", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/omega.csv", "omega", quantile, 0.95), 
    by = "parms"
  )
# tau table ----------------------------------------------------------------
tau_table =
  full_join(
  parms_table("Data/tau.csv", "tau", mean), 
  parms_table("Data/tau.csv", "tau", sd), 
  by = "parms"
  ) %>%
  full_join(
    parms_table("Data/tau.csv", "tau", quantile, 0.05), 
    by = "parms"
  ) %>%
  full_join(
    parms_table("Data/tau.csv", "tau", quantile, 0.95), 
    by = "parms"
  )


All = beta_table %>% 
  bind_rows(theta_table, 
            ptheta_table, 
            gamma_table, 
            pgamma_table, 
            omega_table, 
            alpha_table, 
            tau_table)
