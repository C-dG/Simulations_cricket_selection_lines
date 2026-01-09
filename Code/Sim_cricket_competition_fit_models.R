#===============================================================================
# Fit models
#===============================================================================

# Save local package information:
#install.packages("renv")
#renv::init()

# Load local package information:
renv::restore()

# Load packages
library(parallel)
library(future)
library(future.apply)
library(rstan)

# Compile Stan model for faster fitting
stan_mod1 <- rstan::stan_model(file = "Models/Comp_animal_mod_2treatment.stan") 
stan_mod2 <- rstan::stan_model(file = "Models/Comp_animal_mod_3treatment.stan") 

# Get specs of datasets 
repl_levels <- seq(5,25, by = 5)
rn_rep_levels <- c(1,2,3)
trt_levels <- c(2,3)

specs <- as.data.frame(expand.grid(repl_levels, rn_rep_levels, trt_levels))
colnames(specs) <- c("n_rep", "n_rn_rep", "n_trt")

specs_trt2 <- specs[specs$n_trt == 2,]
specs_trt3 <- specs[specs$n_trt == 3,]
  
#===============================================================================
# Define stan_func

### Treatment with 2 levels ###
stan_func_mod1 <- function(data) {
  library(rstan)
  
  md <- sampling(stan_mod1, data = data,
                 pars = c("beta", "sigma2_G_mu_trt1", "sigma2_G_mu_trt2",
                          "h2_mu_trt1", "h2_mu_trt2", "Omega_G"),
                 chains = 1, iter = 5000, warmup = 2000, 
                 thin = 1, cores = 1, save_warmup = FALSE)
  
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  dat$n_div <- rep(rstan::get_num_divergent(md), nrow(dat)) # Get n divergencies
  dat$fit_time <- rep(sum(rstan::get_elapsed_time(md)), nrow(dat)) # Get elapsed time of model fit 

  return(dat)
}

### Treatment with 3 levels ### 
stan_func_mod2 <- function(data) {
  library(rstan)
  
  md <- sampling(stan_mod2, data = data,
                 pars = c("beta_back", "sigma2_G_mu", "sigma2_G_beta", 
                          "h2_mu", "h2_beta","Omega_G"),
                 chains = 1, iter = 5000, warmup = 2000, 
                 thin = 1, cores = 1, save_warmup = FALSE)
  
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  dat$n_div <- rep(rstan::get_num_divergent(md), nrow(dat)) # Get n divergencies
  dat$fit_time <- rep(sum(rstan::get_elapsed_time(md)), nrow(dat)) # Get elapsed time of model fit 
  
  return(dat)
}

#===============================================================================

# Fit the models

# Run stan models for simulated datalist 

# Specify parallel worker (cores)
plan(multisession, workers = parallelly::availableCores()-1)

### 2 treatment datasets
lapply(seq_len(nrow(specs_trt2)), function(j) {
  stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt2_rn_rep%d.RDS",
                               specs_trt3[j, "n_rep"], specs_trt3[j, "n_rn_rep"]))
  
  output_list_mod_trt2 <- future_lapply(
    stan_data,
    stan_func_mod1,
    future.globals = list(stan_mod1 = stan_mod1),
    future.seed = TRUE
  )
  
  output_mod_trt2 <- do.call(rbind, output_list_mod_trt2)
  write.csv(output_mod_trt2,
            sprintf("Output/output_nrep%d_trt2_rn_rep%d.csv",
                    specs_trt3[j, "n_rep"], specs_trt3[j, "n_rn_rep"]),
            row.names = FALSE)
  gc()
})


### 3 treatment datasets
lapply(seq_len(nrow(specs_trt3)), function(j) {
  stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                               specs_trt3[j, "n_rep"], specs_trt3[j, "n_rn_rep"]))
  
  output_list_mod_trt3 <- future_lapply(
    stan_data,
    stan_func_mod2,
    future.globals = list(stan_mod2 = stan_mod2),
    future.seed = TRUE
  )
  
  output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
  write.csv(output_mod_trt3,
            sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                    specs_trt3[j, "n_rep"], specs_trt3[j, "n_rn_rep"]),
            row.names = FALSE)
  gc()
})

#===============================================================================

# Treatment 3

# Copy pasted per "iteration"

# Specify parallel worker (cores)
parallelly::availableCores()-2
plan(multisession, workers = 5)

# j == 1 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                               specs_trt3[1, "n_rep"], specs_trt3[1, "n_rn_rep"]))
  
output_list_mod_trt3 <- future_lapply(
    stan_data,
    stan_func_mod2,
    future.globals = list(stan_mod2 = stan_mod2),
    future.seed = TRUE)
  
output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
  write.csv(output_mod_trt3,
            sprintf("Output/output_nrep%d_trt3_rn_rep%d_test.csv",
                    specs_trt3[1, "n_rep"], specs_trt3[1, "n_rn_rep"]),
            row.names = FALSE)
gc()

# j == 2 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[2, "n_rep"], specs_trt3[2, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[2, "n_rep"], specs_trt3[2, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 3 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[3, "n_rep"], specs_trt3[3, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[3, "n_rep"], specs_trt3[3, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 4 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[4, "n_rep"], specs_trt3[4, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[4, "n_rep"], specs_trt3[4, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 5 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[5, "n_rep"], specs_trt3[5, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[5, "n_rep"], specs_trt3[5, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 6 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[6, "n_rep"], specs_trt3[6, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[6, "n_rep"], specs_trt3[6, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 7 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[7, "n_rep"], specs_trt3[7, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[7, "n_rep"], specs_trt3[7, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 8 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[8, "n_rep"], specs_trt3[8, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[8, "n_rep"], specs_trt3[8, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 9 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[9, "n_rep"], specs_trt3[9, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[9, "n_rep"], specs_trt3[9, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 10 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[10, "n_rep"], specs_trt3[10, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[10, "n_rep"], specs_trt3[10, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 11 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[11, "n_rep"], specs_trt3[11, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[11, "n_rep"], specs_trt3[11, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 12 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[12, "n_rep"], specs_trt3[12, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[12, "n_rep"], specs_trt3[12, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 13

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[13, "n_rep"], specs_trt3[13, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[13, "n_rep"], specs_trt3[13, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 14 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[14, "n_rep"], specs_trt3[14, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[14, "n_rep"], specs_trt3[14, "n_rn_rep"]),
          row.names = FALSE)
gc()

# j == 15 

stan_data <- readRDS(sprintf("Data/stan_data_nrep%d_trt3_rn_rep%d.RDS",
                             specs_trt3[15, "n_rep"], specs_trt3[15, "n_rn_rep"]))

output_list_mod_trt3 <- future_lapply(
  stan_data,
  stan_func_mod2,
  future.globals = list(stan_mod2 = stan_mod2),
  future.seed = TRUE)

output_mod_trt3 <- do.call(rbind, output_list_mod_trt3)
write.csv(output_mod_trt3,
          sprintf("Output/output_nrep%d_trt3_rn_rep%d.csv",
                  specs_trt3[15, "n_rep"], specs_trt3[15, "n_rn_rep"]),
          row.names = FALSE)
gc()


#===============================================================================
