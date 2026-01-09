#===============================================================================
# Simulations competition treatment half-sib full-sib cricket lines 
#===============================================================================

# Load in packages
library(nadiv)
library(pedtools)

library(mvnfast)
library(Matrix)
library(matrixcalc)

library(rstan)
library(brms)

#===============================================================================
# Set parameters for simulation structure

# variable parameters
n_rep <- 5     # Number of replications
n_trt <- 3      # Number of levels of competitive treatment 
n_rn_rep <- 1   # Number of reaction norm repeats per dam

# Should not be changed
n_pair <- 6                 # Number of parental pairs in replicate
n_dam_f2 <- 2               # Number of f2 dams used

#===============================================================================
# Set parameters for simulation values

K <- 2 # Number of I-level traits (intercept & slope)
# Total number of observations (unique F3 females)
N <- n_rep * (n_trt * n_rn_rep) * (n_pair * n_dam_f2)

# Variance of treatment 
var_trt <- var(rep(1:n_trt, N))

# Genetic/breeding values
var_g_mu <- 0.2
var_g_beta <- 0.2 * var_trt
r_g <- 0.5*(sqrt(var_g_mu * var_g_beta)) # genetic correlation slope & intercept

# PE values 
var_PE_mu <- 0.2
var_PE_beta <- 0.2 * var_trt
r_PE <- 0.5*(sqrt(var_PE_mu * var_PE_beta)) # PE correlation slope & intercept

# population-level parameters
mu <- 0
beta <- 1
#var_e <- 1-(var_g_mu + var_g_beta + var_PE_mu + var_PE_beta)
var_e <- 1-(var_g_mu + var_g_beta)

#===============================================================================
# Create structure at the end of the breeding design (F3)

# Number of offspring used in experiments
n_off_rep <- (n_trt * n_rn_rep) * (n_pair * n_dam_f2)   

# Number of dams within replicate
n_dam_rep <- n_pair * n_dam_f2

# Number of offspring per dam within a replicate
n_off_dam <- n_trt * n_rn_rep

df <- data.frame(#ID = 1:N,
                 ID = rep(1:n_off_rep, times = n_rep),
                 Dam_ID = rep(1:n_dam_rep, each = n_off_dam, times = n_rep),
                 #Sire_ID = rep(1:n_pair, each = n_pair, times = n_rep), #1 male per 
                 Sire_ID = rep(1:n_dam_rep, each = n_off_dam, times = n_rep), # 2 males per
                 #Sire_ID = rep(rep(1:n_pair, each = n_off_dam, times = n_dam_f2), times = n_rep), #2 male but unique per sister 
                 Rep_ID = rep(rep(1:n_rep, each = n_off_rep)),
                 Treatment = rep(rep(1:n_trt, times = n_dam_rep), times = n_rep)
                 )

#===============================================================================
# Generate within replicate pedigree structure

# P (wild-born) 3 females, 3 males (parentage unknown)
df_p <- data.frame(
  ID = paste0(1:6, sep = "_", "P"),
  Dam_ID = rep(NA, 6),
  Sire_ID = rep(NA, 6))

# F1
df_f1 <- data.frame(
  ID = paste0(1:9, sep = "_", "f1"),
  Dam_ID = rep(df_p[1:3, "ID"], each = 3),
  Sire_ID = rep(df_p[4:6, "ID"], each = 3))

# F3
df_f3 <- df[df$Rep_ID == 1, c("ID", "Dam_ID", "Sire_ID")]

df_f3$ID <- paste0(df_f3$ID, sep = "_", "F3")
df_f3$Dam_ID <- paste0(df_f3$Dam_ID, sep = "_", "F2_D")
df_f3$Sire_ID <- paste0(df_f3$Sire_ID, sep = "_", "F2_S")

# F2
# Sire ID's are assumed to be unknown (unrelated to replicate)

# F2 Sires have 2 F2 Dam partners (but not sisters)
# df_f2 <- data.frame(
#   ID = c(unique(df_f3$Dam_ID), unique(df_f3$Sire_ID)),
#   Dam_ID = c(rep(df_f1[c(4,7,1,8,2,5),"ID"], each = 2),
#              rep(NA,6)), # F3 Sires are unrelated
#   Sire_ID = c(rep(df_f1[c(3,6,9),"ID"], each = 4),
#               rep(NA,6))) # F3 Sires are unrelated

# F2 Sisters have a unique male
df_f2 <- data.frame(
  ID = c(unique(df_f3$Dam_ID), unique(df_f3$Sire_ID)),
  Dam_ID = c(rep(df_f1[c(4,7,1,8,2,5),"ID"], each = 2),
             rep(NA,12)), # F3 Sires are unrelated
  Sire_ID = c(rep(df_f1[c(3,6,9),"ID"], each = 4),
              rep(NA,12))) # F3 Sires are unrelated

# Combine into one pedigree 
ped <- rbind(df_f3, df_f2, df_f1, df_p)

# 1: M, 2: F (...)
ped$Sex <- ifelse(ped$ID %in% ped$Dam_ID, 2,
           ifelse(ped$ID %in% ped$Sire_ID, 1, 2))

# Create pedigree image
plot(pedtools::ped(ped$ID, ped$Sire_ID, ped$Dam_ID, ped$Sex))

#===============================================================================
# Create A matrices (relatedness) 

### Based on full pedigree
A <- as.matrix(nadiv::makeA(prepPed(rbind(df_f3, df_f2, df_f1, df_p))))

table(A)
table(A)/(nrow(A)*nrow(A))

# Replicate for all replicates 
A_rep_full <- rep(list(A), n_rep) # Makes a list of all matrices
A_full <- as.matrix(bdiag(A_rep_full)) #Binds matrices into one

# Rename rows and columns to match ID's in df
row.names(A_full) <- paste0(row.names(A),"_", rep(1:n_rep, each = nrow(A)))
colnames(A_full) <- paste0(row.names(A),"_", rep(1:n_rep, each = nrow(A)))

#### Subset for F3 focal individuals (is in correct order!)
A_rep <- A[df_f3$ID, df_f3$ID]

# Table of pairwise relatedness classes within replicate
table(A_rep)
table(A_rep)/(nrow(A_rep)*nrow(A_rep))

# Replicate for all replicates 
A_rep_list <- rep(list(A_rep), n_rep) # Makes a list of all matrices
A_F3 <- as.matrix(bdiag(A_rep_list)) #Binds matrices into one

### Maternal pedigree 
A_rep_dam <- A[unique(df_f3$Dam_ID), unique(df_f3$Dam_ID)]

# Table of pairwise relatedness classes within replicate
table(A_rep_dam)
table(A_rep_dam)/(nrow(A_rep_dam)*nrow(A_rep_dam))

# Replicate for all replicates 
A_rep_dam_list <- rep(list(A_rep_dam), n_rep) # Makes a list of all matrices
A_F2 <- as.matrix(bdiag(A_rep_dam_list)) #Binds matrices into one

# Rename rows and columns to match focal ID's in pedigree
row.names(A_F2) <- 1:nrow(A_F2)
colnames(A_F2) <- 1:nrow(A_F2)

# Check whether matrices are positive definite
matrixcalc::is.positive.definite(A_full)
matrixcalc::is.positive.definite(A_F3)
matrixcalc::is.positive.definite(A_F2)

#Visualise the relatedness matrix
ggpedigree::ggRelatednessMatrix(
  A,
  interactive = FALSE,
  config = list(
    color_scale_midpoint = 0.50,
    plot_title = "",
    axis_text_size = 12
  ))
# 
ggpedigree::ggRelatednessMatrix(
  A_rep,
  interactive = FALSE,
  config = list(
    color_scale_midpoint = 0.50,
    plot_title = "",
    axis_text_size = 12
  ))
#
ggpedigree::ggRelatednessMatrix(
  A_F3,
  interactive = FALSE,
  config = list(
    color_scale_midpoint = 0.50,
    plot_title = "",
    axis_text_size = 12
  ))
# 
# ggpedigree::ggRelatednessMatrix(
#   A_F2,
#   interactive = FALSE,
#   config = list(
#     color_scale_midpoint = 0.50,
#     plot_title = "",
#     axis_text_size = 12
#   ))

#===============================================================================
# Create unique individual indices in df (linked to pedigree)
df$ID <- paste0(df$ID, "_F3_", df$Rep_ID)
df$Dam_ID <- paste0(df$Dam_ID, "_F2_D_", df$Rep_ID)
df$Sire_ID <- paste0(df$Sire_ID, "_F2_S_", df$Rep_ID)

# Center treatment
df$Treatment <- as.vector(scale(as.numeric(df$Treatment), scale = FALSE))

# Rename rows and columns to match focal ID's in pedigree
row.names(A_F3) <- paste0(row.names(A_rep),"_", rep(1:n_rep, each = nrow(A_rep)))
colnames(A_F3) <- paste0(row.names(A_rep),"_", rep(1:n_rep, each = nrow(A_rep)))

table(df$ID == row.names(A_F3)) # Check whether ID match with df


# Same for F2
row.names(A_F2) <- paste0(row.names(A_rep_dam),"_", rep(1:n_rep, each = nrow(A_rep_dam)))
colnames(A_F2) <- paste0(row.names(A_rep_dam),"_", rep(1:n_rep, each = nrow(A_rep_dam)))

table(unique(df$Dam_ID) == row.names(A_F2)) # Check whether ID match with df

#===============================================================================
# Simulate values
G_cov <- matrix(c(var_g_mu, r_g,
              r_g, var_g_beta), 
              nrow = K, ncol = K)  

# # Scale G-matrix with A-matrix (F3 pedigree)
# Sigma_G <- G_cov %x% A_full # Kronecker product of G & A
# 
# # Simulate breeding values
# BV <- rmvn(n = 1, mu = rep(0,K * nrow(A_full)), sigma = Sigma_G)
# 
# # Reshape to N × K matrix
# BV <- as.data.frame(matrix(BV, nrow = nrow(A_full), ncol = K))
# 
# colnames(BV) <- c("BV_mu", "BV_beta")
# cor(BV) # Check whether the correlation matches set values
# 
# # Subset BV's for individuals that will be phenotyped
# BV$ID <- row.names(A_full)
# BV_F3 <- BV[BV$ID %in% df$ID,]
# df <- merge(df, BV_F3, by = "ID") # Merge into df


# Scale G-matrix with A-matrix (F3 pedigree)
Sigma_G <- G_cov %x% A_F3 # Kronecker product of G & A

# Simulate breeding values
BV <- rmvn(n = 1, mu = rep(0,K * nrow(A_F3)), sigma = Sigma_G)

# Reshape to N × K matrix
BV <- as.data.frame(matrix(BV, nrow = nrow(A_F3), ncol = K))

colnames(BV) <- c("BV_mu", "BV_beta")
cor(BV) # Check whether the correlation matches set values

# Merge Breeding values into df
BV$ID <- row.names(A_F3)
df <- merge(df, BV, by = "ID")

# # Simulate PE values 
# # Based on Dam_ID, because offspring does not have repeated measures
# Sigma_PE <- matrix(c(var_PE_mu, r_PE,
#                      r_PE, var_PE_beta),
#                    nrow = K, ncol = K)
# 
# PE <- as.data.frame(rmvn(n = length(unique(df$Dam_ID))
#                          , mu = rep(0,K), sigma = Sigma_PE))
# colnames(PE) <- c("PE_mu", "PE_beta")
# cor(PE) # Check whether the correlation matches set values
# 
# # Match with Dam_ID in df
# PE$Dam_ID <- unique(df$Dam_ID)
# df <- merge(df, PE, by = "Dam_ID")

Sigma_PE <- matrix(c(var_PE_mu, r_PE,
                     r_PE, var_PE_beta),
                   nrow = K, ncol = K)

PE <- as.data.frame(rmvn(n = length(unique(df$ID))
                         , mu = rep(0,K), sigma = Sigma_PE))
colnames(PE) <- c("PE_mu", "PE_beta")
cor(PE) # Check whether the correlation matches set values

# Match with Dam_ID in df
PE$ID <- row.names(A_F3)
df <- merge(df, PE, by = "ID")


# df$PE_mu <- unlist(PE[match(df$Dam_ID, PE$ID),"PE_mu"])
# df$PE_beta <- unlist(PE[match(df$Dam_ID, PE$ID),"PE_beta"])

# Residual variance
df$Res <- rnorm(nrow(df), 0, sqrt(var_e)) 

# Phenotypic equation
df$z <- (mu + df$BV_mu + df$PE_mu) +
  ((beta + df$BV_beta + df$PE_beta)*df$Treatment) +
  df$Res

# Without PE
df$z <- (mu + df$BV_mu) +
  ((beta + df$BV_beta)*df$Treatment) +
  df$Res

hist(df$z)
var(df$z)
str(df)

# Correlations should be ~0.5 
cor(df$BV_mu, df$BV_beta)
cor(df$PE_mu, df$PE_beta)

var(df$PE_beta)

# Order df ID's  according to A
df <- df[match(row.names(A_F3), df$ID),] 
table(df$ID == row.names(A_F3))

# Give focal's ID that matches with row and column ID of A 
df$ID_index <- as.integer(factor(df$ID, levels = row.names(A_F3)))

# Same for dams
df$Dam_index <- as.integer(factor(df$Dam_ID, levels = row.names(A_F2)))

#===============================================================
# Transform data to have 6 repeats per individual
nrep_df <- 6
  
df2 <- do.call("rbind", replicate(nrep_df, df, simplify = FALSE))

# Re-assign error and treatment
df2$Res <- rnorm(nrow(df2), 0, sqrt(var_e)) 
df2$Treatment <- rep(unique(df$Treatment), each = length(unique(df$ID)), times = nrep_df/n_trt)

# Recalculate phenotypes
df2$z <- (mu + df2$BV_mu + df2$PE_mu) +
  ((beta + df2$BV_beta + df2$PE_beta)*df2$Treatment) +
  df2$Res

# Correlations should be ~0.5 
cor(df2$BV_mu, df2$BV_beta)
cor(df2$PE_mu, df2$PE_beta)

#===============================================================================
# Model

# BRMS

# Only F3 pedigree
md_brms <- brm(
  formula = z ~ Treatment + (Treatment|Dam_ID) + (Treatment|gr(ID, cov = A)),
  data   = df,
  data2  = list(A = A_F3),
  family = gaussian(),
  chains = 3, cores = 6,
  warmup = 2000, iter = 5000
)

md_brms <- brm(
  formula = z ~ Treatment + (Treatment|gr(ID, cov = A)),
  data   = df,
  data2  = list(A = A_F3),
  family = gaussian(),
  chains = 1, cores = 4,
  warmup = 2000, iter = 5000
)

md_brms <- brm(
  formula = z ~ 1 + (1|gr(ID, cov = A)),
  data   = df,
  data2  = list(A = A_F3),
  family = gaussian(),
  chains = 1, cores = 4,
  warmup = 2000, iter = 5000
)

# Full pedigree
md_brms <- brm(
  formula = z ~ Treatment + (Treatment|Dam_ID) + (Treatment|gr(ID, cov = A)),
  data   = df,
  data2  = list(A = A_full),
  family = gaussian(),
  chains = 3, cores = 6,
  warmup = 2000, iter = 5000
)

summary(md_brms)

# Stan
# Estimate breeding values for intercept & slope + correlation
# Estimate brood-level permanent environment effects

# ID G, Dam PE
stan_data <- list(No = nrow(df),
                  z = as.vector(scale(df$z)),
                  animal = df$ID_index,
                  Na = length(unique(df$ID)),
                  dam = as.integer(as.factor(df$Dam_ID)),
                  Ndam = length(unique(df$Dam_ID)),
                  Trt = df$Treatment,
                  A = A_F3)

# ID G, No PE
stan_data <- list(No = nrow(df),
                  z = as.vector(scale(df$z)),
                  sd_z = sd(df$z),
                  animal = df$ID_index,
                  Na = length(unique(df$ID)),
                  Trt = df$Treatment,
                  A = A_F3)

# With repeats G ID, PE ID
stan_data <- list(No = nrow(df2),
                  z = as.vector(scale(df2$z)),
                  animal = df2$ID_index,
                  Na = length(unique(df2$ID)),
                  #dam = as.integer(as.factor(df2$Dam_ID)),
                  #Ndam = length(unique(df2$Dam_ID)),
                  Trt = df2$Treatment,
                  A = A_F3)

# Dam level G, Dam PE 
stan_data <- list(No = nrow(df),
                  z = as.vector(scale(df$z)),
                  animal = df$Dam_index,
                  Na = length(unique(df$Dam_ID)),
                  #dam = as.integer(as.factor(df$Dam_ID)),
                  #Ndam = length(unique(df$Dam_ID)),
                  Trt = df$Treatment,
                  A = A_F2)

# 3 treatment level
write(temp <- "
data {
// Number of clusters
int<lower=1> No;  // total number of observations
int<lower=1> Na; // number of individuals
int<lower=1> Ndam; // number of dams

// Response variable
vector[No] z; // Phenotypic observations

// Fixed effects
vector[No] Trt ; // Treatment with n-levels

// Random effects
array[No] int<lower=1> animal; // individual identity for each observation
array[No] int<lower=1> dam; // dam identity for each observation
matrix[Na,Na] A; // Relatedness matrix
}

transformed data {
matrix[Na,Na] LA = cholesky_decompose(A); // Decompose A for easier fitting     
}

parameters {
// Population-level effects
real mu; // Population intercept
real beta; // Population slope

// Random effects
vector<lower=0>[2] sd_I; // Individual-level standard deviation
matrix[Ndam,2] Iz; // Standardised individual-level effects

vector<lower=0>[2] sd_G; // Breeding values standard deviation
matrix[Na,2] Gz; // Standardised breeding values

real<lower=0> sigma_e; // Residual standard deviation

cholesky_factor_corr[2] LI; // PE correlation matrix
cholesky_factor_corr[2] LG; // G correlation matrix
}

transformed parameters {
matrix[Ndam,2] I = Iz * diag_pre_multiply(sd_I, LI)'; // I/PE matrix (PE values)
matrix[Na,2] G = LA * Gz * diag_pre_multiply(sd_G, LG)'; // G matrix (breeding values)
}

model {
vector[No] e_z; // Predicted model values

// Model equation
e_z = mu + I[dam,1] + G[animal,1] + (beta + I[dam,2] + G[animal,2]) .* Trt;

// Priors
mu ~ normal(0,1); // Weakly informative prior, implies i ~ normal(0,sd_Z)
beta ~ normal(0,1);

to_vector(Iz) ~ normal(0,1);
to_vector(Gz) ~ normal(0,1); 

sd_I ~ exponential(2);
sd_G ~ exponential(2);
sigma_e ~ exponential(2);

LI ~ lkj_corr_cholesky(1.5); // Prior PE correlations
LG ~ lkj_corr_cholesky(1.5); // Prior G correlations

// Model likelihood
z ~ normal(e_z, sigma_e);
}

generated quantities {
// Variances
real sigma2_PE_mu = sd_I[1]^2;
real sigma2_PE_beta = sd_I[2]^2;
real sigma2_G_mu = sd_G[1]^2;
real sigma2_G_beta = sd_G[2]^2;

real sigma2_e = sigma_e^2;
real sigma2_P = sigma2_PE_mu + sigma2_PE_beta + 
                sigma2_G_mu + sigma2_G_beta + sigma2_e;

// Heritability 
real h2_mu  = sigma2_G_mu/ sigma2_P;
real h2_beta  = sigma2_G_beta/ sigma2_P; 

// Repeatability
real R_mu  = sigma2_PE_mu/ sigma2_P;
real R_beta  = sigma2_PE_beta/ sigma2_P; 

// Correlations 
matrix[2,2] Omega_I = LI * LI';
matrix[2,2] Omega_G = LG * LG';

// Covariances
matrix[2,2] D_I = diag_matrix(sd_I);
matrix[2,2] D_G = diag_matrix(sd_G);

matrix[2,2] Sigma_I = D_I * Omega_I * D_I;
matrix[2,2] Sigma_G = D_G * Omega_G * D_G;
}",  
file = "Comp_animal_mod.stan")



# Only G estimation, focal level
write(temp <- "
data {
// Number of clusters
int<lower=1> No;  // total number of observations
int<lower=1> Na; // number of individuals

// Response variable
vector[No] z; // Phenotypic observations

// Fixed effects
vector[No] Trt ; // Treatment with n-levels

// Random effects
array[No] int<lower=1> animal; // individual identity for each observation
matrix[Na,Na] A; // Relatedness matrix
}

transformed data {
matrix[Na,Na] LA = cholesky_decompose(A); // Decompose A for easier fitting     
}

parameters {
// Population-level effects
real mu; // Population intercept
real beta; // Population slope

// Random effects
vector<lower=0>[2] sd_G; // Breeding values standard deviation
matrix[Na,2] Gz; // Standardised breeding values

real<lower=0> sigma_e; // Residual standard deviation

cholesky_factor_corr[2] LG; // G correlation matrix
}

transformed parameters {
matrix[Na,2] G = LA * Gz * diag_pre_multiply(sd_G, LG)'; // G matrix (breeding values)
}

model {
vector[No] e_z; // Predicted model values

// Model equation
e_z = mu + G[animal,1] + (beta + G[animal,2]) .* Trt;

// Priors
mu ~ normal(0,1); // Weakly informative prior, implies i ~ normal(0,sd_Z)
beta ~ normal(0,1);

to_vector(Gz) ~ normal(0,1); 
sd_G ~ exponential(2);

sigma_e ~ exponential(2);

LG ~ lkj_corr_cholesky(1.5); // Prior G correlations

// Model likelihood
z ~ normal(e_z, sigma_e);
}

generated quantities {
// Variances
real sigma2_G_mu = sd_G[1]^2;
real sigma2_G_beta = sd_G[2]^2;

real sigma2_e = sigma_e^2;
real sigma2_P = sigma2_G_mu + sigma2_G_beta + sigma2_e;

// Heritability 
real h2_mu  = sigma2_G_mu/ sigma2_P;
real h2_beta  = sigma2_G_beta/ sigma2_P; 

// Correlations 
matrix[2,2] Omega_G = LG * LG';

// Covariances
matrix[2,2] D_G = diag_matrix(sd_G);
matrix[2,2] Sigma_G = D_G * Omega_G * D_G;
}",  
file = "Comp_animal_mod.stan")



# PE at offpsring level (repetaed measures)
write(temp <- "
data {
// Number of clusters
int<lower=1> No;  // total number of observations
int<lower=1> Na; // number of individuals

// Response variable
vector[No] z; // Phenotypic observations

// Fixed effects
vector[No] Trt ; // Treatment with n-levels

// Random effects
array[No] int<lower=1> animal; // individual identity for each observation
matrix[Na,Na] A; // Relatedness matrix
}

transformed data {
matrix[Na,Na] LA = cholesky_decompose(A); // Decompose A for easier fitting     
}

parameters {
// Population-level effects
real mu; // Population intercept
real beta; // Population slope

// Random effects
matrix[Na,2] PEz; // Standardised individual-level effects
matrix[Na,2] Gz; // Standardised breeding values

real<lower=0> sigma_e; // Residual standard deviation

vector<lower=0>[2] sd_I; // Standard deviation individual-level 

vector<lower=0, upper =1>[2] h2; // varG/ varP

cholesky_factor_corr[2] LPE; // PE correlation matrix
cholesky_factor_corr[2] LG; // G correlation matrix
}

transformed parameters {
vector<lower=0>[2] sd_PE; // Standard deviation PE
vector<lower=0>[2] sd_G; // Standard deviation breeding values

sd_G = sd_I .* sqrt(h2);
sd_PE = sd_I .* sqrt(1 - h2);

matrix[Na,2] PE = PEz * diag_pre_multiply(sd_PE, LPE)'; // PE matrix 
matrix[Na,2] G = LA * Gz * diag_pre_multiply(sd_G, LG)'; // G matrix (breeding values)
matrix[Na,2] I = G + PE; // Total individual-level effects
}

model {
vector[No] e_z; // Predicted model values

// Model equation
e_z = mu + I[animal,1] + (beta + I[animal,2]) .* Trt;

// Priors
mu ~ normal(0,1); // Weakly informative prior, implies i ~ normal(0,sd_Z)
beta ~ normal(0,1);

to_vector(PEz) ~ normal(0,1);
to_vector(Gz) ~ normal(0,1); 

sd_I ~ exponential(2);
//sd_PE ~ exponential(2);
//sd_G ~ exponential(2);
sigma_e ~ exponential(2);

to_vector(h2) ~ beta(1.2,1.2);

LPE ~ lkj_corr_cholesky(1.5); // Prior PE correlations
LG ~ lkj_corr_cholesky(1.5); // Prior G correlations

// Model likelihood
z ~ normal(e_z, sigma_e);
}

generated quantities {
// Variances
real sigma2_PE_mu = sd_PE[1]^2;
real sigma2_PE_beta = sd_PE[2]^2;
real sigma2_G_mu = sd_G[1]^2;
real sigma2_G_beta = sd_G[2]^2;

real sigma2_e = sigma_e^2;
real sigma2_P = sigma2_PE_mu + sigma2_PE_beta + 
                sigma2_G_mu + sigma2_G_beta + sigma2_e;

// Heritability 
real h2_mu  = sigma2_G_mu/ sigma2_P;
real h2_beta  = sigma2_G_beta/ sigma2_P; 

// Repeatability
real R_mu  = sigma2_PE_mu/ sigma2_P;
real R_beta  = sigma2_PE_beta/ sigma2_P; 

// Correlations 
matrix[2,2] Omega_PE = LPE * LPE';
matrix[2,2] Omega_G = LG * LG';

// Covariances
matrix[2,2] D_I = diag_matrix(sd_I);
matrix[2,2] D_G = diag_matrix(sd_G);

matrix[2,2] Sigma_I = D_I * Omega_PE * D_I;
matrix[2,2] Sigma_G = D_G * Omega_G * D_G;
}",  
file = "Comp_animal_mod.stan")

# Compile Stan model for faster fitting
stan_animal_mod <- rstan::stan_model(file = "Comp_animal_mod.stan") 

gc()

# Fit the stan model
stan_fit <- rstan::sampling(stan_animal_mod, 
            data = stan_data, 
            chains = 1, 
            warmup = 2000, iter = 5000, 
            thin = 1, 
            cores = parallelly::availableCores()-1,
            save_warmup = FALSE,
)

rstan::get_elapsed_time(stan_fit)

pars <- c("mu", "beta", 
          "sigma2_G_mu", "sigma2_G_beta", 
          "sigma2_PE_mu", "sigma2_PE_beta",
          "sigma2_e", "sigma2_P",
          "h2_mu", "h2_beta", "R_mu", "R_beta", 
          "Omega_G", "Omega_I"
)

pars <- c("mu", "beta", 
          "sigma2_G_mu", "sigma2_G_beta", 
          "sigma2_e", "sigma2_P",
          "h2_mu", "h2_beta", 
          "Omega_G"
)

pars <- c("mu", "beta", 
          "sigma2_G_mu", "sigma2_G_beta", 
          "sigma2_PE_mu", "sigma2_PE_beta",
          "sigma2_e", "sigma2_P",
          "h2_mu", "h2_beta", "R_mu", "R_beta", 
          "Omega_G", "Omega_PE" , "h2", "sd_I"
          )

round(summary(stan_fit, pars = pars)$summary[,c(1,4,6,8,9,10)],3)

shinystan::launch_shinystan(shinystan::as.shinystan(stan_fit, pars = pars))

# 2 treatment model
# Only G estimation, focal level
write(temp <- "
data {
// Number of clusters
int<lower=1> No;  // total number of observations
int<lower=1> Na; // number of individuals

// Response variable
vector[No] z; // Phenotypic observations
real sd_z; // Standard deviation z

// Fixed effects
vector[No] Trt ; // Treatment with n-levels

// Random effects
array[No] int<lower=1> animal; // individual identity for each observation
matrix[Na,Na] A; // Relatedness matrix
}

transformed data {
matrix[Na,Na] LA = cholesky_decompose(A); // Decompose A for easier fitting

// Split data per treatment level 

}

parameters {
// Population-level effects
vector[2] mu; // Population intercept for treatment 1 and 2

// Random effects
vector<lower=0>[2] sd_G; // Breeding values standard deviation
matrix[Na,2] Gz; // Standardised breeding values

vector<lower=0>[2] sigma_e; // Residual standard deviation

cholesky_factor_corr[2] LG; // G correlation matrix
}

transformed parameters {
matrix[Na,2] G = LA * Gz * diag_pre_multiply(sd_G, LG)'; // G matrix (breeding values)
}

model {
// Model equation

for (i in 1:No) {
  if (Trt[i] == -0.5)
    z[i] ~ normal(mu[1] + G[animal[i], 1], sigma_e[1]);
  else
    z[i] ~ normal(mu[2] + G[animal[i], 2], sigma_e[2]);
}

// Priors
mu ~ normal(0,1); // Weakly informative prior, implies i ~ normal(0,sd_Z)

to_vector(Gz) ~ normal(0,1); 
sd_G ~ exponential(2);

sigma_e ~ exponential(2);

LG ~ lkj_corr_cholesky(1.5); // Prior G correlations
}

generated quantities {
// beta
real beta = (abs(mu[1]) + mu[2])*sd_z;

// Variances
real sigma2_G_mu_trt1 = sd_G[1]^2;
real sigma2_G_mu_trt2 = sd_G[2]^2;

real sigma2_e1 = sigma_e[1]^2;
real sigma2_e2 = sigma_e[2]^2;

// Heritability 
real h2_mu_trt1  = sigma2_G_mu_trt1 / (sigma2_G_mu_trt1 + sigma2_e1);
real h2_mu_trt2  = sigma2_G_mu_trt2 / (sigma2_G_mu_trt2 + sigma2_e2);

// Correlations 
matrix[2,2] Omega_G = LG * LG';

// Covariances
matrix[2,2] D_G = diag_matrix(sd_G);
matrix[2,2] Sigma_G = D_G * Omega_G * D_G;
}",  
file = "Comp_animal_mod_2treatment.stan")

# Compile Stan model for faster fitting
stan_animal_mod_2treatment <- rstan::stan_model(file = "Comp_animal_mod_2treatment.stan") 

gc()

# Fit the stan model
stan_fit <- rstan::sampling(stan_animal_mod_2treatment, 
                            data = stan_data, 
                            chains = 1, 
                            warmup = 2000, iter = 5000, 
                            thin = 1, 
                            cores = parallelly::availableCores()-1,
                            save_warmup = FALSE,
)

pars <- c("mu", "beta",
          "sigma2_G_mu_trt1", "sigma2_G_mu_trt2", 
          "sigma2_e1", "sigma2_e2", 
          "h2_mu_trt1", "h2_mu_trt2", 
          "Omega_G"
)

round(summary(stan_fit, pars = pars)$summary[,c(1,4,6,8,9,10)],3)

#===============================================================================
# Process output

# Relative bias
# Proportion of overlap with 0 (same with 0 distribution for variances)


#===============================================================================


