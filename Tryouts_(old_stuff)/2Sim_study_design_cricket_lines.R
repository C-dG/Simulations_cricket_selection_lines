#===============================================================================
# Simulations competition treatment half-sib full-sib cricket lines 
#===============================================================================

# Load in packages
library(mvnfast)
library(Matrix)
library(stan)
library(brms)

#===============================================================================
# Set parameters for simulation structure

# variable parameters
n_rep <- 25     # Number of replications
n_trt <- 3      # Number of levels of competitive treatment 
n_rn_rep <- 1   # Number of reaction norm repeats per dam

# Should not be changed
n_pair <- 6                 # Number of parental pairs in replicate
n_dam_f2 <- 2               # Number of f2 dams used

# Number of offspring used in experiments
n_off <- (n_trt * n_rn_rep) * (n_pair * n_dam_f2)   

#===============================================================================
# Set parameters for simulation values

K <- 2 # Number of I-level traits (intercept & slope)
N <- n_rep * (n_trt * n_rn_rep) * (n_pair * n_dam_f2)

# Genetic/breeding values
var_g_mu <- 0.2
var_g_beta <- 0.2
r_g <- 0.5*(sqrt(var_g_mu * var_g_beta)) # genetic correlation slope & intercept

# PE values 
var_PE_mu <- 0.2
var_PE_beta <- 0.2
r_PE <- 0.5*(sqrt(var_PE_mu * var_PE_beta)) # PE correlation slope & intercept

# population-level parameters
mu <- 0
beta <- 1
var_e <- 1-(var_g_mu + var_g_beta + var_PE_mu + var_PE_beta)

#===============================================================================
# Create structure at the end of the breeding design (F3)

n_off_rep <- n_pair * n_off

df <- data.frame(ID = rep(1:n_off_rep, times = n_rep),
                 P_pair_ID = rep(1:n_pair, each = n_off, times = n_rep),
                 Rep_ID = rep(rep(1:n_rep, each = n_off_rep)),
                 Treatment =  rep(rep(1:n_trt, times = n_pair), times = n_rep))

df$Treatment <- df$Treatment - 1 # Ensure lowest treatment-level is 0 (no competition)

#===============================================================================
# Construct edge list of relatedness within replicate

A_rep_list <- as.data.frame(expand.grid(1:n_off_rep, 1:n_off_rep))
colnames(A_rep_list) <- c("ID1", "ID2")

df1 <- df[df$Rep_ID == 1, ]

# Match with parental pair ID
A_rep_list$P_ID1 <- df1[match(A_rep_list$ID1, df1$ID),"P_pair_ID"]
A_rep_list$P_ID2 <- df1[match(A_rep_list$ID2, df1$ID),"P_pair_ID"]

# Assign relatedness to combinations of parental ID's
A_rep_list$R <- ifelse(A_rep_list$P_ID1 == A_rep_list$P_ID2, 0.5,
                       ifelse((paste(pmin(A_rep_list$P_ID1, A_rep_list$P_ID2), 
                                     pmax(A_rep_list$P_ID1, A_rep_list$P_ID2), sep = ",") %in%
                                 c("1,2", "3,4", "5,6")),0.25,0.125))

A_rep_list$P_comb <- paste0(A_rep_list$P_ID1, sep = "_", A_rep_list$P_ID2)

# Turn edge list into a matrix
A_rep <- matrix(A_rep_list$R, nrow = n_off_rep, ncol = n_off_rep)
diag(A_rep) <- 1 #Change diagonal to R = 1 (ID's are 100% related to themselves)

# create unique ID's within replicate
#df$ID <- paste0(df$ID, sep = "_", df$Rep_ID)
df$ID <- seq_along(df$ID)
df$P_pair_ID <- paste0(df$P_pair_ID, sep = "_", df$Rep_ID)

# n_id * n_id A matrix
A_rep_list <- rep(list(A_rep), n_rep) # Makes a list of all matrices
A <- as.matrix(bdiag(A_rep_list)) #Binds matrices into one

# label matrix with focal ID
rownames(A) <- df$ID
colnames(A) <- df$ID

#===============================================================================
# Simulate values
G <- matrix(c(var_g_mu, r_g,
              r_g, var_g_beta), nrow = K, ncol = K)  

# Scale G-matrix with A-matrix
Sigma <- kronecker(G, A)

# Simulate breeding values
BV <- rmvn(n = 1, mu = rep(0,K * N), sigma = Sigma)

# Reshape to N × K matrix
BV <- as.data.frame(matrix(BV, nrow = N, ncol = K, byrow = TRUE))
colnames(BV) <- c("BV_mu", "BV_beta")

# Simulate PE values
Sigma_PE <- matrix(c(var_PE_mu, r_PE,
                     r_PE, var_PE_beta), nrow = K, ncol = K)  

PE <- as.data.frame(rmvn(n = N, mu = rep(0,K), sigma = Sigma_PE))
colnames(PE) <- c("PE_mu", "PE_beta")

# Residual variance
Res <- rnorm(nrow(df), 0, var_e) 

# Phenotypic equation
df <- cbind(df, BV, PE, Res)

df$Treatment <- as.numeric(df$Treatment)

df$z <- (mu + df$BV_mu + df$PE_mu) +
  ((beta + df$BV_beta + df$PE_beta)*df$Treatment) +
  df$Res

#===============================================================================
# Model

# BRMS
md_brms <- brm(z ~ Treatment + (Treatment|P_pair_ID), data = df)

md_brms <- brm(
  formula = z ~ Treatment + (Treatment|P_pair_ID) + (Treatment|gr(ID, cov = A)),
  data   = df,
  data2  = list(A = A),
  family = gaussian(),
  chains = 4, cores = 4
)

summary(md_brms)

# Stan
stan_data_test <- list(No = nrow(df),
                       z = df$z,
                       animal = df$ID,
                       Na = unique(df$ID),
                       A = A)

stan_data_test <- list(No = nrow(df),
                       z = df$z_1,
                       animal = df$individual,
                       opponent1 = df$opponent1,
                       opponent2 = df$opponent2,
                       Na = unique(df$n_ind),
                       A = A[(ncol(A) - n.ind +1):nrow(A),(nrow(A) - n.ind + 1):ncol(A)])

write(temp <- "data {
  int<lower=1> No;  // total number of observations
  vector[No] z;  // focal trait (i.e. response variable)
  int animal[No]; // individual identity for each observation
  int<lower=1> opponent1[No]; // individual identity for each observation
  int<lower=1> opponent2[No]; // individual identity for each observation
  int<lower=1> Na; // number of individuals
  cov_matrix[Na] A; // additive genetic relatedness matrix
}
transformed data {
  real mean_z;  // sample mean of observed trait values
  real sd_z;  // sample standard deviation of observed trait values
  matrix[Na,Na] LA; // lower-triangular cholesky factor of A
  mean_z = mean(z);
  sd_z = sd(z);
  LA = cholesky_decompose(A);
}
parameters {
  real mu; // overall intercept
  vector<lower=0>[2] sd_G ; // additive genetic standard deviation
  vector<lower=0>[2] sd_I; // permanent individual standard deviation
  real<lower=0> sd_R; // temporary indvidual standard deviation (i.e. residual variance)
  matrix[Na,2] Gz; // standardised breeding values
  matrix[Na,2] iz; // standardised permanent individual effects
  cholesky_factor_corr[2] LG; // Genetic correlations
  cholesky_factor_corr[2] LI; // Genetic correlations
}
transformed parameters {
  //vector<lower=0>[2] var_G; // additive genetic variance
  //real<lower=0> var_I; // permanent individual variance
  real<lower=0> var_R; // temporary individual variance (i.e. residual variance)
  //var_G = square(sd_G);
  //var_I = square(sd_I);
  var_R = square(sd_R);
  // Make G-matrix
	matrix[Na,2] G; //  Unscaled blups intercept and res_impact for each individual
  G = LA * Gz * diag_pre_multiply(sd_G, LG)';
  // Make I-matrix 
  matrix[Na,2] i; // permanent individual effects
  i = iz * diag_pre_multiply(sd_I, LI)';
}
model {
  vector[No] z_exp; // expected phenotypic values
  // Priors
  mu ~ normal(mean_z, 2.5 * sd_z); // Weakly informative prior on the intercept
  to_vector(sd_G) ~ exponential(1/sd_z); // Other possible alternatives e.g. student_t(3,0,sd_z), normal(0,sd_z), cauchy(0,sd_z)
  to_vector(sd_I) ~ exponential(1/sd_z); // Ditto, and note that these distributions are truncated to be positive (given the declaration of sd)
  sd_R ~ exponential(1/sd_z); // Ditto
  // Random effect definition (i.e. lower levels of the hierarchy of effects)
  to_vector(Gz) ~ normal(0,1); // implies a ~ normal(0,sd_A)
  to_vector(iz) ~ normal(0,1); // implies i ~ normal(0,sd_I)
  // Expected trait values for each observation
  // Partition the additive genetic breeding values into a_DGE & a_IGE
  for (o in 1:No)
    z_exp[o] = mu + G[animal[o],1] + G[opponent1[o],2] + G[opponent2[o],2] + i[animal[o],1] + i[opponent1[o],2] + i[opponent2[o],2];
  // likelihood function for observed data
  z ~ normal(z_exp,sd_R);
  LG ~ lkj_corr_cholesky(2); //Prior genetic correlations
  LI ~ lkj_corr_cholesky(2); //Prior genetic correlations
}
generated quantities {
  //Variances DGE & IGE
  real var_DGE = sd_G[1]^2;
  real var_IGE = sd_G[2]^2;
  real var_G = var_DGE + var_IGE;
  //Variances
  real var_DIE = sd_I[1]^2;
  real var_IIE = sd_I[2]^2;
  real var_I = var_DIE + var_IIE;
  // Heritabilities
  real<lower=0> var_P = var_G + var_I + var_R;
  real<lower=0> h2_DGE = var_DGE/var_P; // heritability DGE
  real<lower=0> h2_IGE = var_IGE/var_P; // heritability IGE
  real<lower=0> evolvability_DGE = var_DGE/square(mu);// i.e. mean standardized additive genetic variance
  real<lower=0> evolvability_IGE = var_IGE/square(mu);// i.e. mean standardized additive genetic variance
  // Repeatabilities
  real<lower=0> rep_DIE = var_DIE/var_P; // rep DIE
  real<lower=0> rep_IIE = var_IIE/var_P; // rep IIE
  // Correlation genetic level
  matrix[2, 2]  Omega_G;
  Omega_G = LG * LG';
  real cor_G1 = Omega_G[1,2];
  // Correlation PE level
  matrix[2, 2]  Omega_P;
  Omega_P = LI * LI';
  real cor_P1 = Omega_P[1,2];
}", file = "Animal_mod_edited.stan")
#===============================================================================
# Process output

# Relative bias
# Proportion of overlap with 0 (same with 0 distribution for variances)


#===============================================================================


