
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
}
