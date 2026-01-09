
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
// beta
real beta_back = beta*sd_z;

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
}
