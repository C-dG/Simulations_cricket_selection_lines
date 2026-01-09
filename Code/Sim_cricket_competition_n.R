#===============================================================================
# Simulations competition treatment half-sib full-sib cricket lines 
#===============================================================================

# Save local package information:
#install.packages("renv")
#renv::init()

# Load local package information:
renv::restore()

# Load in packages
library(nadiv)
library(pedtools)

library(mvnfast)
library(Matrix)
library(matrixcalc)

library(rstan)

#===============================================================================
# Set parameters for n simulations

n_sim <- 100 # Number of simulations per parameter set

# Specify number of replicates & number of repeats
repl_levels <- seq(5,25, by = 5)
rn_rep_levels <- c(1,2,3)
trt_levels <- c(2,3)

#===============================================================================
sim<-function(n_rep = 5, # Number of replications 
              n_trt = 3, # Number of levels of competitive treatment
              n_rn_rep = 1 # Number of reaction norm repeats per dam
              ) {
  
# Set parameters for simulation structure

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
var_e <- 1-(var_g_mu + var_g_beta + var_PE_mu + var_PE_beta)

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
  Sire_ID = rep(1:n_dam_rep, each = n_off_dam, times = n_rep), # 2 males per
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

# F2 (Sisters have a unique male)
df_f2 <- data.frame(
  ID = c(unique(df_f3$Dam_ID), unique(df_f3$Sire_ID)),
  Dam_ID = c(rep(df_f1[c(4,7,1,8,2,5),"ID"], each = 2),
             rep(NA,12)), # F3 Sires are unrelated
  Sire_ID = c(rep(df_f1[c(3,6,9),"ID"], each = 4),
              rep(NA,12))) # F3 Sires are unrelated

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

# Check whether matrices are positive definite
matrixcalc::is.positive.definite(A_F3)

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

#===============================================================================
# Simulate values
G_cov <- matrix(c(var_g_mu, r_g,
                  r_g, var_g_beta), 
                nrow = K, ncol = K)  

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

# Simulate PE values 
Sigma_PE <- matrix(c(var_PE_mu, r_PE,
                     r_PE, var_PE_beta),
                   nrow = K, ncol = K)

PE <- as.data.frame(rmvn(n = length(unique(df$ID))
                         , mu = rep(0,K), sigma = Sigma_PE))
colnames(PE) <- c("PE_mu", "PE_beta")
cor(PE) # Check whether the correlation matches set values

# Match with ID in df
PE$ID <- row.names(A_F3)
df <- merge(df, PE, by = "ID")

# Residual variance
df$Res <- rnorm(nrow(df), 0, sqrt(var_e)) 

# Phenotypic equation
df$z <- (mu + df$BV_mu + df$PE_mu) +
  ((beta + df$BV_beta + df$PE_beta)*df$Treatment) +
  df$Res

# Correlations should be ~0.5 
cor(df$BV_mu, df$BV_beta)
cor(df$PE_mu, df$PE_beta)

# Order df ID's  according to A
df <- df[match(row.names(A_F3), df$ID),] 
table(df$ID == row.names(A_F3))

# Give focal's ID that matches with row and column ID of A 
df$ID_index <- as.integer(factor(df$ID, levels = row.names(A_F3)))

# Return Stan data as list (R functions only return one object)
list(No = nrow(df),
     z = as.vector(scale(df$z)),
     sd_z = sd(df$z),
     animal = df$ID_index,
     Na = length(unique(df$ID)),
     Trt = df$Treatment,
     A = A_F3)

}

#===============================================================================
# Simulate data and save for Stan models

# Different combinations of specifications
specs <- as.data.frame(expand.grid(repl_levels, rn_rep_levels, trt_levels))
colnames(specs) <- c("n_rep", "n_rn_rep", "n_trt")
n_specs <- nrow(specs)

# Simulate & save datasets
for (j in 1:n_specs) {
  stan_data <- vector("list", length = n_sim)
  # Simulate data nsim times per n_specs
  for (i in 1:n_sim) {
      stan_data[[i]] <- sim(
        n_rep = specs[j, "n_rep"], # Number of replications 
        n_trt = specs[j, "n_trt"], # Number of levels of competitive treatment
        n_rn_rep = specs[j, "n_rn_rep"] # Number of reaction norm repeats per dam
    )
  }
  saveRDS(stan_data,paste0("Data/stan_data_", 
                           "nrep",specs[j, "n_rep"], "_",
                           "trt", specs[j, "n_trt"], "_",
                           "rn_rep",specs[j, "n_rn_rep"],
                           ".RDS"))
  rm(stan_data)
  gc()
}

#===============================================================================
