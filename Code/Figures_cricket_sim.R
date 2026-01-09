#===============================================================================
# Figures power analysis
#===============================================================================

# Save local package information:
#install.packages("renv")
#renv::init()

# Load local package information:
renv::restore()

library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)

#===============================================================================
# Create sample size tables 

n_rep <- matrix(rep(seq(5,25, by = 5),3),
       nrow = 3, ncol = 5, byrow = TRUE)

n_srn_rep <- matrix(rep(c(1,2,3), each = 5),
       nrow = 3, ncol = 5, byrow = TRUE)

# 2 treatment levels
n_mat_trt2 <- as.data.frame(((12 * 2) * n_rep) * n_srn_rep)
row.names(n_mat_trt2) <- c("1", "2", "3")
colnames(n_mat_trt2) <- c("5", "10", "15", "20", "25")

n_mat_trt2

# 3 treatment levels
n_mat_trt3 <- as.data.frame(((12 * 3) * n_rep) * n_srn_rep)
row.names(n_mat_trt3) <- c("1", "2", "3")
colnames(n_mat_trt3) <- c("5", "10", "15", "20", "25")

n_mat_trt3

#===============================================================================
# Figures model output

### Specify simulated values

# Variance of treatment 
var_trt <- var(rep(1:n_trt, 160))

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

### Load in data
file_names <- paste0("./Output/",(dir("./Output/"))) 

df <- as.data.frame(do.call(rbind,lapply(file_names, function(x){
  d <- read.csv(x)
  d$file <- x 
  d
})))

### Wrangle data

# Separate file name into parameters
df$file <- str_remove_all(df$file, paste(c("./Output/output_", ".csv"), collapse = "|"))

df <- separate(df, file, into = c("Nrep","Treatment","remove","N_rn"), sep = "_",remove = TRUE)

df <- df %>% select(-"remove")

df$Nrep <- str_remove_all(df$Nrep, "nrep")
df$Treatment <- str_remove_all(df$Treatment, "trt")
df$N_rn <- str_remove_all(df$N_rn, "rep")

df$Nrep <- factor(df$Nrep, levels = c("5","10","15","20","25"))

# Remove unnecessary rows 
df <- df %>% filter(!Parameter %in% c("Omega_G[1,1]","Omega_G[2,1]","Omega_G[2,2]"))

# Assess power based on credible intervals overlapping zero
df$Power_YN <- as.numeric(df$X2.5. > 0 | df$X97.5. < 0) 

power_df <- df %>% group_by(Nrep, Treatment, N_rn,Parameter) %>% 
  reframe(Power = sum(Power_YN)/n())

# Get dispersion/uncertainty of credible interval
df$CI_range <- abs(df$X97.5. - df$X2.5.)

### Make figures
power_G_corr <- power_df %>% filter(Parameter == "Omega_G[1,2]")
df_G_corr <- df %>% filter(Parameter == "Omega_G[1,2]")
df_G_corr$X50. <- df_G_corr$X50 -0.5

# Power
png(file="Figures/Power_fig.png", width = 9*1000, height = 6*1000, res = 1000)

ggplot(power_G_corr, aes(y = Power, x = Nrep, color = N_rn)) +
  geom_point(size = 4) +
  geom_line(aes(group = N_rn), linewidth = 1.2) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.025), breaks = seq(0,1,0.2)) +
  xlab("Number of replicates") + 
  geom_hline(yintercept = 0.8, color = "grey34", linetype = "dashed") +
  theme_classic(base_size = 17) + 
  labs(color = "N RN repeats") + 
  scale_color_brewer(palette = "Dark2")

dev.off()

# Bias
png(file="Figures/Bias_fig.png", width = 9*1000, height = 6*1000, res = 1000)

ggplot(df_G_corr, aes(y = X50., x = Nrep, color = N_rn)) +
  geom_boxplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5,0.5), breaks = seq(-0.5, 0.5, 0.25)) +
  ylab("Bias - Genetic correlation") + 
  xlab("Number of replicates") + 
  geom_hline(yintercept = 0, color = "grey34", linetype = "dashed") +
  geom_hline(yintercept = mean(df_G_corr$X50.), color = "darkred", linetype = "dashed") +
  theme_classic(base_size = 17) + 
  labs(color = "N RN repeats") + 
  scale_color_brewer(palette = "Dark2")

dev.off()

# Dispersion
png(file="Figures/Dispersion_fig.png", width = 9*1000, height = 6*1000, res = 1000)

ggplot(df_G_corr, aes(y = CI_range, x = Nrep, color = N_rn)) +
  geom_boxplot() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.5), breaks = seq(0, 1.5, 0.25)) +
  ylab("Uncertainty - Genetic correlation") + 
  xlab("Number of replicates") + 
  geom_hline(yintercept = 0, color = "grey34", linetype = "dashed") +
  theme_classic(base_size = 17) + 
  labs(color = "N RN repeats") + 
  scale_color_brewer(palette = "Dark2")

dev.off()

#===============================================================================
