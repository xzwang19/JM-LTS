
library(rstan)
library(data.table)
library(coda)
library(loo)
library(MASS)
library(dplyr)

#### data generation config ####
random = 1
set.seed(random)

ASSOC = "re" ##### data structure: "re" for shared random effects structure; other options: "cv" for current value structure and "sp" for independent structure
METH = "mcdp_cv" ###### analysis model: "mcdp_cv" for MCDP+CV; other options: "mcdp_re" for MCDP+RE and "mcdp_sp" for separate models

#################################

source(paste0("/data/sim_", ASSOC, ".R")) #### load function to generate data
source("/functions/func_sum.R") ### load function to summarize data


N = 1000 ### sample size
timeinter = 4 ### time interval between visits

datas = get(paste0("sim_", ASSOC))(N, timeinter)

df = datas$df
dat = datas$dl
setDT(df)

print(paste0("Proportion of event is ", mean(dat$ddelta)))

####Parametrization of HMC algorithm-----
nb_iter=6000 #number of iterations
nb_warmup=4000 #number of warm up iterations
nb_chains=3 #number of chains
nb_cores=nb_chains #one chain per core

#### Optional: initial values for MCDP+CV
# init_fun <- function() list(
#   sigma_rshift = 8.5,
#   sigma_b = c(0.45, 0.5),
#   sigma_e = c(0.3, 0.2),
#   shape = 3,
#   
#   alpha = c(-10, -0.5, -0.05, 0, 0),
#   
#   mu1 = c(2.2, 0.3, 0.1),
#   mu2 = c(2.1, 0.15, 0.1),
#   fshift = -2,
#   
#   # random effects start at 0
#   z_rshift = matrix(0, nrow = 1, ncol = dat$N),
#   z_b = matrix(0, nrow = 2, ncol = dat$N),
#   
#   muu1 = -0.06,
#   muu2 = -0.02
# )

fit1=stan(file=paste0("/models/", METH, ".stan"),data=dat,iter=nb_iter,
          warmup=nb_warmup, chains=nb_chains,cores=nb_cores, seed = random) ### optional: init
fit1_sum <- func_sum(fit1, dat, METH, random)




