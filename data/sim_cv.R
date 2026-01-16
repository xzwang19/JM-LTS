###### parameter setting #####

sim_cv <- function(N, timeinter){
  
  ##### parameter for longitudinal
  #### fixed effects
  library(truncnorm)
  library(MASS)
  library(dplyr)
  
  N= N
  
  l1 <- -0.06
  g1 <- 2.15
  v1 <- 0.3
  z1 <- 0.1
  
  l2 <- -0.02
  g2 <- 2.10
  v2 <- 0.15
  z2 <- 0.1
  
  fshift <- -2
  
  #### random effects for intercepts
  b.cov <- matrix(c(0.45^2, 0.45*0.5*0.6, 0.45*0.5*0.6, 0.5^2), ncol = 2)
  
  #### random effects for latent time shift
  sigma.rshift <- 8.5
  
  #### random error
  sigma.1 <- 0.3
  sigma.2 <- 0.2
  
  ##### parameters for survival part
  shape1 = 2.8
  alpha1 = -10
  alpha2 = -0.2
  alpha3 = 0.03
  alpha4 = -1.1
  alpha5 = 0.15

  
  #### simulation of longitudinal part #####
  
  f.1 = function(t, sex, eduyr, l1, g1, v1, z1, fshift, rshift, b1) {
    
    y.1 = l1*exp((t+sex*fshift+rshift)/exp(g1)) + z1*eduyr + v1 + b1
    return(y.1)
  }
  
  f.2 = function(t, sex, eduyr, l2, g2, v2, z2, fshift, rshift, b2) {
    
    y.2 = l2*exp((t+sex*fshift+rshift)/exp(g2)) + z2*eduyr + v2 + b2
    return(y.2)
  }
  
  ###### generate y.1 and y.2 ########
  
  longi.y <- data.frame(ID = NA, value = NA, type = NA, t = NA, sex = NA, eduyr = NA, fshiftx = NA, rshift = NA, b0 = NA) ### longitudinal format
  
  for (i in 1:N){
    
    sex = rbinom(1, 1, 0.5)
    eduyr = round(runif(1, -5, 5),2)
    t0 = round(runif(1, min = -5, max = 5),2)
    bs <- mvrnorm(n = 1, rep(0, 2), b.cov)
    b.10 <- bs[1]
    b.20 <- bs[2]
    rshift = rnorm(1, mean = 0, sd = sigma.rshift)
    fshiftx = fshift*sex
    
    for (t in seq(t0, t0 + 20, by = timeinter)) {
      mean.y.1 = f.1(t, sex, eduyr, l1, g1, v1, z1, fshift, rshift, b.10)
      epsilon.1 = rnorm(1, 0, sigma.1)
      y.1 = mean.y.1 + epsilon.1
      
      mean.y.2 = f.2(t, sex, eduyr, l2, g2, v2, z2, fshift, rshift, b.20)
      epsilon.2 = rnorm(1, 0, sigma.2)
      y.2 = mean.y.2 + epsilon.2
      
      longi.y <- rbind(longi.y, c(i, y.1, 1, t, sex, eduyr, fshiftx, rshift, b.10), c(i, y.2, 2, t, sex, eduyr, fshiftx, rshift, b.20))   
    }
  }
  
  
  longi.y <- na.omit(longi.y)
  setDT(longi.y)
  
  ######## simulation of survival part ############
  
  ######### generate survival data #############
  
  n <- N
  
  ########### data initialization #########
  
  pass01<-rep(0,n) 
  
  obs<-c(1:n);time0<-rep(0,n);
  
  ## database: initialization at 0
  database<-data.frame(obs,TIME0=time0,ddelta=pass01)
  
  j <- 0
  for (k in 1:n){
    j <- j+1
    
    years1 <- max(longi.y[ID==j]$t)
    b.10 <- longi.y[type==1 & ID==j]$b0[1]
    b.20 <- longi.y[type==2 & ID==j]$b0[1]
    rshift <- longi.y[type==1 & ID==j]$rshift[1]
    sex <- longi.y[type==1 & ID==j]$sex[1]
    eduyr <- longi.y[type==1 & ID==j]$eduyr[1]
    
    m1 <- function(t){
      
      y.1 <- l1*exp((t+sex*fshift+rshift)/exp(g1)) + v1 + z1*eduyr + b.10
      
      return(y.1)
    }
    
    m2 <- function(t){
      
      y.1 <- l2*exp((t+sex*fshift+rshift)/exp(g2)) + v2 + z2*eduyr + b.20
      
      return(y.1)
    }
    
    haz01 <- function(t){
      
      hazard <- shape1*t**(shape1-1)*exp(alpha1 + alpha2*sex + alpha3*eduyr + alpha4*m1(t) + alpha5*m2(t))
      return(hazard)
    }
    
    invS01 <- function (t, u) {
      integrate(haz01, lower = 0, upper = t)$value + log(u)
    }
    
    ### t01 generation    
    u <- runif(1)
    
    Up <- 30
    tries <- 5
    Root01 <- try(uniroot(invS01, interval = c(1e-05, Up), u = u)$root, TRUE) 
    while(inherits(Root01, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 30
      Root01 <- try(uniroot(invS01, interval = c(1e-05, Up), u = u)$root, TRUE)
    }
    t1 <- if (!inherits(Root01, "try-error")) Root01 else NA
    
    if(t1<=years1){
      bh<-c(j, t1, 1)
    }
    
    if(years1 < t1){
      bh<-c(j, years1, 0)
    }
    database[j,]<-bh
  }
  
  setDT(database)
  database$ID <- database$obs
  SURVTTE <- database
  sum(SURVTTE$ddelta)
  
  ####### combine longitudinal with survival ##########
  
  longsurv <- inner_join(longi.y, SURVTTE, by="ID")
  setDT(longsurv)
  longsurv <- longsurv[t<=TIME0]
  
  ######### create a dat list #######
  
  data = longsurv
  data1 = longsurv[type==1]
  data2 = longsurv[type==2]
  
  nb_obs = length(data1$ID)#total number of observations
  N = length(unique(data1$ID))#number of patients
  
  nb_time = c()#number of longitudinal measurements per patient
  
  for (i in 1:N) {
    id <- unique(data1$ID)[i]
    nb_time[i] = length(data1[which(data1$ID == id), ]$t)
    
  }
  
  nb_max_time = max(nb_time) #maximum number
  
  measures1 = data1$value 
  measures2 = data2$value
  
  TIMEL = data1$t # timescale is age since 60
  setDT(data1)
  
  TIMES = data1[!duplicated(ID)]$TIME0 
  
  ddelta = data1[!duplicated(ID)]$ddelta
  
  sex = data1[!duplicated(ID)]$sex
  eduyr = data1[!duplicated(ID)]$eduyr
  
  ###integration parameters
  
  x<-c(-0.960289856,
       -0.796666477,
       -0.52553241,
       -0.183434642,
       0.183434642,
       0.52553241,
       0.796666477,
       0.960289856
  )
  w<-c(0.101228536,
       0.222381034,
       0.313706646,
       0.362683783,
       0.362683783,
       0.313706646,
       0.222381034,
       0.101228536
  )
  
  dat = list(
    nb_obs = nb_obs,
    N = N,
    nb_time = nb_time,
    nb_max_time = nb_max_time,
    measures1 = measures1,
    measures2 = measures2,
    sex = sex,
    eduyr = eduyr,
    TIMEL = TIMEL,
    TIMES = TIMES,
    ddelta = ddelta,
    #  covar = covar,
    order = 8 ,
    points = x,
    weights = w
  )
  
  return(list("df" = longsurv, "dl" = dat))
}



