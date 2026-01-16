//block defining functions
functions{
  
  /*Longitudinal sub-model: first outcome*/
  real YCOG(real t, real sex, real eduyr, real muu, vector mu, real fshift, real rshift, real b){
    real ycog; // mean outcome variable in longitudinal process
    
    ycog=muu*exp((t + sex*fshift + rshift)/exp(mu[1])) + mu[2] + mu[3]*eduyr + b;
    
    return ycog; 
    
  }

  
  /*Hazard function */
  real ht(real T, real sex, real eduyr, real shape, vector alpha){
    real hr; // hazard function
    
    hr=shape*pow(T, (shape-1))*exp(alpha[1] + alpha[2]*sex + alpha[3]*eduyr);
    return(hr);
  }
  
  
   /* Cumulative hazard function */
  real cht(real T, real sex, real eduyr, real shape, vector alpha) {

    real H; // cumulative hazard

    H = pow(T, shape)*exp(alpha[1] + alpha[2]*sex + alpha[3]*eduyr);
    return H;
  }

}



data{
  
  //longitudinal data
  int nb_obs;
  int<lower=1> N;
  int nb_time[N];
  int<lower=0> nb_max_time;
  vector[nb_obs] measures1;
  vector[nb_obs] measures2;
  vector[nb_obs] TIMEL;
  vector[N] sex;
  vector[N] eduyr;
  //vector[N] covar;

  
  // Numerical computation of the integrate of hazard (Gauss Legendre)
  int order;
  vector[order] points;
  vector[order] weights;

  //survival data
  vector[N] TIMES;
  vector[N] ddelta;

}


parameters{
  //population parameters
  real<lower=-1, upper=0> muu1;
  vector[3] mu1; ### add a covariate
  real<lower=-1, upper=0> muu2;
  vector[3] mu2; ### add a covariate
  real fshift;
  
  real<lower=0> shape;
  vector[3] alpha;
  
  real<lower=0> sigma_rshift;
  vector<lower=0>[2] sigma_e;
  cholesky_factor_corr[2] Lcorr_b;
  vector<lower=0>[2] sigma_b;
	
  matrix[1, N] z_rshift;
  matrix[2, N] z_b;
  

}

transformed parameters{
  
  matrix[N, 2] b; // random time shift
  matrix[2, 2] sigma_b_L;
	matrix[N, 1] rshift; // random intercepts

	for (i in 1:N){
	  rshift[i,1] = sigma_rshift*z_rshift[1,i];
	}                // Random effects Nx2 dimension
  
  sigma_b_L = diag_pre_multiply(sigma_b, Lcorr_b);
	b = (sigma_b_L * z_b)';

  
}


model{

  int w;
  real ycog1k; // mean longitudinal outcome
  real ycog2k;
  
  w=0;

  /*Model log-likelihood*/
  for(i in 1:N){
  
    //each time for each patient in longitudinal process
    for(t in 1:nb_time[i]){
      ycog1k=YCOG(TIMEL[w+t], sex[i], eduyr[i], muu1, mu1, fshift, rshift[i,1], b[i,1]);
      measures1[w+t]~normal(ycog1k, sigma_e[1]);//density of longitudinal process 
      ycog2k=YCOG(TIMEL[w+t], sex[i], eduyr[i], muu2, mu2, fshift, rshift[i,1], b[i,2]);
      measures2[w+t]~normal(ycog2k, sigma_e[2]);//density of longitudinal process 
    }
  
    w=w+nb_time[i];

    if(ddelta[i]==1){//density of survival process in case of death
      target+=log(ht(TIMES[i], sex[i], eduyr[i], shape, alpha)) - cht(TIMES[i], sex[i], eduyr[i], shape, alpha);
    }

    if(ddelta[i]==0){//density of survival process in case of censoring
      target+=-cht(TIMES[i], sex[i], eduyr[i], shape, alpha);

    }

  }
  
  
  // random effect distribution //
  
    to_vector(z_b)~normal(0,1);
    to_vector(z_rshift) ~ normal(0,1);


  // /*priors definition*/
  muu1 ~ uniform(-1,0);
  mu1[1] ~ normal(0,10);
  mu1[2] ~ normal(0,10);
  mu1[3] ~ normal(0,10);
  //mu1[4] ~ normal(0,10);

  muu2 ~ uniform(-1,0);
  mu2[1] ~ normal(0,10);
  mu2[2] ~ normal(0,10);
  mu2[3] ~ normal(0,10);
  //mu2[4] ~ normal(0,10);
  
  fshift ~ normal(0,10);

  //Lcorr_b ~ lkj_corr_cholesky(1);
  sigma_rshift ~ cauchy(0,5);
  Lcorr_b ~ lkj_corr_cholesky(1);
  sigma_b ~ cauchy(0,0.5);
  sigma_e ~ cauchy(0,0.5);
  
  shape ~ lognormal(0,1.5);
  alpha[1] ~ normal(-10,10);
  alpha[2] ~ normal(0,10);
  alpha[3] ~ normal(0,10);

}

generated quantities{
 matrix[2, 2] CORR_b;

	int w;
  real ycogn1;
  real ycogn2;
  vector[N] log_lik;
  real normloglik1;
  real normloglik2;
  real sumnormloglik;


// generate estimated covariance matrix

	CORR_b = multiply_lower_tri_self_transpose(Lcorr_b);


// generate piecewise loglikelihood
  w=0;

    for(i in 1:N){
    //random effects of patient i
    sumnormloglik = 0;

    //each time for each patient
    for(t in 1:nb_time[i]){
      ycogn1=YCOG(TIMEL[w+t], sex[i], eduyr[i], muu1, mu1, fshift, rshift[i,1], b[i,1]);
      ycogn2=YCOG(TIMEL[w+t], sex[i], eduyr[i], muu2, mu2, fshift, rshift[i,1], b[i,2]);
      normloglik1=normal_lpdf(measures1[w+t] | ycogn1, sigma_e[1]);//density of longitudinal process
      normloglik2=normal_lpdf(measures2[w+t] | ycogn2, sigma_e[2]);//density of longitudinal process
      sumnormloglik = sumnormloglik + normloglik1 + normloglik2;
    }

		if(ddelta[i]==1){//density of survival process in case of death
      log_lik[i] = sumnormloglik + log(ht(TIMES[i], sex[i], eduyr[i], shape, alpha)) - cht(TIMES[i], sex[i], eduyr[i], shape, alpha);
    }

    if(ddelta[i]==0){//density of survival process in case of censoring
      log_lik[i] = sumnormloglik -cht(TIMES[i], sex[i], eduyr[i], shape, alpha);
    }

    w=w+nb_time[i];

}
	
}








