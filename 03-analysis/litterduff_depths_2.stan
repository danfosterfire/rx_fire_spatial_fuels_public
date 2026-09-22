
data {
  int N; // number of observations
  int P; // number of plots
  int L; // number of distinct locations within each plot
  int Q; // number of covariates: time*compartment + intercept
  
  array[N] int Y; // observed responses
  array[N] int plot_id; // plot IDs
  array[N] int location_id; // location IDs
  
  matrix[P,Q] X; // covariates: time*compartment at the plot level
  
  
  array[L] vector[2] coords; // plot-relative coordinates of each subsample
}

parameters {
  vector[Q] beta; // fixed effect coefficients 
  vector[Q] beta_alpha;
  vector[Q] beta_rho;
  vector[Q] beta_kappa;
  
  matrix[L,P] zGP; // standard-normal deviates for GP effect for each loction,plot,group
  
  real<lower=0> sigmaPlot; // standard deviation of plot effect
  
  vector[P] z_plot; // standard normal deviates for plot effect
  
  
}

transformed parameters {

  vector[P] plot_effect; // atual plot effects
  
  for (i in 1:P){
    plot_effect[i] = sigmaPlot * z_plot[i];
  }
}

model {
  // variable declarations
  vector[N] logmu; // linear predictor
  matrix[L,P] GP; // actual GP effects for each plot:location
  array[P] matrix[L,L] K; // covariance matrix for within-plot GP
  array[P] matrix[L,L] KL; // cholesky decomposition of K
  vector[P] alpha;
  vector[P] rho;
  vector[P] intercept;
  vector[P] kappa;
  
  // priors
  beta ~ normal(0, 2);
  sigmaPlot ~ normal(0, 2);
  beta_kappa ~ normal(0,2);
  beta_alpha ~ normal(0,2);
  beta_rho ~ normal(0, 2);
  
  // random effect realizations
  z_plot ~ std_normal();
  to_vector(zGP) ~ std_normal();
  
  // plot level parameters
  alpha = exp(X * beta_alpha);
  rho = exp(X * beta_rho);
  intercept = X * beta;
  kappa = X * beta_kappa;
  
  // gaussian processses
  for (plot in 1:P){
    K[plot] = gp_exp_quad_cov(coords, alpha[plot], rho[plot]);
    K[plot] = add_diag(K[plot], 1e-9);
    KL[plot] = cholesky_decompose(K[plot]);
    GP[,plot] = KL[plot] * zGP[,plot];
  }
  
  // linear predictors, at observation level
  for (i in 1:N){
    logmu[i] = 
      intercept[plot_id[i]] +
      GP[location_id[i], plot_id[i]] +
      plot_effect[plot_id[i]];
  }
  
  // likelihood, at observation level
  for (i in 1:N){
    Y[i] ~ neg_binomial_2_log(logmu[i], kappa[plot_id[i]]);
  }
}



