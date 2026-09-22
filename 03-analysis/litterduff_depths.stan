
data {
  int N; // number of observations
  int P; // number of plots
  int L; // number of distinct locations within each plot
  
  array[N] int Y; // observed responses
  array[N] int plot_id; // plot IDs
  array[N] int location_id; // location IDs
  
  matrix[N,6] X; // fixed effects design matrix, giving dummy variables for 
                 // COMP (A, B, or C) crossed with TIMESTEP (pre and post), with 
                 // A:PRE the intercept
  
  array[N] int pre; // indicator 1 if pre 0 if post
  array[N] int post; // timestep explanatory variable; 1 if post, 0 otherwise
  
  array[L] vector[2] coords; // plot-relative coordinates of each subsample
}

parameters {
  vector[6] beta; // fixed effect coefficients 
  matrix[L,P] zGP_pre; // standard-normal deviates for GP effect prefire
  matrix[L,P] zGP_post; // standard normal deviates for GP effect postfire
  
  real<lower=0> alpha_pre; // standard deviation of prefire GP
  real<lower=0> rho_pre; // length-scales of prefire GP
  
  real<lower=0> alpha_post; // standard deviation of postfire GP
  real<lower=0> rho_post; // length-scale of postfire GP
  
  real<lower=0> sigmaPlot; // standard deviation of plot effect
  vector[P] z_plot; // standard normal deviates for plot effect
  
  real<lower=0> kappa; // NB dispersion parameter
  
  
}

transformed parameters {

  vector[P] plot_effect; // atual plot effects

  
  plot_effect = sigmaPlot * z_plot;
  
}

model {
  // variable declarations
  matrix[L,P] GP_pre; // actual GP effects for each plot:location
  matrix[L,P] GP_post; // actual GP effects for postfire
  matrix[L,L] K_pre; // covariance matrix for within-plot GP (prefire)
  matrix[L,L] KL_pre; // cholesky decomposition of K_pre
  matrix[L,L] K_post; // covariance matrix for within-plot GP (postfire)
  matrix[L,L] KL_post; // cholesky decomposition of K_post
  
  vector[N] logmu; // linear predictor
  
  // priors
  beta ~ normal(0, 5);
  sigmaPlot ~ normal(0, 5);
  kappa ~ cauchy(0,5);
  alpha_pre ~ normal(0, 5);
  alpha_post ~ normal(0, 5);
  rho_pre ~ inv_gamma(5,40);
  rho_post ~ inv_gamma(5,40);
  
  // random effect realizations
  z_plot ~ std_normal();
  to_vector(zGP_pre) ~ std_normal();
  to_vector(zGP_post) ~ std_normal();
  
  
  // gaussian processse
  K_pre = gp_exp_quad_cov(coords, alpha_pre, rho_pre);
  K_pre = add_diag(K_pre, 1e-9);
  KL_pre = cholesky_decompose(K_pre);
  GP_pre = KL_pre * zGP_pre;
  
  K_post = gp_exp_quad_cov(coords, alpha_post, rho_post);
  K_post = add_diag(K_post, 1e-9);
  KL_post = cholesky_decompose(K_post);
  GP_post = KL_post * zGP_post;
  
  // linear predictors
  for (i in 1:N){
    logmu[i] = 
      X[i,] * beta +
      GP_pre[location_id[i], plot_id[i]] * pre[i] +
      GP_post[location_id[i], plot_id[i]] * post[i] +
      plot_effect[plot_id[i]];
    
  }
  
  
  // likelihood
  Y ~ neg_binomial_2_log(logmu, kappa);
}

//generated quantities {
//  vector[N] log_lik;
  
//  for (i in 1:N){
//    log_lik[i] = neg_binomial_2_lpmf(Y[i]|exp(logmu[i]),kappa);
//  }
//}
