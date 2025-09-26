data {
  int<lower=1> N; // number of samples
  int<lower = 1> K_total; // total number of components
  array[N, 3] int<lower=1> index_components; // number of components
  vector<lower = 0, upper = 1>[K_total] w; // mass fractions of components (sum to 1)
  vector<lower = 0>[N] gamma_mirs_obs_p1;
  vector<lower = 0>[N] gamma_mirs_obs_p2;
  real b_intercept_p1;
  real<lower = 0.0> b_intercept_p2;
  
  // parameters
  vector<lower = 0>[K_total] D_p1; 
  vector<lower = 0>[K_total] D_p2;
  vector<lower = 0>[N] phi_p1;
  vector<lower = 0>[N] phi_p2;
  vector<lower = 0>[N] phi_p3;
}

parameters {
  vector<lower=0, upper=1>[K_total] D; // degrees of decomposition for each component, bounded between 0 and 1
  vector<lower=0>[N] phi; // observation noise
  vector<lower = 0.0, upper = 1.0>[N] gamma_mirs_obs;
  real b_intercept;
}

transformed parameters {
  vector<lower=0, upper = 1>[N] gamma_mirs;
  {
    vector[K_total] gamma_mirs_component;
    gamma_mirs_component = logit(D) - b_intercept;
    for(n in 1:N) {
      int index_from = index_components[n, 2];
      int index_to = index_components[n, 3];
      gamma_mirs[n] = inv_logit(dot_product(w[index_from:index_to], gamma_mirs_component[index_from:index_to]) + b_intercept);
    }
  }
}

model {
  // Priors
  D ~ beta(D_p1, D_p2); // example prior, can be changed based on domain knowledge
  phi ~ gamma(phi_p1, phi_p2 .* phi_p3); // small noise
  b_intercept ~ normal(b_intercept_p1, b_intercept_p2);

  // Likelihood
  // The overall A is a weighted sum of component A's
  {
    gamma_mirs_obs_p1 ~ beta(gamma_mirs_obs .* gamma_mirs_obs_p2, (1.0 - gamma_mirs_obs) .* gamma_mirs_obs_p2);
    vector[N] phi_scaled = phi .* phi_p3;
    gamma_mirs_obs ~ beta(gamma_mirs .* phi_scaled, (1.0 - gamma_mirs) .* phi_scaled);
  }

}
