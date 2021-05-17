data {
  int<lower=0> N;
  int PZ; // Dimension of instruments
  int PX; // Dimension of exogenous
  int PE; // Dimension of endogenous
  
  matrix[N,PZ] Z; // Instruments
  matrix[N,PX] X; // Exogenous
  vector[N] E; // Endogenous
  
  vector[N] y; // Dependent
}
parameters {
  // First stage params
  real alpha_fs;
  vector[PX] beta_fs; //  Exogenous params
  vector[PZ] gamma_fs; // Instrument params
  real<lower=0> sigma_fs;
  
  real alpha_ss;
  vector[PX] beta_ss;
  real gamma_ss;
  real<lower=0> sigma_ss;
}
model {
  vector[N] mu_fs;
  vector[N] mu_ss;
  
  for (n in 1:N){
    mu_fs[n] = alpha_fs + X[n] * beta_fs + Z[n] * gamma_fs;
    mu_ss[n] = alpha_ss + X[n] * beta_ss + mu_fs[n] * gamma_ss;
  }
  
  E ~ normal(mu_fs, sigma_fs);
  y ~ normal(mu_ss, sigma_ss);
  
}
generated quantities{
  vector[N] y_hat;
  vector[N] mu_ss;
  vector[N] mu_fs;
  
  for (n in 1:N){
    mu_fs[n] = alpha_fs + X[n] * beta_fs + Z[n] * gamma_fs;
    mu_ss[n] = alpha_ss + X[n] * beta_ss + mu_fs[n] * gamma_ss;
    y_hat[n] = normal_rng(mu_ss[n], sigma_ss);
  }
}
