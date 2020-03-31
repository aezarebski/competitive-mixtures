functions {
  /*
  Return the model prediction for the recipient given a donor proportion and the
  model parameter theta.
  */
  real f(real theta, real d) {
    return d / (d + (1 - d) * exp(theta));
  }

  /*
  Return a clipped x into [a, b]
  */
  real clipper(real x, real a, real b) {
    if (x > b) {
      return b;
    } else if (x < a) {
      return a;
    } else {
      return x;
    }
  }
}

data {
  int<lower=1> num_obs;
  real<lower=0,upper=1> d_hat[num_obs];
  real<lower=0,upper=1> r_hat[num_obs];
  real<lower=0> sigma_pyro;
  real<lower=0> sigma_tau;
}

transformed data {

}

parameters {
  real theta;
  real<lower=0,upper=1> d[num_obs];
  real tau[num_obs];
}

transformed parameters {
  real<lower=0,upper=1> r[num_obs];
  for (ix in 1:num_obs)
    r[ix] = clipper(f(theta, d[ix]) + tau[ix], 0, 1);
}

model {
  // Prior distribution
  target += normal_lpdf(theta | 0, 1);
  for (ix in 1:num_obs) {
    target += uniform_lpdf(d[ix] | 0, 1);
    target += normal_lpdf(tau[ix] | 0, sigma_tau);
    if (d_hat[ix] == 0) {
      target += normal_lcdf(0 | d[ix], sigma_pyro);
    } else if (d_hat[ix] == 1) {
      target += normal_lccdf(1 | d[ix], sigma_pyro);
    } else {
      target += normal_lpdf(d_hat[ix] | d[ix], sigma_pyro);
    }
    if (r_hat[ix] == 0) {
      target += normal_lcdf(0 | r[ix], sigma_pyro);
    } else if (r_hat[ix] == 1) {
      target += normal_lccdf(1 | r[ix], sigma_pyro);
    } else {
      target += normal_lpdf(r_hat[ix] | r[ix], sigma_pyro);
    }
  }
}

generated quantities {

}
