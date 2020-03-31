functions {
  real[] petrie2015(
    real t,
    real[] y,
    real[] theta,
    real[] x_r,
    int[] x_i) {
    /*
    Return the differential for the model given in Petrie *et al* (2015).

    Note that this differential is defined in terms of the unscaled parameters
    variables

    y = (T,          // 1
        L_A,        // 2
        L_B,        // 3
        I_A,        // 4
        I_B,        // 5
        VTCID_A,    // 6
        VTCID_B,    // 7
        VRNA_A,     // 8
        VRNA_B)     // 9
    theta = (beta_A,    // 1
            beta_B,    // 2
            k_A,       // 3
            k_B,       // 4
            delta_A,   // 5
            delta_B,   // 6
            p_A,       // 7
            p_B,       // 8
            ch_A,      // 9
            ch_B,      // 10
            dinf_A,    // 11
            dinf_b,    // 12
            xi_A,      // 13
            xi_B)      // 14
    */
    real dydt[9];
    // dTdt = - (beta_A * VTCID_A - beta_B * VTCID_B) T
    dydt[1] = - (theta[1] * y[6] + theta[2] * y[7]) * y[1];
    // dL_Adt = beta_A * T * VTCID_A - k_A * L_A
    dydt[2] = theta[1] * y[1] * y[6] - theta[3] * y[2];
    // dL_Bdt = beta_B * T * VTCID_B - k_B * L_B
    dydt[3] = theta[2] * y[1] * y[7] - theta[4] * y[3];
    // dI_Adt = k_A * L_A - delta_A * I_A
    dydt[4] = theta[3] * y[2] - theta[5] * y[4];
    // dI_Bdt = k_B * L_B - delta_B * I_B
    dydt[5] = theta[4] * y[3] - theta[6] * y[5];
    // dVTCID_Adt = p_A * I_A - (ch_A + dinf_A) * VTCID_A
    dydt[6] = theta[7] * y[4] - (theta[9] + theta[11]) * y[6];
    // dVTCID_Bdt = p_B * I_B - (ch_B + dinf_B) * VTCID_B
    dydt[7] = theta[8] * y[5] - (theta[10] + theta[12]) * y[7];
    // dVRNA_Adt = xi_A * p_A * I_A - ch_A * VRNA_A
    dydt[8] = theta[13] * theta[7] * y[4] - theta[9] * y[8];
    // dVRNA_Bdt = xi_B * p_B * I_B - ch_B * VRNA_B
    dydt[9] = theta[14] * theta[8] * y[5] - theta[10] * y[9];
    return dydt;
  }

  real[,,] model_solution(
    int num_ferrets,
    int[] is_pure,
    int max_num_obs,
    real t_0,
    real[,] tiv_0,
    real[] theta,
    real[] obs_times,
    int obs_dim,
    real[] x_r,
    int[] x_i) {
    /*
    The TIV model from Petrie *et al* (2015)

    Note that this function returns the ODE in the natural coordinates, this is
    the actual values for the TCID and RNA values and the proportions for the
    PYRO.
    */
    real model_sol[num_ferrets, obs_dim, max_num_obs];
    real return_array[num_ferrets, obs_dim, max_num_obs];
    real tiv_sol[max_num_obs,9];
    for (ferret in 1:num_ferrets) {
      tiv_sol = integrate_ode_bdf(
        petrie2015,
        tiv_0[ferret,:],
        t_0,
        obs_times,
        theta,
        x_r,
        x_i);
      for (obs in 1:max_num_obs) {
        /*
        These values are due to the particular specification of the state vector
        which is linked to the initial condition and the differential.
        */
        model_sol[ferret, 1, obs] = tiv_sol[obs, 6] + tiv_sol[obs, 7];
        model_sol[ferret, 2, obs] = tiv_sol[obs, 8] + tiv_sol[obs, 9];
        model_sol[ferret, 3, obs] = tiv_sol[obs, 9] / model_sol[ferret, 2, obs];
      }
    }
    /*
    Numerical issues can lead to small, but negative, solutions. However, for
    all positive times the solution should be positive. Since this only effects
    pathological parameter values we assign machine precision to these
    non-positive values.
     */
    for (ix in 1:num_ferrets) {
      for (jx in 1:obs_dim) {
	for (kx in 1:max_num_obs) {
	  if (model_sol[ix, jx, kx] > 0) {
	    return_array[ix, jx, kx] = model_sol[ix, jx, kx];
	  } else {
	    return_array[ix, jx, kx] = machine_precision();
	  }
	}
      }
    }
    return return_array;
  }


  real ln_mu(real a, real b) {
    return 0.5 * (a + b) / log10(e());
  }

  real ln_sigma(real a, real b) {
    return 0.25 * (b - a) / log10(e());
  }
}

data {
  real init_target_cells;
  int num_ferrets;
  int num_obs[num_ferrets];
  /*
  The ferrets can be infected with one of two strains or a mixture of the two.
  The data `is_pure` describes the mixture of virus used to infect a given
  ferret.
  -1 indicates pure wild
  0 indicates a mixture
  +1 indicates pure mutant
  */
  int<lower=-1,upper=1> is_pure[num_ferrets];
  int num_pure_wild;
  int num_pure_mutant;
  int num_mix;
  int max_rna_complete;
  real rna_complete[num_ferrets, max_rna_complete];
  int rna_complete_ixs[num_ferrets, max_rna_complete];
  int rna_ixs[num_ferrets];
  int max_tcid_complete;
  real tcid_complete[num_ferrets, max_tcid_complete];
  int tcid_complete_ixs[num_ferrets, max_tcid_complete];
  int max_tcid_censored;
  real tcid_censored[num_ferrets, max_tcid_censored];
  int tcid_censored_ixs[num_ferrets, max_tcid_censored];
  int tcid_ixs[2, num_ferrets];
  int max_pyro_complete;
  real pyro_complete[num_ferrets, max_pyro_complete];
  int pyro_complete_ixs[num_ferrets, max_pyro_complete];
  int max_pyro_left_censored;
  real pyro_left_censored[num_ferrets, max_pyro_left_censored];
  int pyro_left_censored_ixs[num_ferrets, max_pyro_left_censored];
  int max_pyro_right_censored;
  real pyro_right_censored[num_ferrets, max_pyro_right_censored];
  int pyro_right_censored_ixs[num_ferrets, max_pyro_right_censored];
  int pyro_ixs[3, num_ferrets];
  /*
  The standard deviation of the measurement error for each of the measurement
  types is defined by:
  */
  real sigma_tcid;
  real sigma_rna;
  real sigma_pyro;  // Calibration suggests this should be 0.05
}
transformed data {
  /*
  Bare in mind that the transformed data block is only run once after the data
  is read in and not at every sample. Therefore we can afford to write less
  efficient code here because it is only evaluated once.
  */
  real t_0;
  int obs_dim;
  real obs_times[max(num_obs)];
  real obs_thresh_tcid;
  real obs_thresh_rna;
  real zero_infection_states[8];
  real x_r[0];
  int x_i[0];
  real init_target_cells_array[num_ferrets];
  real zero_ic_array[num_ferrets];
  int pure_wild_jx;
  int pure_wild_ixs[num_pure_wild];
  int pure_mutant_jx;
  int pure_mutant_ixs[num_pure_mutant];
  int num_pure_ferrets;
  int mix_jx;
  int mix_ixs[num_mix];
  t_0 = 0;  // Starting time for the ODE
  obs_dim = 3; // TCID, RNA and pyrosequencing
  for (obs in 1:max(num_obs)) {
    obs_times[obs] = obs;
  }
  /*
  Minimum values for a proper measurement.
  The `obs_thresh_rna` was a guess based on the minimium value observed in the
  data.
  */
  obs_thresh_tcid = 0.5;
  obs_thresh_rna = 0.5;
  /*
  Define an array for use in the transformed solver.
  */
  for (ix in 1:8)
    zero_infection_states[ix] = 0;

  num_pure_ferrets = num_pure_wild + num_pure_mutant;
  /*
  This loop constructs the index arrays for subsetting each type of ferret
  infection. The arrays `pure_wild_ixs`, `pure_mutant_ixs` and `mix_ixs` are
  used to index just these ferrets in the data.
  */
  pure_wild_jx = 1;
  pure_mutant_jx = 1;
  mix_jx = 1;
  for (ferret in 1:num_ferrets) {
    /*
    These arrays are used in the transformed parameters block when constructing
    the initial condition.
    */
    init_target_cells_array[ferret] = init_target_cells;
    zero_ic_array[ferret] = 0;
    /*
    Construct the arrays for indexing ferrets with a particular mixture
    */
    if (is_pure[ferret] == -1) {
      pure_wild_ixs[pure_wild_jx] = ferret;
      pure_wild_jx = pure_wild_jx + 1;
    } else if (is_pure[ferret] == 0) {
      mix_ixs[mix_jx] = ferret;
      mix_jx = mix_jx + 1;
    } else {
      pure_mutant_ixs[pure_mutant_jx] = ferret;
      pure_mutant_jx = pure_mutant_jx + 1;
    }
  }
}

parameters {
  /*
  Each ferret has its own initial condition

  To improve the efficiency of the sampler we seperate the initial condition
  into two arrays: one for the ferrets infected with a single strain and one for
  the ferrets infected with a mixture.

  The parameters of the system are stored in `theta`. The same parameters are
  shared by all of the ferrets.
  */
  real<lower=0> VRNA_0_pure[num_pure_ferrets];
  real<lower=0> VTCID_0_pure[num_pure_ferrets];
  real<lower=0> VRNA_0_mix[num_mix, 2];
  real<lower=0> VTCID_0_mix[num_mix, 2];
  real<lower=0> theta[14];  // Parameters shared by ferrets
}

transformed parameters {
  /*
  Construct the initial condition for each of the ferrets depending on whether
  they are a mixture ferret or not.

  Solve the system for each of the ferrets.
  */
  real<lower=0> initial_condition[num_ferrets, 9];
  real<lower=0> model_sol[num_ferrets, obs_dim, max(num_obs)];
  // Initial condition shared by all ferrets
  initial_condition[:, 1] = init_target_cells_array;
  initial_condition[:, 2] = zero_ic_array;
  initial_condition[:, 3] = zero_ic_array;
  initial_condition[:, 4] = zero_ic_array;
  initial_condition[:, 5] = zero_ic_array;
  /*
  Construct the components of the initial condition specific to each of the
  infections
  */
  initial_condition[pure_wild_ixs, 6] = VTCID_0_pure[1:num_pure_wild];
  initial_condition[pure_wild_ixs, 8] = VRNA_0_pure[1:num_pure_wild];
  initial_condition[pure_wild_ixs, 7] = zero_ic_array[1:num_pure_wild];
  initial_condition[pure_wild_ixs, 9] = zero_ic_array[1:num_pure_wild];

  initial_condition[mix_ixs, 6:7] = VTCID_0_mix;
  initial_condition[mix_ixs, 8:9] = VRNA_0_mix;

  initial_condition[pure_mutant_ixs, 7] = VTCID_0_pure[(num_pure_wild+1):num_pure_ferrets];
  initial_condition[pure_mutant_ixs, 9] = VRNA_0_pure[(num_pure_wild+1):num_pure_ferrets];
  initial_condition[pure_mutant_ixs, 6] = zero_ic_array[1:num_pure_mutant];
  initial_condition[pure_mutant_ixs, 8] = zero_ic_array[1:num_pure_mutant];
  // Solve the actual system.
  model_sol = model_solution(
    num_ferrets,
    is_pure,
    max(num_obs),
    t_0,
    initial_condition,
    theta,
    obs_times,
    obs_dim,
    zero_infection_states,
    x_i);
}

model {
  /*
  Prior distribution

  This follows from the idea that if
  log10(X) ~ uniform(a, b)
  then
  X ~lognormal((a+b)/(2*log10(e)), (b-a)/(4*log10(e)))
  */
  VTCID_0_pure ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  VRNA_0_pure ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  VRNA_0_mix[:,1] ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  VRNA_0_mix[:,2] ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  VTCID_0_mix[:,1] ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  VTCID_0_mix[:,2] ~ lognormal(ln_mu(-4,0), ln_sigma(-4,0));
  theta[1:2] ~ lognormal(ln_mu(-9,-1), ln_sigma(-9,-1));
  theta[3:4] ~ lognormal(ln_mu(0,1.4), ln_sigma(0,1.4));
  theta[5:6] ~ lognormal(ln_mu(0.38,0.9), ln_sigma(0.38,0.9));
  theta[7:8] ~ lognormal(ln_mu(-4,2), ln_sigma(-4,2));
  theta[9:10] ~ lognormal(ln_mu(-2,3), ln_sigma(-2,3));
  theta[11:12] ~ lognormal(ln_mu(0.48,0.50), ln_sigma(0.48,0.50));
  theta[13:14] ~ lognormal(ln_mu(0,6), ln_sigma(0,6));
  /*
  Ferret specific parts of the target
  */
  for (ferret in 1:num_ferrets) {
    /*
    TCID sampling statements.

    Note that Section 5.2 of the stan docs state the target += ... statment will
    add the sum of the components when a container argument is supplied.
    */
    tcid_complete[ferret, 1:tcid_ixs[1, ferret]]
      ~ normal(
          log10(
            model_sol[
              ferret,
              1,
              tcid_complete_ixs[ferret, 1:tcid_ixs[1, ferret]]
            ]),
          sigma_tcid
        );
    target += normal_lcdf(
      tcid_censored[ferret, 1:tcid_ixs[2, ferret]] |
      log10(
        model_sol[
          ferret,
          1,
          tcid_censored_ixs[ferret, 1:tcid_ixs[2, ferret]]
        ]),
      sigma_tcid
    );
    /*
    RealTime RNA sampling statements.

    There is no censor statement because measurements below the measurement
    threshold are listed as -4 rather than the minimum read.
    */
    rna_complete[ferret, 1:rna_ixs[ferret]]
      ~ normal(
        log10(model_sol[
          ferret,
          2,
          rna_complete_ixs[ferret, 1:rna_ixs[ferret]]
        ]),
        sigma_rna
      );
    if (is_pure[ferret] == 0) {
      /*
      Pyrosequencing sampling statements

      Note that these sampling statements should only be called when we are 
      considering a ferret that has a mixture and there are actuall measurements
      (hence the if statments dependent upon `pyro_ixs`.
      */
      if (pyro_ixs[1, ferret] >= 1) {
        pyro_complete[ferret, 1:pyro_ixs[1, ferret]]
          ~ normal(
              model_sol[
                ferret,
                3,
                pyro_complete_ixs[ferret, 1:pyro_ixs[1, ferret]]
              ],
              sigma_pyro
            );
      }
      if (pyro_ixs[2, ferret] >= 1) {
        target += normal_lcdf(
          pyro_left_censored[ferret, 1:pyro_ixs[2, ferret]] |
            model_sol[
              ferret,
              3,
              pyro_left_censored_ixs[ferret, 1:pyro_ixs[2, ferret]]
            ],
          sigma_pyro
          );
      }
      if (pyro_ixs[3, ferret] >= 1) {
        target += normal_lccdf(
          pyro_right_censored[ferret, 1:pyro_ixs[3, ferret]] |
            model_sol[
              ferret,
              3,
              pyro_right_censored_ixs[ferret, 1:pyro_ixs[3, ferret]]
            ],
          sigma_pyro
          );
      }
    }
  }
}

generated quantities {
  real sol[num_ferrets, obs_dim, max(num_obs)] = model_sol;
}
