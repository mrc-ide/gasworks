template <typename real_type, typename rng_state_type>
__host__ __device__
real_type ll_nbinom(real_type data, real_type model, real_type kappa,
                    real_type exp_noise,
                    rng_state_type& rng_state) {
  if (std::isnan(data)) {
    return 0;
  }
  real_type mu = model +
    dust::random::exponential<real_type>(rng_state, exp_noise);
  return dust::density::negative_binomial_mu(data, kappa, mu, true);
}

template <typename real_type>
__host__ __device__
real_type ll_norm(real_type data, real_type model, real_type sd) {
  if (std::isnan(data)) {
    return 0;
  }
  return dust::density::normal(data, model, sd, true);
}

// [[odin.dust::compare_data(scarlet_fever_cases            = real_type)]]
// [[odin.dust::compare_data(igas_inc                       = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_04    = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_05_14 = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_15_44 = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_45_64 = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_65_74 = real_type)]]
// [[odin.dust::compare_data(daily_scarlet_fever_rate_75    = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate         = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_04      = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_05_14   = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_15_44   = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_45_64   = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_65_74   = real_type)]]
// [[odin.dust::compare_data(daily_pharyngitis_rate_75      = real_type)]]
// [[odin.dust::compare_function]]
template <typename T>
typename T::real_type
compare(const typename T::real_type * state,
        const typename T::data_type& data,
        const typename T::internal_type internal,
        std::shared_ptr<const typename T::shared_type> shared,
        typename T::rng_state_type& rng_state) {
  typedef typename T::real_type real_type;

  const size_t n_group = odin(n_group);

  const real_type ll_sf_cases = ll_nbinom(data.scarlet_fever_cases,
                                          odin(scarlet_fever_cases),
                                          odin(k_hpr), odin(exp_noise),
                                          rng_state);
  const real_type ll_igas = ll_nbinom(data.igas_inc, odin(igas_inc),
                                      odin(k_hpr), odin(exp_noise),
                                      rng_state);

  // Scarlet fever daily rate (always disaggregated)
  const real_type obs_daily_scarlet_fever_rate[] =
    { data.daily_scarlet_fever_rate_04,
      data.daily_scarlet_fever_rate_05_14,
      data.daily_scarlet_fever_rate_15_44,
      data.daily_scarlet_fever_rate_45_64,
      data.daily_scarlet_fever_rate_65_74,
      data.daily_scarlet_fever_rate_75
    };
  const real_type mod_daily_scarlet_fever_rate[] =
    { odin(daily_scarlet_fever_rate_04),
      odin(daily_scarlet_fever_rate_05_14),
      odin(daily_scarlet_fever_rate_15_44),
      odin(daily_scarlet_fever_rate_45_64),
      odin(daily_scarlet_fever_rate_65_74),
      odin(daily_scarlet_fever_rate_75)
    };
  real_type ll_sf_rate = 0.0;
  for (size_t i = 0; i < n_group; ++i) {
    ll_sf_rate += ll_norm(obs_daily_scarlet_fever_rate[i],
                          mod_daily_scarlet_fever_rate[i],
                          odin(k_gp));
  }

  // NOTE: better to arrange the data so we only need to look to see
  // if the aggregated value is better; requires tweaks to your data
  // processing function.
  const bool ll_pharyngitis_aggregated =
    std::isnan(data.daily_pharyngitis_rate_04) &&
    std::isnan(data.daily_pharyngitis_rate_05_14) &&
    std::isnan(data.daily_pharyngitis_rate_15_44) &&
    std::isnan(data.daily_pharyngitis_rate_45_64) &&
    std::isnan(data.daily_pharyngitis_rate_65_74) &&
    std::isnan(data.daily_pharyngitis_rate_75);

  real_type ll_pharyngitis = 0.0;

  if (ll_pharyngitis_aggregated) {
    ll_pharyngitis = ll_norm(data.daily_pharyngitis_rate *
                             odin(etiologic_fraction),
                             odin(daily_gas_pharyngitis_rate),
                             odin(k_gp));
  } else {
    const real_type obs_daily_pharyngitis_rate[] =
      { data.daily_pharyngitis_rate_04,
        data.daily_pharyngitis_rate_05_14,
        data.daily_pharyngitis_rate_15_44,
        data.daily_pharyngitis_rate_45_64,
        data.daily_pharyngitis_rate_65_74,
        data.daily_pharyngitis_rate_75
      };
    const real_type mod_daily_gas_pharyngitis_rate[] =
      { odin(daily_gas_pharyngitis_rate_04),
        odin(daily_gas_pharyngitis_rate_05_14),
        odin(daily_gas_pharyngitis_rate_15_44),
        odin(daily_gas_pharyngitis_rate_45_64),
        odin(daily_gas_pharyngitis_rate_65_74),
        odin(daily_gas_pharyngitis_rate_75)
      };
    const real_type phi_S[] =
      { odin(etiologic_fraction_04),
        odin(etiologic_fraction_05_14),
        odin(etiologic_fraction_15_44),
        odin(etiologic_fraction_45_64),
        odin(etiologic_fraction_65_74),
        odin(etiologic_fraction_75)
      };
    for (size_t i = 0; i < n_group; ++i) {
      ll_pharyngitis += ll_norm(obs_daily_pharyngitis_rate[i] * phi_S[i],
                                mod_daily_gas_pharyngitis_rate[i],
                                odin(k_gp));
    }
  }

  return ll_pharyngitis + ll_sf_rate + ll_sf_cases + ll_igas;
}
