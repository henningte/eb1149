#### Metadata on models ####

#' Defines which models to compute and which settings to use
#' 
#' @export
dd_get_model_info <- function(dd_data_model_preprocessed, dd_mir_preprocessing_settings) {
  
  d <- dd_data_model_preprocessed
  d_spectral_preprocessing_settings <- dd_mir_preprocessing_settings
  
  tidyr::expand_grid(
    id_preprocessing = d_spectral_preprocessing_settings$id_preprocessing, # compute one model for each preprocessing variant
  ) |>
    dplyr::mutate(
      id_model = seq_along(id_preprocessing),
      id_sample_all = 
        purrr::map(id_model, function(i) {
          d[[id_preprocessing[[i]]]] |>
            dplyr::pull(id_sample)
        }),
      sample_size_all = purrr::map_int(id_sample_all, length),
      do_compute_model = TRUE,
      validation_mode = "loo",
      y = 
        purrr::map(seq_along(id_model), function(i) {
          .id_sample_all <- id_sample_all[[i]]
          .id_preprocessing <- id_preprocessing[[i]]
          d[[.id_preprocessing]] |>
            dplyr::filter(id_sample %in% .id_sample_all) |>
            dplyr::select(dplyr::all_of(c("mass_relative_mass"))) |>
            dplyr::mutate(
              mass_relative_mass =
                dplyr::case_when(
                  mass_relative_mass == 1 ~ mass_relative_mass - 0.001, #---note beta distribution assumes values are in (0,1)
                  TRUE ~ mass_relative_mass
                ),
              mass_relative_mass = 1 - mass_relative_mass
            )
        }),
      x = 
        purrr::map(seq_along(id_model), function(i) {
          .id_sample_all <- id_sample_all[[i]]
          .id_preprocessing <- id_preprocessing[[i]]
          d[[.id_preprocessing]] |>
            dplyr::filter(id_sample %in% .id_sample_all) |>
            dplyr::select(id_sample, spectra)
        }),
      y_center = 0,
      y_scale = 1,
      x_train =
        purrr::map(x, function(.x) {
          .x |>
            ir::ir_flatten() |>
            dplyr::select(-1) |>
            t() |>
            scale(center = TRUE, scale = TRUE)
        }),
      x_center = 
        x_train |>
        purrr::map(attr, which = "scaled:center"),
      x_scale = 
        x_train |>
        purrr::map(attr, which = "scaled:scale"),
      prediction_domain =
        purrr::map(seq_along(id_model), function(i) {
          x[[i]] |>
            dplyr::mutate(
              spectra =
                purrr::map(spectra, function(.x) {
                  .x |>
                    dplyr::mutate(
                      y = (y  - x_center[[i]]) / x_scale[[i]]
                    )
                })
            ) |>
            irpeat::irp_as_irp_prediction_domain()
        }),
      likelihood_name = "beta",
      likelihood = list(brms::Beta(link = "logit", link_phi = "log")),
      # prior_hs_par_ratio =
      #   purrr::map_dbl(x_train, function(.x) {
      #     n <- nrow(.x) - 4L # in each case, four observations are removed
      #     p <- ncol(.x) - 2L # two target variables have to be subtracted
      #    p0 <- 5L # prior guess for the number of relevant variables
      #     p0/(p-p0) * 1/sqrt(n) # tau0
      #   }),
      par_ratio =
        purrr::map_dbl(x_train, function(.x) {
          p <- ncol(.x)
          p0 <- 5L
          p0/p
        }),
      priors =
        purrr::map(id_model, function(i) {
          c(
            brms::prior_string(paste0("horseshoe(df = 1, par_ratio = ", par_ratio[[i]], ", df_global = 1, autoscale = TRUE)"),
                               class = "b"),
            brms::prior_string("normal(0.0, 2.5)", class = "Intercept"),
            brms::prior_string("gamma(7.0, 7.0/90.0)", class = "phi", lb = 0.0)
          )
        })
    )
  
}




#### Compute Models ####

#' Helper function that creates a data frame in the format expected by brms
#'
#' @param dd_model_info One row from `dd_model_info`.
#' 
#' @export
dd_make_brms_data <- function(dd_model_info) {
  
  stopifnot(nrow(dd_model_info) == 1L)
  
  tibble::tibble(
    y = dd_model_info$y[[1]]$mass_relative_mass,
    x = dd_model_info$x_train[[1]]
  )
  
}


#' Compiles a brmsfit model (to avoid recompilation with different datasets)
#'
#' @param dd_model_info One row from `dd_model_info`.
#'
#' @param ... Additional arguments passed to `brms::brm()`.
#' 
#' @export
dd_make_brms_compiled_model <- function(dd_model_info, ...) {
  
  d_modeling <- dd_model_info
  
  brms::brm(
    y ~ ., 
    data = dd_make_brms_data(dd_model_info = d_modeling), 
    family = d_modeling$likelihood[[1]],
    prior = d_modeling$priors[[1]],
    chains = 0,
    ...
  )
  
}


#' Defines brms MCMC settings
#'
#' @export
dd_get_stan_1_mcmc_settings <- function() {
  
  list(
    iter = 4000L,
    warmup = 2000L,
    chains = 4L,
    adapt_delta = 0.99,
    max_treedepth = 12
  )
  
}


#### For cross-validation ####

#' Defines grouping of observations for grouped k-fold CV
#' 
#' @param K Integer. The number of CV folds to create.
#' 
#' @export
dd_make_cv_groups <- function(dd_data_model, K = 10L) {
  
  res <- 
  dd_data_model |>
    dplyr::mutate(
      reuter2019_site =
        dplyr::case_when(
          dd_dataset_label == "reuter2019" ~ stringr::str_extract(sample_label, "^[A-Z]"),
          TRUE ~ NA_character_
        ),
      group_cv =
        paste0(dd_dataset_label, "_", id_dataset, "_", taxon_rank_value, "_", reuter2019_site) |>
        as.factor()
    ) |>
    dplyr::pull(group_cv)
  
  loo::kfold_split_grouped(K = K, x = res)
  
}

#' Defines grouping of observations for stratified k-fold CV
#' 
#' Same as `dd_make_cv_groups`, but uses `loo::kfold_split_stratified()` instead
#' of `loo::kfold_split_grouped()`.
#' 
#' @param K Integer. The number of CV folds to create.
#' 
#' @export
dd_make_cv_groups_2 <- function(dd_data_model, K = 10L) {
  
  res <- 
    dd_data_model |>
    dplyr::mutate(
      reuter2019_site =
        dplyr::case_when(
          dd_dataset_label == "reuter2019" ~ stringr::str_extract(sample_label, "^[A-Z]"),
          TRUE ~ NA_character_
        ),
      group_cv =
        paste0(dd_dataset_label, "_", id_dataset, "_", taxon_rank_value, "_", reuter2019_site) |>
        as.factor()
    ) |>
    dplyr::pull(group_cv)
  
  loo::kfold_split_stratified(K = K, x = res)
  
}


#' Computes the RMSE for cross-validation models
#' 
#' @param A list with the same structure as `dd_stan_1_kfold`.
#' 
#' @export
dd_make_kfold_rmse <- function(kfold, dd_model_info, file) {
  
    purrr::map(seq_along(kfold), function(i) {
      res <- 
        tibble::tibble(
          yhat = 
            brms::kfold_predict(kfold[[i]])$yrep |>
            posterior::as_draws_df() |>
            posterior::as_draws_rvars()
        ) |>
        dplyr::mutate(
          yhat = do.call("c", yhat),
          rmse = dd_rmse_rvar(yhat = yhat, y = dd_model_info$y[[i]]$mass_relative_mass),
          id_model = i
        ) |>
        dplyr::select(id_model, rmse) |>
        dplyr::slice(1)
    }) |>
    dplyr::bind_rows() |>
    saveRDS_rvars(file = file)
  
  file
  
}

#' @describeIn dd_make_kfold_rmse
#' @export
dd_make_fit_rmse <- function(fit, dd_model_info, file) {
  
  purrr::map(seq_along(fit), function(i) {
    res <- 
      tibble::tibble(
        yhat = 
          brms::posterior_predict(fit[[i]]) |>
          posterior::as_draws_df() |>
          posterior::as_draws_rvars()
      ) |>
      dplyr::mutate(
        yhat = do.call("c", yhat),
        rmse = dd_rmse_rvar(yhat = yhat, y = dd_model_info$y[[i]]$mass_relative_mass),
        id_model = i
      ) |>
      dplyr::select(id_model, rmse) |>
      dplyr::slice(1)
  }) |>
    dplyr::bind_rows() |>
    saveRDS_rvars(file = file)
  
  file
  
}


#' Creates a data frame that collects fitted values and predictions for the training data
#' 
#' @export
dd_make_data_model_evaluation_1 <- function(dd_data_model, dd_stan_1_fit, dd_stan_1_kfold, dd_stan_1_kfold_2, dd_stan_1_kfold_cv_folds, dd_stan_1_kfold_cv_folds_2, dd_model_info) {
  
  # helper function to create data frames for cross-validation predictions
  dd_make_plot_15_helper_1 <- function(dd_model_info, kfold_cv_folds, kfold) {
    res <- 
      purrr::map(seq_len(nrow(dd_model_info)), function(i) {
        dplyr::bind_cols(
          dd_data_model |>
            dplyr::select(-spectra) |>
            dplyr::mutate(
              group_cv = kfold_cv_folds
            ),
          tibble::tibble(
            y = dd_model_info$y[[i]]$mass_relative_mass,
            yhat = 
              brms::kfold_predict(kfold[[i]])$yrep |>
              as.data.frame() |>
              purrr::map(function(.x) {
                tibble::tibble(
                  yhat_mean = mean(.x),
                  yhat_lower = posterior::quantile2(.x, probs = 0.025),
                  yhat_upper = posterior::quantile2(.x, probs = 0.975)
                )
              }),
            id_model = i
          ) |>
            tidyr::unnest(yhat)
        )
      })
    do.call("rbind", res)
  }
  
  dd_make_plot_15_helper_2 <- function(dd_model_info, dd_stan_1_fit) {
    res <- 
      purrr::map(seq_len(nrow(dd_model_info)), function(i) {
        dplyr::bind_cols(
          dd_data_model |>
            dplyr::select(-spectra) |>
            dplyr::mutate(
              group_cv = NA_integer_
            ),
          tibble::tibble(
            y = dd_model_info$y[[i]]$mass_relative_mass,
            yhat = 
              brms::posterior_predict(dd_stan_1_fit[[i]]) |>
              as.data.frame() |>
              purrr::map(function(.x) {
                tibble::tibble(
                  yhat_mean = mean(.x),
                  yhat_lower = posterior::quantile2(.x, probs = 0.025),
                  yhat_upper = posterior::quantile2(.x, probs = 0.975)
                )
              }),
            id_model = i
          ) |>
            tidyr::unnest(yhat)
        )
      })
    do.call("rbind", res)
  }
  
  res_1 <- 
    dd_make_plot_15_helper_1(
      dd_model_info = dd_model_info, 
      kfold_cv_folds = dd_stan_1_kfold_cv_folds, 
      kfold = dd_stan_1_kfold
    ) |>
    dplyr::mutate(
      cv_type = "Grouped CV"
    )
  
  res_2 <- 
    dd_make_plot_15_helper_1(
      dd_model_info = dd_model_info, 
      kfold_cv_folds = dd_stan_1_kfold_cv_folds_2, 
      kfold = dd_stan_1_kfold_2
    ) |>
    dplyr::mutate(
      cv_type = "Stratified CV"
    )
  
  res_fit <- 
    dd_make_plot_15_helper_2(
      dd_model_info = dd_model_info, 
      dd_stan_1_fit = dd_stan_1_fit
    ) |>
    dplyr::mutate(
      cv_type = "Fitted values"
    )
  
  do.call("rbind", list(res_1, res_2, res_fit)) 
  
}


#' Summarize information on Monte Carlo errors for target quantities
#' 
#' @export
dd_get_mcse_1 <- function(dd_stan_1_fit, dd_mir_preprocessing_config, dd_data_model, dd_data_pmird_peat_cores, dd_model_info) {
  
  # helper function to compute MCSE
  dd_get_mcse_1_helper_1 <- function(x, summarize = FALSE, na_rm = FALSE) {
    res <- 
      x |>
      dplyr::mutate(
        dd_mcse_mean = posterior::mcse_mean(degree_of_decomposition),
        dd_mcse_sd = posterior::mcse_sd(degree_of_decomposition),
        dd_mcse_lower = posterior::mcse_quantile(degree_of_decomposition, probs = 0.025),
        dd_mcse_upper = posterior::mcse_quantile(degree_of_decomposition, probs = 0.975)
      )
    
    if(summarize) {
      res <- 
        res |>
        dplyr::group_by(degree_of_decomposition_in_pd) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::starts_with("dd_mcse_"),
            function(.x) {
              tibble::tibble(
                min = min(.x, na.rm = na_rm),
                max = max(.x, na.rm = na_rm)
              ) |>
                list()
            }
          ),
          .groups = "drop"
        ) |>
        tidyr::pivot_longer(
          dplyr::starts_with("dd_mcse_"),
          names_to = "variable",
          values_to = "mcse"
        ) |>
        tidyr::unnest("mcse")
    }
    
    res
    
  }
  
  purrr::map(seq_along(dd_stan_1_fit), function(i) {
    res <- 
      list(
        irp_degree_of_decomposition(
          x = dd_data_model, 
          model = dd_stan_1_fit[[i]], 
          config = dd_mir_preprocessing_config[[i]], 
          prediction_domain = dd_model_info$prediction_domain[[i]]
        ) |>
          dd_get_mcse_1_helper_1(summarize = TRUE) |>
          dplyr::mutate(
            dataset_type = "training_data"
          ),
        irp_degree_of_decomposition(
          x = dd_data_pmird_peat_cores, 
          model = dd_stan_1_fit[[i]], 
          config = dd_mir_preprocessing_config[[i]], 
          prediction_domain = dd_model_info$prediction_domain[[i]]
        ) |>
          dd_get_mcse_1_helper_1(summarize = TRUE, na_rm = TRUE) |>
          dplyr::mutate(
            dataset_type = "peat_data"
          )
      ) |> 
      dd_do_call("rbind") |>
      dplyr::mutate(
        id_model = i
      ) |>
      dplyr::relocate("id_model", .before = dplyr::everything())
    res
  }) |>
    dd_do_call("rbind")
  
}


#' Extracts the Rhat for all parameters of a list of `brmsfit` objects
#' 
#' @export
dd_get_rhat_1 <- function(dd_stan_1_fit) {
  
  purrr::map(seq_along(dd_stan_1_fit), function(i) {
    res <- 
      brms::rhat(dd_stan_1_fit[[i]])
    tibble::tibble(
      id_model = i,
      parameter = names(res),
      rhat = res
    )
  }) |>
    dd_do_call("rbind")
  
}


#### For projpred ####

#' Defines which models to compute and which settings to use
#' 
#' @export
dd_get_model_info_2 <- function(dd_data_model_preprocessed, dd_mir_preprocessing_settings) {
  
  d <- dd_data_model_preprocessed
  d_spectral_preprocessing_settings <- dd_mir_preprocessing_settings
  
  tidyr::expand_grid(
    id_preprocessing = d_spectral_preprocessing_settings$id_preprocessing, # compute one model for each preprocessing variant
  ) |>
    dplyr::mutate(
      id_model = seq_along(id_preprocessing),
      id_sample_all = 
        purrr::map(id_model, function(i) {
          d[[id_preprocessing[[i]]]] |>
            dplyr::pull(id_sample)
        }),
      sample_size_all = purrr::map_int(id_sample_all, length),
      do_compute_model = TRUE,
      validation_mode = "loo",
      y = 
        purrr::map(seq_along(id_model), function(i) {
          .id_sample_all <- id_sample_all[[i]]
          .id_preprocessing <- id_preprocessing[[i]]
          d[[.id_preprocessing]] |>
            dplyr::filter(id_sample %in% .id_sample_all) |>
            dplyr::select(dplyr::all_of(c("mass_relative_mass"))) |>
            dplyr::mutate(
              mass_relative_mass =
                dplyr::case_when(
                  mass_relative_mass == 1 ~ mass_relative_mass - 0.001, #---note beta distribution assumes values are in (0,1)
                  TRUE ~ mass_relative_mass
                ),
              mass_relative_mass = 1 - mass_relative_mass,
              mass_relative_mass = brms::logit_scaled(mass_relative_mass)
            )
        }),
      x = 
        purrr::map(seq_along(id_model), function(i) {
          .id_sample_all <- id_sample_all[[i]]
          .id_preprocessing <- id_preprocessing[[i]]
          d[[.id_preprocessing]] |>
            dplyr::filter(id_sample %in% .id_sample_all) |>
            dplyr::select(id_sample, spectra)
        }),
      y_center = 0,
      y_scale = 1,
      x_train =
        purrr::map(x, function(.x) {
          .x |>
            ir::ir_flatten() |>
            dplyr::select(-1) |>
            t() |>
            scale(center = TRUE, scale = TRUE)
        }),
      x_center = 
        x_train |>
        purrr::map(attr, which = "scaled:center"),
      x_scale = 
        x_train |>
        purrr::map(attr, which = "scaled:scale"),
      prediction_domain =
        purrr::map(seq_along(id_model), function(i) {
          x[[i]] |>
            dplyr::mutate(
              spectra =
                purrr::map(spectra, function(.x) {
                  .x |>
                    dplyr::mutate(
                      y = (y  - x_center[[i]]) / x_scale[[i]]
                    )
                })
            ) |>
            irpeat::irp_as_irp_prediction_domain()
        }),
      likelihood_name = "Gaussian",
      likelihood = list(gaussian(link = "identity")),
      # prior_hs_par_ratio =
      #   purrr::map_dbl(x_train, function(.x) {
      #     n <- nrow(.x) - 4L # in each case, four observations are removed
      #     p <- ncol(.x) - 2L # two target variables have to be subtracted
      #    p0 <- 5L # prior guess for the number of relevant variables
      #     p0/(p-p0) * 1/sqrt(n) # tau0
      #   }),
      par_ratio =
        purrr::map_dbl(x_train, function(.x) {
          p <- ncol(.x)
          p0 <- 5L
          p0/p
        }),
      priors =
        purrr::map(id_model, function(i) {
          c(
            brms::prior_string(paste0("horseshoe(df = 1, par_ratio = ", par_ratio[[i]], ", df_global = 1, autoscale = TRUE)"),
                               class = "b"),
            brms::prior_string("normal(0.0, 2.5)", class = "Intercept"),
            brms::prior_string("normal(0.0, 0.5)", class = "sigma", lb = 0.0)
          )
        })
    )
  
}

#' Helper function that creates a data frame in the format expected by brms
#'
#' @param dd_model_info One row from `dd_model_info`.
#' 
#' @export
dd_make_brms_data_2 <- function(dd_model_info) {
  
  stopifnot(nrow(dd_model_info) == 1L)
  
  tibble::tibble(
    y = dd_model_info$y[[1]]$mass_relative_mass
  ) |> 
    dplyr::bind_cols(
      as.data.frame(dd_model_info$x_train[[1]])
    )
  
}


#' Compiles a brmsfit model (to avoid recompilation with different datasets)
#'
#' @param dd_model_info One row from `dd_model_info`.
#'
#' @param ... Additional arguments passed to `brms::brm()`.
#' 
#' @export
dd_make_brms_compiled_model_2 <- function(dd_model_info, ...) {
  
  d_modeling <- dd_model_info
  
  brms::brm(
    y ~ ., 
    data = dd_make_brms_data_2(dd_model_info = d_modeling), 
    family = d_modeling$likelihood[[1]],
    prior = d_modeling$priors[[1]],
    chains = 0,
    ...
  )
  
}

