# Functions to estimate gamma_MIRS in multi component systems when gamma for each component is known

#### utils ####

#' Function factory for prediction functions from project eb1149 (degree of decompsoition). Predicts gamma_MIRS on the latent scale (logit-transformed average)
#'
#' @param model A [`brmsfit`](brms::brm) object.
#'
#' @param config A list with configuration parameters for the prediction.
#'
#' @param prediction_domain A list with two elements:
#' \describe{
#'   \item{`train`}{An
#'   [`irp_prediction_domain`](irpeat::new_irp_prediction_domain) object
#'   representing the prediction domain for the training data.}
#'   \item{`test`}{An
#'   [`irp_prediction_domain`](irpeat::new_irp_prediction_domain) object
#'   representing the prediction domain for the testing data.}
#' }
#'
#' @param target_variable_name A character value representing the column name
#' for the column with predicted values.
#'
#' @param irpeatmodels_required_version A character value with format "x.y.z"
#' representing the minimum version of the 'irpeatmodels' package required
#' to make predictions with the function which is generated.
#'
#' @param .f_check_packages A function which checks whether the correct version
#' of 'irpeatmodels' is installed and of other packages which may be required.
#' Must take the argument `irpeatmodels_required_version` as input.
#'
#' @return A function that makes predictions with a model computed in
#' project eb1079.
#'
#' @keywords Internal
#' @noRd
irp_function_factory_eb1149 <- function(model, config, prediction_domain, target_variable_name, irpeatmodels_required_version = "0.0.0") {
  
  .f_check_packages <- function() {
    
    irpeat:::check_irpeatmodels(version = irpeatmodels_required_version)
    if(! requireNamespace("brms", versionCheck = list(op = ">=", version = "2.20.1"), quietly = TRUE)) {
      rlang::abort(paste0("You have to install the 'brms' package (>=", version,") to use this function."))
    }
    if(! requireNamespace("brms", versionCheck = list(op = ">=", version = "2.20.1"), quietly = TRUE)) {
      rlang::abort(paste0("You have to install the 'brms' package (>=", version,") to use this function."))
    }
    if(! requireNamespace("posterior", versionCheck = list(op = ">=", version = "1.5.0"), quietly = TRUE)) {
      rlang::abort(paste0("You have to install the 'brms' package (>=", version,") to use this function."))
    }
    
  }
  
  function(x, do_summary = FALSE, summary_function_mean = mean, summary_function_sd = stats::sd, check_prediction_domain = "train", return_as_list = FALSE) {
    
    .f_check_packages()
    stopifnot(inherits(x, "ir"))
    stopifnot(is.logical(do_summary) && length(do_summary) == 1)
    stopifnot(is.logical(return_as_list) && length(return_as_list) == 1)
    if(do_summary && return_as_list) {
      stop("Both `do_summary` and `return_as_list` are set to `TRUE`, but only one of both must be `TRUE`.")
    }
    
    x_or <- x
    x <- irpeat:::irp_preprocess_eb1149(x = x, config = config)
    x_in_pd <- irpeat::irp_is_in_prediction_domain(x = x, prediction_domain = prediction_domain)
    
    newdata <-
      tibble::tibble(
        x =
          x |>
          ir::ir_flatten() |>
          dplyr::select(-1) |>
          t()
      )
    
    yhat <-
      brms::posterior_linpred(object = model, newdata = newdata) |>
      posterior::rvar()
    
    posterior::draws_of(yhat) <- units::set_units(posterior::draws_of(yhat), value = "1", mode = "standard")
    
    if(return_as_list) {
      yhat <-
        posterior::draws_of(yhat) |>
        as.data.frame() |>
        purrr::map(function(.x) .x)
    } else if(do_summary) {
      yhat <-
        purrr::map(yhat, function(.x) {
          quantities::set_quantities(
            x = summary_function_mean(.x),
            unit = "1",
            errors = summary_function_sd(.x),
            mode = "standard"
          )
        })
      yhat <- do.call("c", yhat)
    }
    
    b_Intercept <- 
      as.data.frame(model) |> 
      dplyr::pull(b_Intercept) |> 
      posterior::rvar()
    if(do_summary) {
      b_Intercept <- mean(b_Intercept)
    }
    
    
    x_or[[target_variable_name]] <- structure(yhat, b_Intercept = b_Intercept)
    x_or[[paste0(target_variable_name, "_in_pd")]] <- x_in_pd$is_in_prediction_domain
    
    x_or
    
  }
  
}

#' Preprocesses MIRS to compute mixtures
#' 
#' @param x An `ir` object.
#' 
#' @return `x` in preprocessed form.
#' 
#' @export
dd_preprocess_for_mixing <- function(x) {
  x |>
    ir::ir_interpolate() |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
} 


#' Computes properties of multi-component systems from individual components
#' 
#' @param x A list of `ir` objects, where each row is a litter type. Each data 
#' frame needs to have the columns `scale_factor` and `gamma`.
#' 
#' @export
dd_make_mixture <- function(x) {
  
  stopifnot(is.list(x))
  stopifnot(all(purrr::map_lgl(x, function(.x) inherits(.x, "ir"))))
  stopifnot(all(purrr::map_lgl(x, function(.x) all(c("scale_factor", "gamma") %in% colnames(.x)))))
  
  # preprocess spectra for mixing
  x <- 
    purrr::map(x, dd_preprocess_for_mixing) |>
    purrr::map(function(.x) {
      .x * .x$scale_factor
    })
  
  x <- 
    purrr::map2(x, seq_along(x), function(.x, .y) {
      .x |>
        dplyr::mutate(
          id_mixture = seq_along(gamma),
          id_component = .y,
          scale_factor_0 = scale_factor / (1 - gamma)
        )
    })
  
  # summarize parameters of components
  res_components_properties <- 
    purrr::map(x, function(.x) {
      .x |>
        dplyr::select(id_mixture, id_component, id_sample, taxon_rank_value, taxon_organ, scale_factor, scale_factor_0, gamma, spectra)
    }) |>
    dplyr::bind_rows()
  
  res_mixture <- 
    purrr::reduce(x, .f = magrittr::add) |>
    dplyr::select(id_mixture, spectra)
  
  res_mixture_properties <- 
    res_components_properties |>
    dplyr::group_by(id_mixture) |>
    dplyr::summarise(
      m_1 = sum(scale_factor),
      m_0 = sum(scale_factor_0),
      gamma = 1 - sum(scale_factor) / sum(scale_factor_0),
      .groups = "drop"
    )
  
  list(
    components = res_components_properties,
    mixture = 
      dplyr::left_join(
        res_mixture, 
        res_mixture_properties,
        by = "id_mixture"
      )
  )
  
}



#### prediction functions ####

#' @export
irp_degree_of_decomposition_1_logit <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_1_brms,
    config = irpeatmodels::model_degree_of_decomposition_1_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_1_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_1_logit",
    irpeatmodels_required_version =
      "0.0.0"
  )

#' @export
irp_degree_of_decomposition_2_logit <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_2_brms,
    config = irpeatmodels::model_degree_of_decomposition_2_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_2_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_2_logit",
    irpeatmodels_required_version =
      "0.0.0"
  )


#' @export
irp_degree_of_decomposition_3_logit <-
  irp_function_factory_eb1149(
    model = irpeatmodels::model_degree_of_decomposition_3_brms,
    config = irpeatmodels::model_degree_of_decomposition_3_config,
    prediction_domain = irpeatmodels::model_degree_of_decomposition_3_prediction_domain,
    target_variable_name =
      "degree_of_decomposition_3_logit",
    irpeatmodels_required_version =
      "0.0.0"
  )



#### Estimating gamma_MIRS from gamma ####

#' Model that predicts gamma_MIRS from gamma for individual components
#' 
#' @export
dd_make_model_gamma_mirs_from_gamma <- function(dd_data_model, id_model) {
  
  prediction_function <- sym(paste0("irp_degree_of_decomposition_", id_model, "_logit"))
  target_variable <- paste0("degree_of_decomposition_", id_model, "_logit")
  
  d_litter <- 
    dd_data_model |>
    eval(prediction_function)(
      do_summary = TRUE, 
      summary_function_mean = mean, 
      summary_function_sd = posterior::sd
    )
  
  d_litter <- 
    d_litter |>
    dplyr::mutate(
      b_Intercept = attr(d_litter[, as.character(target_variable), drop = TRUE], "b_Intercept"),
      degree_of_decomposition_no_intercept = as.numeric(get(rlang::as_name(target_variable))) - b_Intercept,
      gamma = 1 - mass_relative_mass,
      gamma = gamma + 0.001
    )
  
  d_litter |>
    dplyr::select(degree_of_decomposition_no_intercept, gamma, b_Intercept) %>%
    structure(class = c("gamma_model", class(.)))
  
}

#' Predicts gamma_MIRS from gamma
#' 
#' @export
dd_predict_with_model_gamma_mirs_from_gamma <- function(object, ..., newdata = NULL) {
  
  if(is.null(newdata)) {
    newdata <- object
  }
  
  with(object, brms::logit_scaled(newdata$gamma) - b_Intercept)
  
}

#' For a mixture, predicts gamma_MIRS from the known gamma of the components
#' 
#' @param x An object as returned by `dd_make_mixture()`.
#' 
#' @param dd_model_gamma_mirs_from_gamma_1 One model from `dd_model_gamma_mirs_from_gamma_1`.
#' 
#' @export
dd_mixture_estimate_gamma_mirs_from_gamma <- function(x, dd_model_gamma_mirs_from_gamma_1, b_Intercept) {
  
  res <- 
    x$components |>
    dplyr::mutate(
      gamma_mirs_hat = 
        dd_predict_with_model_gamma_mirs_from_gamma(
          dd_model_gamma_mirs_from_gamma_1,
          newdata = x$components
        )
    ) |>
    dplyr::group_by(id_mixture) |>
    dplyr::summarise(
      gamma_mirs_hat = sum(gamma_mirs_hat * scale_factor/sum(scale_factor)) + b_Intercept,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      gamma_mirs_hat = brms::inv_logit_scaled(gamma_mirs_hat)
    )
  
  res |>
    dplyr::pull(gamma_mirs_hat)
  
}



#' Preprocesses data from the pmird database for peat cores for application of the mixing model
#' 
#' @export
dd_make_data_preparation_for_mixing_model_1 <- function(dd_data_pmird_peat_cores_macrofossils, dd_data_pmird_peat_cores_irpeat) {
  
  d_macrofossil <- 
    dd_data_pmird_peat_cores_macrofossils |> 
    dplyr::filter(macrofossil_volume_fraction > 0)
  
  d_macrofossil_summary <- 
    d_macrofossil |>
    dplyr::group_by(core_label, id_sample, sample_depth_upper, sample_depth_lower) |>
    dplyr::summarise(
      n_macrofossil_types = length(taxon_rank_value),
      sum_macrofossil_volume_fraction = sum(macrofossil_volume_fraction),
      .groups = "drop"
    )
  
  d_macrofossil <- 
    d_macrofossil |>
    dplyr::left_join(
      d_macrofossil_summary |>
        dplyr::select(id_sample, sum_macrofossil_volume_fraction),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      w = macrofossil_volume_fraction / sum_macrofossil_volume_fraction
    )
  
  d_peat <- 
    d_macrofossil_summary |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat |>
        dplyr::select(core_label, sample_depth_upper, sample_depth_lower, bulk_density, bulk_density_1, spectra),
      by = dplyr::join_by(core_label, within(y$sample_depth_lower, y$sample_depth_upper, x$sample_depth_lower, x$sample_depth_upper))
    ) |>
    dplyr::mutate(
      spectra_is_empty =
        purrr::map_lgl(spectra, is.null)
    ) |>
    dplyr::filter(! spectra_is_empty) |>
    ir::ir_as_ir() |>
    irpeat::irp_degree_of_decomposition_1(do_summary = TRUE, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_2(do_summary = TRUE, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_3(do_summary = TRUE, summary_function_sd = posterior::sd)
  
  d_macrofossil <- 
    d_macrofossil |>
    dplyr::filter(id_sample %in% d_peat$id_sample)
  
  d_macrofossil_summary <- 
    d_macrofossil_summary |>
    dplyr::filter(id_sample %in% d_peat$id_sample)
  
  list(
    d_macrofossil = d_macrofossil,
    d_macrofossil_summary = d_macrofossil_summary,
    d_peat = d_peat
  )
  
}

#' Helper function to compute phi for a beta distribution from mu and sd
compute_beta_phi <- function(mean, sd) {
  # Check for valid input
  if (any(mean <= 0) || any(mean >= 1)) {
    stop("Mean must be between 0 and 1 (exclusive).")
  }
  if (any(sd <= 0)) {
    stop("Standard deviation must be positive.")
  }
  
  # Compute the concentration parameter (common factor)
  var <- sd^2
  common_factor <- (mean * (1 - mean) / var) - 1
  
  # Compute alpha and beta
  alpha <- mean * common_factor
  beta <- (1 - mean) * common_factor
  
  # Compute phi (total concentration)
  phi <- alpha + beta
  
  phi
}


#' Prepares Stan data for the mixing model
#' 
#' @export
dd_make_stan_data_mixing_model_1 <- function(dd_data_preparation_for_mixing_model_1, id_model) {
  
  gamma_mirs_model <- 
    switch(
      id_model,
      "1" = irpeatmodels::model_degree_of_decomposition_1_brms,
      "2" = irpeatmodels::model_degree_of_decomposition_2_brms,
      "3" = irpeatmodels::model_degree_of_decomposition_3_brms
    )
 
  b_Intercept <- 
    as.data.frame(gamma_mirs_model) |> 
    dplyr::pull(b_Intercept)
  
  d <- dd_data_preparation_for_mixing_model_1
  
  index_components <- 
    d$d_macrofossil_summary |>
    dplyr::select(n_macrofossil_types) |>
    dplyr::mutate(
      start = c(1, cumsum(n_macrofossil_types[-length(n_macrofossil_types)]) + 1),
      end = cumsum(n_macrofossil_types)
    )
  
  gamma_mirs_obs <- d$d_peat[[paste0("degree_of_decomposition_", id_model)]]
  gamma_mirs_obs_mu <- as.numeric(gamma_mirs_obs)
  gamma_mirs_obs_sd <- as.numeric(errors::errors(gamma_mirs_obs))
  
  
  stan_data <- 
    dplyr::lst(
      N = nrow(d$d_peat),
      K_total = nrow(d$d_macrofossil),
      w = d$d_macrofossil$w,
      gamma_mirs_obs_p1 = gamma_mirs_obs_mu,
      gamma_mirs_obs_p2 = compute_beta_phi(gamma_mirs_obs_mu, gamma_mirs_obs_sd),
      b_intercept_p1 = mean(b_Intercept),
      b_intercept_p2 = sd(b_Intercept),
      index_components = 
        index_components |>
        as.matrix(),
      D_p1 = rep((0.5) * 2, K_total),
      D_p2 = rep((1 - 0.5) * 2, K_total),
      phi_p1 = rep(5, N),
      phi_p2 = rep(5/200, N),
      phi_p3 = rep(200, N)
    )
  
  stan_data
  
}


#' Fits the Stan mixing model
#' 
#' @export
dd_make_fit_stan_mixing_model_1 <- function(dd_mm2_cmd, dd_stan_data_mixing_model_1, file) {
  
  stan_data <- dd_stan_data_mixing_model_1
  
  stan_fit <- 
    dd_mm2_cmd$sample(
      data = tidybayes::compose_data(stan_data),
      iter_sampling = 1000,
      iter_warmup = 2000,
      chains = 4,
      parallel_chains = 1,
      seed = 646,
      sig_figs = 10,
      max_treedepth = 16,
      adapt_delta = 0.8,
      threads_per_chain = 1,
      save_warmup = TRUE,
      refresh = 100
    )
  
  stan_fit$save_object(file)
  file
  
}


#' Adds predicted gamma for components to `dd_data_preparation_for_mixing_model_1` and reconstruct NPP (initial mass rates)
#'
#' @export
dd_make_res_reconstructions_1 <- function(dd_data_preparation_for_mixing_model_1, dd_fit_stan_mixing_model_1, file = "mixing_model/dd_res_reconstructions.rds") {
  
  res <- dd_data_preparation_for_mixing_model_1
  
  res$d_peat$bulk_density_1_summarized <- res$d_peat$bulk_density_1
  res$d_peat <- 
    res$d_peat |>
    irpeat::irp_bulk_density_1(do_summary = FALSE, return_as_list = FALSE)
  a <- as.data.frame(res$d_peat)
  index <- seq_len(nrow(posterior::draws_of(a$bulk_density_1[[1]])))
  index <- posterior::rvar(x = index %in% sample(index, size = 4000, replace = FALSE))
  a$bulk_density_1 <- a$bulk_density_1[index]
  res$d_peat$bulk_density_1 <- a$bulk_density_1
  
  res_gamma <- 
    dd_fit_stan_mixing_model_1 |>
    purrr::map2(seq_along(dd_fit_stan_mixing_model_1), function(.x, .y) {
      .x |>
        tidybayes::spread_rvars(
          D[id_component]
        ) |>
        dplyr::select(D) |>
        setNames(nm = paste0("component_gamma_", .y))
    }) |>
    dplyr::bind_cols()
    
  res$d_macrofossil <- 
    res$d_macrofossil |>
    dplyr::bind_cols(res_gamma) |>
    dplyr::left_join(
      res$d_peat |>
        dplyr::select(id_sample, bulk_density_1),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      bulk_density_1 = bulk_density_1 * 10000/1000, #---note: ignore uncertainty for now
      component_m_1 = w * bulk_density_1 * (sample_depth_lower - sample_depth_upper),
      dplyr::across(
        dplyr::starts_with("component_gamma_"),
        function(.x) {
          component_m_1/(1 - .x)
        },
        .names = "component_m_011{.col}"
      )
    ) |>
    dplyr::rename_with(
      .fn = function(.x) stringr::str_remove(.x, pattern = "11component_gamma"),
      .cols = dplyr::starts_with("component_m_011component_gamma_")
    )
  
  res$d_macrofossil_summary <- 
    res$d_macrofossil_summary |>
    dplyr::left_join(
      res$d_macrofossil |>
        dplyr::select(id_sample, dplyr::starts_with("component_gamma_"), dplyr::starts_with("component_m_0_")) |>
        dplyr::group_by(id_sample) |>
        dplyr::summarize(
          dplyr::across(
            dplyr::starts_with("component_m_0_"),
            function(.x) {
              posterior::rvar_sum(.x)
            },
            .names = "total_{.col}"
          ),
          dplyr::across(
            dplyr::starts_with("component_gamma_"),
            function(.x) {
              cur_id <- 
                dplyr::cur_column() |>
                stringr::str_extract("\\d{1}$")
              m_0 <- dplyr::cur_data()[[paste0("component_m_0_", cur_id)]]
              m_0_sum <- dplyr::cur_data()[[paste0("total_component_m_0_", cur_id)]]
              posterior::rvar_sum(.x * m_0/m_0_sum)
            }
          ),
          .groups = "drop"
        ) |>
        dplyr::rename_with(
          .fn = function(.x) stringr::str_remove(.x, pattern = "^component_"),
          .cols = dplyr::starts_with("component_gamma_")
        ) |>
        dplyr::rename_with(
          .fn = function(.x) stringr::str_remove(.x, pattern = "^total_component_"),
          .cols = dplyr::starts_with("total_component_m_0_")
        ),
      by = "id_sample"
    )
  
  # initial mass rate
  res$d_macrofossil <- 
    res$d_macrofossil |>
    dplyr::group_split(core_label) |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          age =
            tibble::tibble(
              age_upper = suppressMessages(dd_get_age_for_depth(sample_depth_upper, core_label = unique(core_label), postbomb = 1, ssize = 4000, bacon_is_loaded = FALSE)),
              age_lower = suppressMessages(dd_get_age_for_depth(sample_depth_lower, core_label = unique(core_label), postbomb = 1, ssize = 4000, bacon_is_loaded = TRUE))
            ) |>
            dplyr::mutate(
              dplyr::across(dplyr::everything(), function(.x) (lubridate::year(sampling_date) - 1950) + .x)
            )
        ) |>
        tidyr::unnest("age")
    }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      duration = age_lower - age_upper,
      dplyr::across(
        dplyr::starts_with("component_m_0_"),
        function(.x) {
          .x/duration
        },
        .names = "npp_{.col}"
      )
    ) |>
    dplyr::rename_with(
      .fn = function(.x) stringr::str_remove(.x, pattern = "component_m_0_"),
      .cols = dplyr::starts_with("npp_component_m_0_")
    )
  
  res$d_macrofossil_summary <- 
    res$d_macrofossil_summary |>
    dplyr::left_join(
      res$d_macrofossil |>
        dplyr::select(id_sample, age_upper, age_lower, duration) |>
        dplyr::filter(! duplicated(id_sample)),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("m_0_"),
        function(.x) {
          .x/duration
        },
        .names = "npp_{.col}"
      )
    ) |>
    dplyr::rename_with(
      .fn = function(.x) stringr::str_remove(.x, pattern = "m_0_"),
      .cols = dplyr::starts_with("npp_m_0_")
    )
  
  # to analyze the relation between saturated hydraulic conductivity and degree of decomposition
  res$d_macrofossil_summary <- 
    res$d_macrofossil_summary |>
    dplyr::left_join(
      res$d_peat |>
        irpeat::irp_saturated_hydraulic_conductivity_1(do_summary = FALSE) |>
        irpeat::irp_volume_fraction_solids_1(do_summary = FALSE) |>
        dplyr::select(core_label, id_sample, saturated_hydraulic_conductivity_1, saturated_hydraulic_conductivity_1_in_pd, volume_fraction_solids_1, volume_fraction_solids_1_in_pd),
      by = c("core_label", "id_sample")
    )
  
  saveRDS_rvars(res, file) 
  file
  
} 



 