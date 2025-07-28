#' Defines settings for the spectral preprocessing for this project
#' 
#' @export
dd_get_mir_preprocessing_settings <- function() {
  
  spectral_preprocessing_n <- 3L
  
  d_spectral_preprocessing_settings <- 
    tibble::tibble(
      do_interpolate = TRUE,
      interpolate_start = rep(list(NULL), spectral_preprocessing_n),
      interpolate_dw = 1,
      do_clip = TRUE,
      clip_range = list(tibble::tibble(start = 650, end = 4000)),
      do_interpolate_region = FALSE,
      interpolate_region_range = rep(list(NULL), spectral_preprocessing_n),
      do_bc = TRUE,
      bc_method = "rubberband",
      bc_cutoff = 5, 
      bc_do_impute = TRUE,
      do_smooth = c(FALSE, TRUE, TRUE),
      smooth_method = "sg",
      smooth_p = 3,
      smooth_n = 31,
      smooth_m = c(0, 1, 2),
      smooth_ts = 1,
      smooth_k = NA_real_,
      do_normalise = TRUE,
      normalise_method = "snv",
      do_bin = TRUE,
      bin_width = 10,
      bin_new_x_type = "mean",
      do_scale = FALSE, # ---note: don't scale now because this can only be done with the final data subset
      scale_center = NA,
      scale_scale = NA,
      id_preprocessing = seq_along(scale_scale)
    )  
  
  d_spectral_preprocessing_settings
  
}


#' Preprocess MIRS for model development
#' 
#' @export
dd_make_data_model_preprocessed <- function(dd_data_model, dd_mir_preprocessing_settings) {
  
  
  res <- 
    dd_data_model |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_clip(range = dd_mir_preprocessing_settings$clip_range[[1]]) |>
    ir::ir_interpolate_region(range = tibble::tibble(start = c(2300), end = c(2380))) #---note: to remove influence of CO2 peak on normalization
    
  #res <- dd_data_model
  d_spectral_preprocessing_settings <- dd_mir_preprocessing_settings
  
  d_spectral_preprocessing_settings |>
    dplyr::mutate(
      spectra =
        purrr::map(id_preprocessing, function(i) {
          
          # do the preprocessing with `irpeat::ir_preprocess()`
          res <- 
            res |>
            irpeat::irp_preprocess(
              do_interpolate = TRUE,
              interpolate_start = interpolate_start[[i]][[1]],
              interpolate_dw = interpolate_dw[[i]],
              do_clip = TRUE,
              clip_range = clip_range[[i]],
              do_interpolate_region = FALSE,
              interpolate_region_range = interpolate_region_range[[i]][[1]],
              do_bc = do_bc[[i]],
              bc_method = bc_method[[i]],
              bc_cutoff = bc_cutoff[[i]],
              bc_do_impute = bc_do_impute[[i]], 
              do_smooth = do_smooth[[i]],
              smooth_method = smooth_method[[i]],
              smooth_p = smooth_p[[i]],
              smooth_n = smooth_n[[i]],
              smooth_m = smooth_m[[i]],
              smooth_ts = smooth_ts[[i]],
              smooth_k = smooth_k[[i]],
              do_normalise = do_normalise[[i]],
              normalise_method = normalise_method[[i]],
              do_bin = do_bin[[i]],
              bin_width = bin_width[[i]],
              bin_new_x_type = bin_new_x_type[[i]],
              do_scale = do_scale[[i]],
              scale_center = scale_center[[i]],
              scale_scale = scale_scale[[i]],
              do_return_as_ir = TRUE
            ) |>
            ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000)))
          
        })
    ) |>
    dplyr::pull(spectra)
  
}


#' Creates a list with preprocessing configurations that is compatible with irpeat
#' 
#' @export
dd_make_mir_preprocessing_config <- function(dd_mir_preprocessing_settings, dd_model_info) {
  
  purrr::map(seq_len(nrow(dd_model_info)), function(i) {
    
    res <- 
      dd_mir_preprocessing_settings |>
      dplyr::slice(dd_model_info$id_preprocessing[[i]]) |>
      dplyr::select(-id_preprocessing) |>
      dplyr::mutate(
        do_scale = TRUE,
        do_bc = TRUE
      ) |>
      unclass()
    
    res$scale_center <- dd_model_info$x_center[[i]]
    res$scale_scale <- dd_model_info$x_scale[[i]]
    
    list(
      irp_preprocess = res
    )
    
  })
  
}


#' Preprocesses spectra according to a preprocessing config object
#' 
#' @param x An ir object to be preprocessed.
#' 
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#' 
#' @export
irp_preprocess_eb1149 <- function (x, config) {
  
  res <- x
  
  irpeat::irp_preprocess(
    x = res, 
    do_interpolate = config$irp_preprocess$do_interpolate, 
    interpolate_start = config$irp_preprocess$interpolate_start[[1]], 
    interpolate_dw = config$irp_preprocess$interpolate_dw, 
    do_clip = config$irp_preprocess$do_clip[[1]], 
    clip_range = config$irp_preprocess$clip_range[[1]], 
    do_interpolate_region = FALSE, 
    interpolate_region_range = config$irp_preprocess$interpolate_region_range[[1]], 
    do_bc = config$irp_preprocess$do_bc, 
    bc_method = config$irp_preprocess$bc_method, 
    bc_cutoff = config$irp_preprocess$bc_cutoff, 
    bc_do_impute = config$irp_preprocess$bc_do_impute, 
    do_smooth = config$irp_preprocess$do_smooth, 
    smooth_method = config$irp_preprocess$smooth_method, 
    smooth_p = config$irp_preprocess$smooth_p, 
    smooth_n = config$irp_preprocess$smooth_n, 
    smooth_m = config$irp_preprocess$smooth_m, 
    smooth_ts = config$irp_preprocess$smooth_ts, 
    smooth_k = config$irp_preprocess$smooth_k, 
    do_normalise = config$irp_preprocess$do_normalise, 
    normalise_method = config$irp_preprocess$normalise_method, 
    do_bin = config$irp_preprocess$do_bin, 
    bin_width = config$irp_preprocess$bin_width, 
    bin_new_x_type = config$irp_preprocess$bin_new_x_type, 
    do_scale = FALSE, 
    scale_center = config$irp_preprocess$scale_center, 
    scale_scale = config$irp_preprocess$scale_scale, 
    do_return_as_ir = TRUE
  ) |>
    ir::ir_clip(range = data.frame(start = c(650, 2400), end = c(2250, 4000))) |>
    ir::ir_scale(center = config$irp_preprocess$scale_center, scale = config$irp_preprocess$scale_scale)
  
}


#' Makes predictions with a model
#' 
#' @param x An ir object
#' 
#' @param config A list with arguments compatible with `irpeat::irp_preprocess()`.
#' 
#' @export
irp_degree_of_decomposition <- function(x, model, config, prediction_domain) {
  
  x_or <- x
  x <- irp_preprocess_eb1149(x = x, config = config)
  x_in_pd <- irpeat::irp_is_in_prediction_domain(x = x, prediction_domain = prediction_domain)
  
  newdata <- 
    tibble::tibble(
      x = 
        x |>
        ir::ir_flatten() |>
        dplyr::select(-1) |>
        t()
    )
  
  res <- 
    tidybayes::predicted_rvars(
      newdata = newdata,
      object = model,
      value = "yhat"
    )
    
  x_or$degree_of_decomposition <- res$yhat
  x_or$degree_of_decomposition_in_pd <- x_in_pd$is_in_prediction_domain
  
  x_or
  
}




#' Preprocessing pipeline for spectra for plotting
#' 
#' @param x An object of class `ir`.
#' 
#' @export
dd_preprocess_mir_for_plotting <- function(x) {
  
  x |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_clip(range = data.frame(start = 650, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
    
}


