#' Creates a connection to the pmird database
#' 
#' @export
dd_get_pmird_connection <- function() {
  
  RMariaDB::dbConnect(
    drv = RMariaDB::MariaDB(),
    dbname = "pmird",
    username = "root", # ---todo: adjust or get from config file
    password = "coucou",
    host = "mariadb"
  )
  
}

#' Pipe-friendly `do.call()`
#' 
#' @param f Character value. Name of function to apply
#' 
#' @export
dd_do_call <- function(x, f) {
  do.call(f, x)
}


#' Saving and loading tibbles with rvars columns
#'
#' See https://github.com/stan-dev/posterior/issues/307
#'
#' @param object Object to save.
#'
#' @param file Character value. Path to the file.
#'
#' @export
saveRDS_rvars <- function(object, file) {
  
  saveRDS(
    object = object,
    file = file,
    refhook = \(x) if (any(c("vec_proxy", "vec_proxy_equal") %in% names(x))) ""
  )
  
}

#' @export
readRDS_rvars <- function(file) {
  
  readRDS(
    file = file,
    refhook = \(x) new.env()
  )
  
}


#' Gets a model name for a model id
#' 
#' @export
dd_model_names <- function(x, dd_model_info) {
  dd_model_info |>
    dplyr::mutate(
      model_name =
        dplyr::case_when(
          id_model == 1 ~ "1",
          id_model == 2 ~ "2",
          id_model == 3 ~ "3"
        )
    ) |>
    dplyr::slice(match(!!x, id_model)) |>
    dplyr::pull(model_name)
}


#' RMSE
#' 
#' @export
dd_rmse <- function(yhat, y) {
  
  sqrt(mean((yhat - y)^2))
  
}

#' @describeIn dd_rmse
dd_rmse_rvar <- function(yhat, y) {
  
  sqrt(posterior::rvar_mean((yhat - y)^2)) 
  
}



#### Helper functions for the paper ####

#' Formats author list
#' 
#' @export
dd_make_author_list <- function(x) {
  
  purrr::map_chr(x, function(.x) {
    paste0(.x$name, "$^{\\text{", .x$affil, "}}$")
  }) |>
    paste(collapse = "  \n")
  
}

#' Formats affiliation list
#' 
#' @export
dd_make_affiliation_list <- function(x) {
  
  purrr::map_chr(x, function(.x) {
    paste0("$^{", .x$number, "}$ ", .x$text)
  }) |>
    paste(collapse = "  \n")
  
}


#' Extracts required bibtex keys for citations of data in the pmird database
#' 
#' @export
dd_make_citations_data <- function(dd_data_model, dd_data_pmird_peat_cores) {
  
  # data citations
  con <- dd_get_pmird_connection()
  pmird_undecomposed_litter_citations <- 
    dd_data_model |> 
    dplyr::filter(dd_dataset_label == "pmird_undecomposed_litter") |> 
    dplyr::pull(id_measurement) |>
    pmird::pm_get_citations(con = con)
  dd_data_pmird_peat_cores_citations <- 
    dd_data_pmird_peat_cores |> 
    dplyr::pull(id_measurement) |>
    pmird::pm_get_citations(con = con)
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  list(
    pmird_undecomposed_litter_citations = pmird_undecomposed_litter_citations,
    dd_data_pmird_peat_cores_citations = dd_data_pmird_peat_cores_citations
  )
  
}


#' Formats `dd_stan_1_kfold_comparison` and alike for creating a table to print in the paper
#' 
#' @param kfold_comparison An object with the same structure as `dd_stan_1_kfold_comparison`.
#' 
#' @param kfold_target_name Name of the target that created the kfold-CV models 
#' used to compute the comparison in `kfold_comparison`.
#' 
#' @export
dd_make_kfold_comparison_dataframe <- function(kfold_comparison, kfold_target_name = "dd_stan_1_kfold", dd_model_info) {
  
  kfold_comparison |>
    as.data.frame() |>
    dplyr::mutate(
      model_name = 
        match(targets::tar_branch_names_raw(kfold_target_name, 1:3), rownames(kfold_comparison)) |>
        dd_model_names(dd_model_info = dd_model_info)
    ) |>
    dplyr::select(model_name, elpd_diff, se_diff)
  
}



#### Output ####

#' Prepares output for the irpeatmodels package
#' 
#' @export
dd_make_output_1 <- function(dd_model_info, dd_stan_1_fit, dd_mir_preprocessing_config) {
  
  stopifnot(nrow(dd_model_info) == length(dd_stan_1_fit) && nrow(dd_model_info) == length(dd_mir_preprocessing_config))
  
  file_names <- character(length = nrow(dd_model_info))
  for(i in seq_len(nrow(dd_model_info))) {
    res <- 
      list(
        brms_model = dd_stan_1_fit[[i]],
        config = dd_mir_preprocessing_config[[i]],
        prediction_domain = dd_model_info$prediction_domain[[i]]
      )
    file_names[[i]] <- paste0("output_for_irpeatmodels/degree_of_decomposition_", dd_model_info$id_model[[i]], ".rds")
    saveRDS(res, file_names[[i]])
  }
  
  file_names
  
}

