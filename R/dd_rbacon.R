#' Creates csv files from `dd_data_pmird_peat_cores_14C` that can be used with rbacon
#' 
#' @export
dd_prepare_rbacon_csv <- function(dd_data_pmird_peat_cores_14C) {
  
  # crate folders if needed
  root_dir <- "data/derived_data/Bacon_runs/"
  core_label <- unique(dd_data_pmird_peat_cores_14C$core_label)
  target_dir <- paste0(root_dir, core_label) 
  for(.x in target_dir) {
    if(! dir.exists(.x)) {
      dir.create(.x)
    }
  }
  
  # create csv files
  res_csv <- 
    dd_data_pmird_peat_cores_14C |>
    dplyr::mutate(
      labID = lab_code_14C,
      age = age_14C,
      error = age_14C_err,
      depth = (sample_depth_upper + sample_depth_lower) / 2,
      cc = 1
    ) |>
    dplyr::arrange(core_label, sample_depth_upper) |>
    dplyr::select(dplyr::all_of(c("core_label", "sampling_date", "labID", "age", "error", "depth", "cc"))) |>
    dplyr::group_split(core_label) |>
    purrr::map(function(.x) {
      res <- 
        dplyr::bind_rows(
        .x |> 
          dplyr::slice(1) |>
          dplyr::mutate(
            age = 1950 - lubridate::year(sampling_date),
            error = 1,
            depth = 0,
            cc = 0
          ),
        .x
      )
      
      add_date_1 <- IntCal::pMC.age(mn = 131.78, sdev = 0.39, ratio = 100, decimals = 0)
      
      # manually add dates (modern ages)
      if(unique(na.omit(res$core_label)) == "eb1018_OD2" && ! "Poz-109925" %in% res$labID) {
        res <- 
          dplyr::bind_rows(
            res,
            tibble::tibble(
              labID = "Poz-109925",
              age = add_date_1[[1]],
              error = add_date_1[[2]],
              depth = 15.5,
              cc = 1
            )
          )
      } 
      
      res |>
        dplyr::select(-core_label, -sampling_date) |>
        dplyr::arrange(depth)
        
    })
  
  # write csv
  res_file <- paste0(target_dir, "/", core_label, ".csv")
  for(i in seq_along(res_csv)) {
    write.csv(res_csv[[i]], file = res_file[[i]], row.names = FALSE)
  }
  
  names(res_file) <- core_label
  res_file
  
}



#' Runs rbacon for cores
#' 
#' @param hiatus.depths,postbomb See `rbacon::Bacon()`.
#' 
#' @export
dd_run_rbacon <- function(dd_rbacon_csv, hiatus.depths, postbomb, dd_stan_1_mcmc_settings) {
  
  core_label <- names(dd_rbacon_csv)
  
  rbacon::Bacon(
    core = core_label, 
    coredir = "data/derived_data/Bacon_runs/", 
    run = TRUE,
    thick = 5,
    prob = 0.95,
    add.bottom = TRUE,
    d.by = 1,
    depth.unit = "cm",
    age.unit = "yr",
    acc.shape = 1.5,
    acc.mean = 
      switch(
        core_label,
        "eb1005_MH1" =,
        "eb1018_OD2" = 10,
        "eb1006_MK1" = 25,
      ),
    mem.strength = 10,
    mem.mean = 0.5,
    boundary = NA,
    hiatus.depths = hiatus.depths,
    hiatus.max = 10000,
    accept.suggestions = FALSE,
    ask = FALSE,
    ssize = with(dd_stan_1_mcmc_settings, (iter - warmup) * chains),
    save.info = TRUE,
    postbomb = postbomb,
    age.min = 
      read.csv(paste0("data/derived_data/Bacon_runs/", core_label, "/", core_label, ".csv")) |>
      dplyr::filter(depth == 0) |>
      dplyr::pull(age)
  )
  
  core_label
  
}



#' Extracts depths from the rbacon age-depth model of a core for specified ages
#' 
#' Wrapper around `Bacon::Bacon.Age.d()` that takes care of loading the correct
#' age-depth model.
#' 
#' @param x Numeric vector. The depth values for which to extract ages.
#' 
#' @param bacon_is_loaded Logical value. Can be set to `TRUE` if the correct 
#' age-depth model is already loaded, to avoid unnecessary computations.
#' 
#' @param ... Arguments passed to `rbacon::Bacon.Age.d()`.
#' 
#' @export
dd_get_age_for_depth <- function(x, core_label, postbomb, ssize, bacon_is_loaded = FALSE, ...) {
  
  # load core from Bacon folder
  if(! bacon_is_loaded) {
    rbacon::Bacon(
      core = core_label, 
      coredir = "data/derived_data/Bacon_runs/", 
      run = FALSE,
      accept.suggestions = TRUE,
      ask = FALSE,
      thick = 5,
      prob = 0.95,
      add.bottom = TRUE,
      d.by = 1,
      depth.unit = "cm",
      age.unit = "yr",
      acc.shape = 1.5,
      acc.mean = 10,
      mem.strength = 10,
      mem.mean = 0.5,
      boundary = NA,
      postbomb = postbomb,
      ssize = ssize
    )
    rbacon::agedepth()
    dev.off()
  }
  
  # get ages
  res <- 
    purrr::map(x, function(.x) {
      rbacon::Bacon.Age.d(.x, ...)
    }) |>
    dd_do_call("cbind") |>
    posterior::rvar()
  
  res
  
}




