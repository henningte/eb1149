#### Functions to load and harmonize data sources ####

#' Extracts data on undecomposed litter from the pmird database
#' 
#' @export
dd_make_data_pmird_undecomposed_litter <- function() {
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dplyr::filter(sample_type %in% c("vegetation", "litter") & ! id_dataset %in% c(1, 24) & is.na(sample_depth_lower) & ! is.na(taxon_rank_value)) |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dplyr::filter(! is.na(mirs_file)) |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res <-
    res |>
    pmird::pm_load_spectra(directory = "data/raw_data/pmird-database") |>
    dplyr::mutate(
      mass_relative_mass = 1.0, #---note: assumed
      dd_dataset_label = "pmird_undecomposed_litter"
    )
  
  res
  
}


#' Extracts data from Reuter.2019 from the pmird database and attaches additional data on mass losses from Pangeae
#' 
#' @export
dd_make_data_reuter2019 <- function() {
  
  requireNamespace("magrittr", quietly = TRUE)
  
  # modified from pmird script:
  d2 <-
    readRDS("data/raw_data/pangaea_metadata_reuter2019.rds")[[2]]$data %>%
    dplyr::slice(-c(1:3)) %>% # smaples without MIRS
    dplyr::select(2:11) %>%
    stats::setNames(c("sample_label", "sampling_latitude", "sampling_longitude", "taxon_organ", "sampling_date", "sampling_site", "decomposition_site", "C", "N", "mass_loss")) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(c("C", "N", "mass_loss")),
        function(x) x/100
      ),
      mass_relative_mass = 
        ifelse(is.na(mass_loss), 1.0, 1.0 - mass_loss),
      decomposition_site =
        dplyr::case_when(
          decomposition_site == "" ~ NA_character_,
          TRUE ~ decomposition_site
        ),
      sampling_date =
        dplyr::case_when(
          sampling_date == "" ~ NA_character_,
          sampling_date != "" ~ paste0(sampling_date, "-01")
        ) %>%
        as.Date(),
      taxon_rank_value = "Phragmites australis", # from Reuter.2020
      taxon_rank_name = "species",
      taxon_organ = 
        dplyr::case_when(
          stringr::str_detect(taxon_organ, "Leaf") ~ "leaves",
          stringr::str_detect(taxon_organ, "Rhizom") ~ "rhizomes"
        ),
      sample_type = "vegetation",
      id_replicate =
        rep_len(1:3, length.out = nrow(.)), # ---note: assuming that also for transplanted samples, the order is the 1, 2, 3 and these replicate ids match.
      sample_treatment =
        dplyr::case_when(
          !is.na(sampling_date) ~ "control",
          TRUE ~ "treatment"
        ),
      comments_samples =
        dplyr::case_when(
          !is.na(sampling_date) ~ "Only year and month of sample collection known. sample_label has the following coding: (1) site where the sample was collected, (2) site substrate in which the sample was incubated during the decomposition experiment, (3) taxon_organ, (4) replicate id. Samples both before incubatin and after incubation are included in the database.",
          TRUE ~ "These are the samples at the end of the decomposition experiment. sample_label has the following coding: (1) site where the sample was collected, (2) site substrate in which the sample was incubated during the decomposition experiment, (3) taxon_organ, (4) replicate id. Samples both before incubatin and after incubation are included in the database."
        ),
      id_sample = seq_len(nrow(.))
    ) %>%
    dplyr::mutate(
      sample_label =
        paste0(
          sampling_site %>% stringr::str_sub(start = 1L, end = 1L), "_",
          decomposition_site %>% stringr::str_sub(start = 1L, end = 1L), "_",
          taxon_organ, "_",
          id_replicate
        )
    ) |>
    dplyr::select(sample_label, mass_relative_mass)
  
  
  ## spectra from pmird
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dplyr::filter(id_dataset == 24) |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dplyr::filter(! is.na(mirs_file)) |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res <-
    res |>
    pmird::pm_load_spectra(directory = "data/raw_data/pmird-database") |>
    dplyr::left_join(
      d2,
      by = "sample_label"
    ) |>
    dplyr::mutate(
      dd_dataset_label = "reuter2019"
    )
  
  res
  
}


#' Prepares data from Arsenault.2024a (eb1064)
#' 
#' @export
dd_make_data_Arsenault2024a <- function() {
  
  dir_source <- "data/raw_data/eb1064/"

  # import
  d_samples <- readRDS(paste0(dir_source, "samples.rds"))
  
  d_mir <- 
    readRDS(paste0(dir_source, "d_mir.rds")) |>
    dplyr::left_join(
      d_samples |>
        dplyr::select(dplyr::starts_with("id_sample")),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      spectra =
        purrr::map(spectra, function(.x) {
          if(is.null(.x)) {
            tibble::tibble(x = spectra[[1]]$x, y = NA_real_)
          } else {
            .x
          }
        })
    ) |>
    ir::ir_as_ir()
  
  d_data <- 
    readRDS(paste0(dir_source, "data.rds")) |>
    dplyr::left_join(
      d_samples |>
        dplyr::select(dplyr::starts_with("id_sample")),
      by = "id_sample"
    )
  
  # wrangle
  res <- 
    d_mir |> 
    dplyr::left_join(
      d_data |> 
        dplyr::filter(attribute_name == "mass_relative_mass") |> 
        dplyr::select(id_sample, value) |> 
        dplyr::rename(mass_relative_mass = "value"), 
      by = "id_sample"
    ) |>
    dplyr::left_join(
      d_data |> 
        dplyr::filter(attribute_name == "C_relative_mass") |> 
        dplyr::select(id_sample, value, error) |> 
        dplyr::rename(
          C = "value",
          C_err = "error"
        ), 
      by = "id_sample"
    ) |>
    dplyr::left_join(
      d_data |> 
        dplyr::filter(attribute_name == "N_relative_mass") |> 
        dplyr::select(id_sample, value, error) |> 
        dplyr::rename(
          N = "value",
          N_err = "error"
        ),
      by = "id_sample"
    ) |>
    dplyr::left_join(
      d_data |> 
        dplyr::filter(attribute_name == "13C") |> 
        dplyr::select(id_sample, value, error) |> 
        dplyr::rename(
          d13C = "value",
          d13C_err = "error"
        ), 
      by = "id_sample"
    ) |>
    dplyr::left_join(
      d_data |> 
        dplyr::filter(attribute_name == "15N") |> 
        dplyr::select(id_sample, value, error) |> 
        dplyr::rename(
          d15N = "value",
          d15N_err = "error"
        ), 
      by = "id_sample"
    ) |>
    dplyr::left_join(
      d_samples, 
      by = c("id_sample", "id_dataset", "id_sample_origin", "id_sample_incubation_start", "id_sample_parent")
    ) |>
    dplyr::select(! dplyr::all_of(c("attribute_name", "value_type"))) |>
    dplyr::filter(sample_type != "sediment") |>
    dplyr::filter(incubation_duration > 0 | (incubation_duration == 0 & ! duplicated(taxon_rank_value))) |>
    dplyr::filter(! is.na(mass_relative_mass) & ! purrr::map_lgl(spectra, function(.x) is.null(.x) || all(is.na(.x$y))))
  
  res <- 
    res |>
    dplyr::mutate(
      dd_dataset_label = "arsenault2024a"
    )
  
  res
  
}


#' Wrapper function that loads and combines data used for model development
#' 
#' @export
dd_make_data_model <- function() {
  
  dplyr::bind_rows(
    dd_make_data_pmird_undecomposed_litter(),
    dd_make_data_reuter2019(),
    dd_make_data_Arsenault2024a()
  ) |>
    dplyr::mutate(
      id_sample_original = id_sample,
      id_sample = seq_along(id_sample_original)
    )
  
}




#' Extracts peat depth profiles from the pmird database (except macrofossils)
#' 
#' @export
dd_make_data_pmird_peat_cores <- function() {
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dplyr::filter((stringr::str_detect(core_label, "(eb1018|eb1006|eb1005)") | id_dataset == 11  | id_dataset == 8) & sample_type == "peat") |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dplyr::filter(! is.na(mirs_file)) |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res <- 
    res |>
    dplyr::filter(! (id_dataset == 8 & core_label == 5)) |>
    dplyr::mutate(
      core_label =
        dplyr::case_when(
          id_dataset == 8 & core_label == "2" ~ "1",
          TRUE ~ core_label
        )
    )
  
  res <-
    res |>
    pmird::pm_load_spectra(directory = "data/raw_data/pmird-database") |>
    dplyr::mutate(
      dd_dataset_label = "pmird_peat_cores"
    )
  
  res
  
}


#' Extracts peat depth profiles from the pmird database (only macrofossils)
#' 
#' @export
dd_make_data_pmird_peat_cores_macrofossils <- function() {
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <- 
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dplyr::select(-taxon_rank_value, -taxon_rank_name, -taxon_organ) |>
    dplyr::filter((stringr::str_detect(core_label, "(eb1018|eb1006|eb1005)") | id_dataset == 11  | id_dataset == 8) & sample_type == "peat") |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(macrofossils, by = "id_measurement") |>
    dm::left_join(macrofossil_types, by = "id_macrofossil_type") |>
    dm::left_join(mactyps_to_taxclas, by = "id_macrofossil_type") |>
    dm::left_join(taxonomic_classifications, by = "id_taxonomic_classification") |>
    dplyr::filter(! is.na(macrofossil_volume_fraction)) |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res
  
}



#' Extracts a peat sample with high mineral content from the pmird database
#' 
#' @export
dd_make_data_pmird_mineral <- function() {
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dplyr::filter(! (is.na(C) | is.na(loss_on_ignition)) & ! is.na(mirs_file)) |>
    dplyr::filter(C < 0.1 | loss_on_ignition < 0.1) |>
    dplyr::filter(sample_label == "PM8.19") |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res <-
    res |>
    pmird::pm_load_spectra(directory = "data/raw_data/pmird-database") |>
    dplyr::mutate(
      dd_dataset_label = "pmird_mineral"
    )
  
  res
  
}


#' Computes or predicts decomposition indicators for an object of class `ir`
#' 
#' @param x An object of class `ir`.
#' 
#' @export
dd_make_data_pmird_add_irpeat <- function(x) {
  
  res <- 
    x |>
    irpeat::irp_bulk_density_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_C_to_N_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_O_to_C_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_H_to_C_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_dgf0_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_nosc_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_d13C_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_d15N_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_macroporosity_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_eac_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_edc_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_nitrogen_content_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_carbon_content_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    dplyr::mutate(
      dgf0_1 = units::set_units(dgf0_1, kJ/mol)
    )
  
  res_hi <- 
    x |>
    #irpeat::irp_carbon_content_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    ir::ir_interpolate(start = NULL, dw = 1) |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area") |>
    irpeat::irp_hi(x1 = 1630, x2 = 1090) |>
    dplyr::select(dplyr::starts_with("hi_") & dplyr::ends_with("1090"))
  
  dplyr::bind_cols(
    res, 
    res_hi
  )
  
}



#' Prepares peat accumulation data from Galka
#' 
#' @export
dd_make_data_pmird_peat_cores_accumulation <- function() {
  
  res <- 
    read.csv("private/MhII CAR. full data Henning.csv") |>
    dplyr::select(4, 5, 15, 17, 20) |>
    setNames(nm = c("sample_depth_upper", "sample_depth_lower", "age", "loss_on_ignition", "bulk_density_ogranic_matter")) |>
    dplyr::filter(! is.na(sample_depth_upper)) |>
    dplyr::mutate(
      dplyr::across(dplyr::everything(), as.numeric),
      core_label = "eb1005_MH1",
      dd_dataset_label = "pmird_peat_cores"
    ) |>
    dplyr::relocate(dplyr::all_of(c("dd_dataset_label", "core_label")), .before = dplyr::everything())
  
}


#' Age-depth models for cores
#' 
#' @export
dd_make_data_age_depth_model <- function() {
 
  files <- c("data/raw_data/d88_1.rds", "data/raw_data/d89_1.rds")
  
  purrr::map(files, function(.x) {
    res <- readRDS(.x)
    tibble::tibble(
      core_label = res$core_label[[1]],
      model = list(approxfun(x = res$sample_depth_mid, y = res$age_bp, method = "linear", yleft = NA_real_, yright = NA_real_))
    )
  }) |>
    dplyr::bind_rows()

}

#' TA-WTD models for cores
#' 
#' @export
dd_make_data_ta_wtd <- function() {
  
  files <- c("data/raw_data/d88_2.rds", "data/raw_data/d89_2.rds")
  
  purrr::map(files, function(.x) {
    res <- readRDS(.x)
    tibble::tibble(
      core_label = res$core_label[[1]],
      model = list(approxfun(x = res$sample_depth_mid, y = res$water_table_depth_testate_amoebae, method = "linear", yleft = NA_real_, yright = NA_real_))
    )
  }) |>
    dplyr::bind_rows()
  
}


#' Raw radiocarbon data for the peat cores
#' 
#' @export
dd_make_data_pmird_peat_cores_14C <- function() {
  
  # connect to database
  con <- dd_get_pmird_connection()
  
  dm_pmird <-
    pmird::pm_get_dm(con, learn_keys = TRUE)
  
  res <-
    dm_pmird |>
    dm::dm_zoom_to(samples) |>
    dplyr::filter((stringr::str_detect(core_label, "(eb1018|eb1006|eb1005)") | id_dataset == 11  | id_dataset == 8)) |>
    dm::left_join(data_to_samples, by = "id_sample") |>
    dm::left_join(data, by = "id_measurement") |>
    dm::left_join(mir_metadata, by = "id_measurement") |>
    dplyr::filter(! is.na(age_14C)) |>
    dm::pull_tbl() |>
    tibble::as_tibble()
  
  on.exit(RMariaDB::dbDisconnect(con))
  
  res <- 
    res |>
    dplyr::filter(! (id_dataset == 8 & core_label == 5)) |>
    dplyr::mutate(
      core_label =
        dplyr::case_when(
          id_dataset == 8 & core_label == "2" ~ "1",
          TRUE ~ core_label
        )
    )
  
  res <-
    res |>
    dplyr::mutate(
      dd_dataset_label = "pmird_peat_cores_radiocarbon"
    )
  
  res
  
}


#' Imports file `ecy2462-sup-0005-npp_moss.csv` from "A peatland productivity and decomposition parameter database"
#' 
#' @export
dd_make_data_bona2018_moss_npp <- function() {
  
  res <- 
    read.csv("data/raw_data/ecy2462-sup-0005-npp_moss.csv")
  colnames(res) <- tolower(colnames(res))
  res
  
}


#' Imports file `Bengtsson_etal_sph_holarctic_growth.csv` from Bengtsson et al. (2021)
#' 
#' @export
dd_make_data_bengtsson2021_moss_npp <- function() {
  
  res <- 
    read.csv("data/raw_data/Bengtsson_etal_sph_holarctic_growth.csv")
  colnames(res) <- tolower(colnames(res))
  res
  
}


#### Simulations ####

#' Simulates the bias one makes when using the MIRS predicted degree of decomposition for a two-component system
#' 
#' @export
dd_make_simulation_2 <- function() {
  
    tidyr::expand_grid(
      mf_1_1 = seq(0, 1.0, length.out = 500),
      m_1_1 = 1.0, 
        #10^seq(log10(1/10000000), log10(10000000), length.out = 50),
      #m_2_1 = c(0, 1), 
        #seq(0, 1.0, length.out = 20),
      gamma_1 = c(seq(0.0, 0.9, by = 0.1), 0.999),
      gamma_2 = c(seq(0.0, 0.9, by = 0.1), 0.999)
    ) |>
    dplyr::filter(! (gamma_1 == 1 & gamma_2 == 1)) |>
    dplyr::mutate(
      m_1_1 =
        dplyr::case_when(
          mf_1_1 == 0 ~ 0,
          TRUE ~ m_1_1,  
        ),
      m_2_1 = 
        dplyr::case_when(
          mf_1_1 == 0 ~ 1, #---note: dummy value
          TRUE ~ m_1_1 / mf_1_1 - m_1_1,  
        ),
      m_1_0 = m_1_1 / (1 - gamma_1),
      m_2_0 = m_2_1 / (1 - gamma_2),
      m_0 = m_1_0 + m_2_0,
      m_1 = m_1_1 + m_2_1,
      gamma_mirs = 
        dplyr::case_when(
          m_1_1 == 0 ~ gamma_2,
          m_2_1 == 0 ~ gamma_1,
          TRUE ~ 1/m_1 * (m_1_1 * gamma_1 + m_2_1 * gamma_2)
        ),
      gamma = 
        dplyr::case_when(
          m_1_0 == 0 ~ gamma_2,
          m_2_0 == 0 ~ gamma_1,
          TRUE ~ 1/m_0 * (m_1_0 * gamma_1 + m_2_0 * gamma_2)
        ),
      bias = gamma_mirs - gamma,
      bias_alternative_formula = ((m_1_1 * gamma_1 + m_2_1 * gamma_2) * m_0 - m_1 * m_0 + m_1^2) / (m_0 * m_1)
    )
  
}


#' Mixes spectra of available litter and estimates the bias in the MIRS-predicted degree of decomposition
#' 
#' @export
dd_make_simulation_3 <- function(dd_data_model) {
  
  dd_preprocess_for_bias_1 <- function(x) {
    x |>
      ir::ir_interpolate() |>
      ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
      ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
      ir::ir_normalise(method = "area")
  } 
  
  dd_data_model <- 
    dd_data_model |>
    dd_preprocess_for_bias_1()
  
  # compare fitted values and predictions
  dd_data_model_check_1 <- 
    dd_data_model |>
    irpeat::irp_degree_of_decomposition_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_2(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_3(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd)
  
  # get all undecomposed samples
  d1 <- 
    dd_data_model |>
    dplyr::filter(mass_relative_mass >= 1.0) |>
    dplyr::filter(! duplicated(taxon_rank_value, fromLast = TRUE)) |>
    dplyr::slice_sample(n = 15, replace = FALSE) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::where(is.integer), as.integer
      )
    ) |>
    dplyr::bind_rows(
      dd_data_model |>
        dplyr::filter(taxon_rank_value %in% (dd_data_model |> dplyr::filter(mass_relative_mass < 1.0) |> dplyr::pull(taxon_rank_value) |> unique())) |>
        dplyr::group_split(taxon_rank_value) |>
        purrr::map(function(.x) {
          .x |> 
            dplyr::arrange(mass_relative_mass) |>
            dplyr::slice(1:2) |>
            dplyr::filter(! duplicated(mass_relative_mass)) |>
            dplyr::mutate(
              dplyr::across(
                dplyr::where(is.integer), as.integer
              )
            )
        }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(
          dplyr::across(
            dplyr::where(is.integer), as.integer
          )
        )
    ) |>
    dplyr::mutate(
      scale_factor = list(10^seq(log10(1/1000), log10(1000), length.out = 15)) #list(unique(c(0.0, 0.1, 0.5, 0.8, seq(0.0, 10.0, length.out = 5))))
    ) |>
    tidyr::unnest("scale_factor") %>%
    magrittr::multiply_by(e2 = .$scale_factor)
  
  
  # get five decomposed samples
  d2 <- 
    dd_data_model |>
    dplyr::filter(! paste0(id_dataset, "_", id_sample) %in% paste0(d1$id_dataset, "_", d1$id_sample)) |>
    dplyr::filter(taxon_rank_value %in% (dd_data_model |> dplyr::filter(mass_relative_mass < 1.0) |> dplyr::pull(taxon_rank_value) |> unique())) |>
    dplyr::group_split(taxon_rank_value) |>
    purrr::map(function(.x) {
      .x |> 
        dplyr::arrange(mass_relative_mass) |>
        dplyr::slice(c(1:2, nrow(.x))) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::where(is.integer), as.integer
          )
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      scale_factor_2 = 1
    ) |>
    rep(nrow(d1))
  
  d1 |>
    dplyr::slice(rep(seq_len(nrow(d1)), each = nrow(d2)/nrow(d1))) %>%
    magrittr::add(d2) |>
    dplyr::mutate(
      scale_factor_1 = scale_factor,
      scale_factor_2 = d2$scale_factor_2,
      gamma_1 = (1 - mass_relative_mass),    
      gamma_2 = (1 - d2$mass_relative_mass),   
      m_1 = scale_factor_1 + scale_factor_2,
      scale_factor_1_0 = scale_factor_1/(1 - gamma_1),
      scale_factor_2_0 = scale_factor_2/(1 - gamma_2),
      m_0 = scale_factor_1_0 + scale_factor_2_0,
      taxon_rank_value2 = d2$taxon_rank_value, 
      taxon_organ2 = d2$taxon_organ,
      id_measurement2 = d2$id_measurement,
      gamma = (m_0 - m_1) / m_0,
      bias_hat =  ((gamma_1 * scale_factor_1 + gamma_2 * scale_factor_2) * m_0 - m_1 * m_0 + m_1^2) / (m_1 * m_0) 
    ) |>
    irpeat::irp_degree_of_decomposition_1(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_2(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd) |>
    irpeat::irp_degree_of_decomposition_3(do_summary = TRUE, summary_function_mean = mean, summary_function_sd = posterior::sd)
  
}


#' Uses the decomposition module of the HPM to simulate the maximum expected difference in the degree of decomposition of a two-component system where all parts decompose under the same environmental conditions
#' 
#' @export
dd_make_simulation_4 <- function() {
  
  # degree of decomposition according to the HPM decomposition module
  dd_make_simulation_4_helper_1 <- function(x, k_0, alpha) {
    
    1 - 1/(1 + (alpha - 1) * k_0 * x)^(1/(alpha - 1))
    
  }
  
  dd_simulation_4 <- 
    tidyr::expand_grid(
      alpha = c(1.0000001, 2, 3, 4),
      k_0_1 = c(0.01, 0.05, 0.1),
      #k_0_2 = 10^seq(log10(0.01), log10(0.8), length.out = 10),
      k_0_2 = c(k_0_1, seq(0.1, 0.8, by = 0.1)),
      decomposition_progress = 
        c(seq(0.2, 1, length.out = 70), seq(1, 40, length.out = 70), seq(70, 1000, length.out = 50)) |>
        unique()
    ) |>
    dplyr::mutate(
      gamma_1 = 
        dd_make_simulation_4_helper_1(
          x = decomposition_progress, 
          k_0 = k_0_1, 
          alpha = alpha
        ),
      gamma_2 = 
        dd_make_simulation_4_helper_1(
          x = decomposition_progress, 
          k_0 = k_0_2, 
          alpha = alpha
        ),
      gamma_difference = gamma_1 - gamma_2
    )
  
  dd_simulation_4
  
  # dd_simulation_4 |>
  #   dplyr::arrange(k_0_1, alpha, k_0_2, decomposition_progress) |>
  #   ggplot(aes(y = decomposition_progress, x = (gamma_1 + gamma_2)/2, group = k_0_2)) +
  #   geom_path(aes(color = k_0_2)) +
  #   facet_grid(k_0_1 ~ round(alpha, 0)) 
  
}




