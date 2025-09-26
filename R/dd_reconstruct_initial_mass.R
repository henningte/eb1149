#' Computes an initial mass accumulation rate for eb1005
#' 
#' @export
dd_reconstruct_initial_mass_accumulation_eb1005 <- function(dd_data_pmird_peat_cores_accumulation, dd_data_pmird_peat_cores_yhat, dd_data_pmird_peat_cores) {
  
  # data wrangling
  res <- 
    dd_data_pmird_peat_cores_yhat |>
    dplyr::left_join(
      dd_data_pmird_peat_cores |>
        dplyr::select(id_sample, core_label, sample_depth_upper, sample_depth_lower),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      index_yhat = 
        purrr::map(seq_along(core_label), function(i) {
          index <- which(dd_data_pmird_peat_cores_accumulation$core_label == core_label[[i]] & dd_data_pmird_peat_cores_accumulation$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_accumulation$sample_depth_lower <= sample_depth_lower[[i]])
          if(length(index) == 0) {
            NA
          } else {
            index
          }
        })
    ) |>
    dplyr::filter(! purrr::map_lgl(index_yhat, function(.x) all(is.na(.x)))) |>
    dplyr::select(-sample_depth_upper, -sample_depth_lower) |>
    tidyr::unnest(index_yhat) |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_accumulation |>
        dplyr::mutate(
          duration = c(NA_real_, age[-1] - age[-length(age)]),
          index_yhat = seq_along(sample_depth_lower)
        ),
      by = c("core_label", "index_yhat")
    ) |>
    dplyr::filter(! is.na(duration))
  
  # reconstruct
  a <- 
    res |>
    dplyr::mutate(
      initial_mass = bulk_density_ogranic_matter / (1.0 - degree_of_decomposition),
      initial_mass_rate = initial_mass / duration * 100 * 100
    )
  
  
  
  res <- 
    dd_data_pmird_peat_cores_irpeat |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    )|>
    dplyr::mutate(
      initial_mass = as.numeric(bulk_density_1) / (1.0 - degree_of_decomposition)
    )
  
  
  res |> 
    dplyr::arrange(sample_depth_upper) |>
    dplyr::select(core_label, sample_depth_upper, initial_mass, degree_of_decomposition, bulk_density_1) |>
    ggplot(aes(y = posterior::quantile2(initial_mass, probs = 0.2), x = sample_depth_upper)) +
    geom_path(aes(color = mean(degree_of_decomposition))) +
    geom_path(aes(y = as.numeric(bulk_density_1)), color = "brown") +
    facet_wrap(~ core_label, scales = "free_y")
  
  res |> dplyr::filter(posterior::quantile2(initial_mass, probs = 0.2) > 0.7) |> dplyr::select(core_label, sample_depth_upper, 195:ncol(res)) |>
    dplyr::mutate(dplyr::across(where(posterior::is_rvar), median)) |> View()
  
}





#' Computes an initial mass accumulation rate for eb1006
#' 
#' @export
dd_reconstruct_initial_mass_accumulation_eb1006 <- function(dd_data_pmird_peat_cores_yhat, dd_data_pmird_peat_cores_irpeat, dd_data_age_depth_model, dd_data_ta_wtd) {
  
  dd_data_pmird_peat_cores_macrofossils <- 
    dd_data_pmird_peat_cores_macrofossils |>
    dplyr::filter(stringr::str_detect(taxon_rank_value, "Sphagnum") | stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum|Ericaceae)") | macrofossil_type %in% c("wood", "amorphous organic matter")) |>
    dplyr::select(core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value, macrofossil_taxon_organ, macrofossil_type, macrofossil_volume_fraction) |>
    dplyr::mutate(
      taxon_rank_value = 
        dplyr::case_when(
          stringr::str_detect(taxon_rank_value, "Sphagnum") & stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_hummock",
          stringr::str_detect(taxon_rank_value, "Sphagnum") & ! stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_lawn_hollow",
          stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum)") & stringr::str_detect(macrofossil_taxon_organ, "root") ~ "sedge_root",
          stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum)") & ! stringr::str_detect(macrofossil_taxon_organ, "root") ~ "sedge_aboveground",
          stringr::str_detect(taxon_rank_value, "Ericaceae") | macrofossil_type == "wood" ~ "shrub",
          is.na(taxon_rank_value) & macrofossil_type == "amorphous organic matter" ~ "uom"
        )
    ) |>
    dplyr::group_by(core_label, taxon_rank_value, sample_depth_upper, sample_depth_lower) |>
    dplyr::summarise(
      macrofossil_volume_fraction = sum(macrofossil_volume_fraction, na.rm = TRUE) / 100,
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "taxon_rank_value",
      values_from = "macrofossil_volume_fraction"
    )
  
  # data wrangling
  dd_data_eb1006_age_depth_model <- 
    dd_data_age_depth_model$model[dd_data_age_depth_model$core_label == "eb1006_MK1"][[1]]
  dd_data_eb1018_age_depth_model <- 
    dd_data_age_depth_model$model[dd_data_age_depth_model$core_label == "eb1018_OD2"][[1]]
  dd_data_eb1006_ta_wtd <- 
    dd_data_ta_wtd$model[dd_data_ta_wtd$core_label == "eb1006_MK1"][[1]]
  dd_data_eb1018_ta_wtd <- 
    dd_data_ta_wtd$model[dd_data_ta_wtd$core_label == "eb1018_OD2"][[1]]
  
  res <- 
    dd_data_pmird_peat_cores_yhat |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat |>
        dplyr::select(id_sample, core_label, sample_depth_upper, sample_depth_lower, bulk_density_1, bulk_density_1_in_pd, Ti),
      by = "id_sample"
    ) |>
    dplyr::filter(core_label %in% c("eb1006_MK1", "eb1018_OD2")) |>
    dplyr::mutate(
      age_upper = 
        dplyr::case_when(
          core_label == "eb1006_MK1" ~ dd_data_eb1006_age_depth_model(v = sample_depth_upper),
          core_label == "eb1018_OD2" ~ dd_data_eb1018_age_depth_model(v = sample_depth_upper)
        ),
      age_lower = 
        dplyr::case_when(
          core_label == "eb1006_MK1" ~ dd_data_eb1006_age_depth_model(v = sample_depth_lower),
          core_label == "eb1018_OD2" ~ dd_data_eb1018_age_depth_model(v = sample_depth_lower)
        ),
      duration = age_lower - age_upper,
      wtd_ta = dd_data_eb1006_ta_wtd((sample_depth_upper + sample_depth_lower) / 2)
    )
  
  # reconstruct
  res$degree_of_decomposition_corrected <- res$degree_of_decomposition
  for(i in seq_len(nrow(res))) {
    b <- res$degree_of_decomposition_corrected[[i]]
    res$degree_of_decomposition_corrected[[i]][b > 0.95] <- 0.95
  }
  a <- 
    res |>
    dplyr::mutate(
      initial_mass = as.numeric(bulk_density_1) * 100 * 100 * (sample_depth_lower - sample_depth_upper) / 1000 / (1.0 - degree_of_decomposition_corrected), # kg/m^2
      initial_mass_rate = initial_mass / duration # kg/m^2/yr
    ) |>
    dplyr::mutate(
      macrofossil_volume_fraction_sphagnum_hummock =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$sphagnum_hummock[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$sphagnum_hummock[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
      macrofossil_volume_fraction_sphagnum_lawn_hollow =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$sphagnum_lawn_hollow[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$sphagnum_lawn_hollow[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
      macrofossil_volume_fraction_sedge_aboveground =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$sedge_aboveground[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$sedge_aboveground[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
      macrofossil_volume_fraction_sedge_root =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$sedge_root[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$sedge_root[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
      macrofossil_volume_fraction_sedge =
        dplyr::case_when(
          is.na(macrofossil_volume_fraction_sedge_root) ~ macrofossil_volume_fraction_sedge_aboveground,
          is.na(macrofossil_volume_fraction_sedge_aboveground) ~ macrofossil_volume_fraction_sedge_root,
          TRUE ~ macrofossil_volume_fraction_sedge_root + macrofossil_volume_fraction_sedge_aboveground
        ),
      macrofossil_volume_fraction_shrub =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$shrub[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$shrub[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
      macrofossil_volume_fraction_uom =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$uom[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$uom[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        })
    )
  
  
  a |> 
    dplyr::arrange(sample_depth_upper) |>
    dplyr::select(core_label, sample_depth_upper, sample_depth_lower, age_upper, duration, initial_mass_rate, degree_of_decomposition, bulk_density_1, wtd_ta, macrofossil_volume_fraction_sphagnum_hummock, macrofossil_volume_fraction_sedge_aboveground, macrofossil_volume_fraction_sedge_root, macrofossil_volume_fraction_shrub, macrofossil_volume_fraction_sphagnum_lawn_hollow) |>
    #dplyr::filter(age_upper >= 3500 & age_upper <= 6000) |>
    ggplot(aes(x = age_upper)) +
    ggdist::stat_lineribbon(aes(ydist = initial_mass_rate)) +
    geom_point(aes(y = median(initial_mass_rate), x = age_upper, color = mean(degree_of_decomposition))) +
    scale_fill_brewer() +
    geom_path(aes(y = as.numeric(bulk_density_1) * 100 * 100 * (sample_depth_lower - sample_depth_upper) / 1000 / duration), color = "brown") +
    geom_path(aes(y = wtd_ta / 100), color = "blue") +
    ggnewscale::new_scale_color() +
    geom_path(
      data =
        a |>
        dplyr::select(core_label, age_upper, dplyr::starts_with("macrofossil_volume_fraction_")) |>
        tidyr::pivot_longer(
          ! dplyr::any_of(c("core_label", "age_upper")),
          names_to = "variable",
          values_to = "value"
        ),
      aes(y = value + 1, color = variable)
    ) +
    facet_wrap(~ core_label, scales = "free") +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank()
    )
  
  # restore initial Ti contents
  a |>
    ggplot(aes(y = Ti, x = age_upper)) +
    geom_point() +
    geom_point(aes(y = as.numeric(as.numeric(bulk_density_1) * 100 * 100 * (sample_depth_lower - sample_depth_upper) / 1000 * Ti) / mean(initial_mass)), color = "red") +
    facet_wrap(~ core_label)
  
  a |>
    dplyr::filter(macrofossil_volume_fraction_sphagnum_hummock > 0.9 | macrofossil_volume_fraction_sedge > 0.8) |>
    ggplot(aes(y = mean(degree_of_decomposition_corrected), x = age_upper)) +
    geom_point(aes(color = macrofossil_volume_fraction_sphagnum_hummock > 0.9))
  
  res <- 
    dd_data_pmird_peat_cores_irpeat |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    )|>
    dplyr::mutate(
      initial_mass = as.numeric(bulk_density_1) / (1.0 - degree_of_decomposition)
    )
  
  
  res |> 
    dplyr::arrange(sample_depth_upper) |>
    dplyr::select(core_label, sample_depth_upper, initial_mass, degree_of_decomposition, bulk_density_1) |>
    ggplot(aes(y = posterior::quantile2(initial_mass, probs = 0.2), x = sample_depth_upper)) +
    geom_path(aes(color = mean(degree_of_decomposition))) +
    geom_path(aes(y = as.numeric(bulk_density_1)), color = "brown") +
    facet_wrap(~ core_label, scales = "free_y")
  
  res |> dplyr::filter(posterior::quantile2(initial_mass, probs = 0.2) > 0.7) |> dplyr::select(core_label, sample_depth_upper, 195:ncol(res)) |>
    dplyr::mutate(dplyr::across(where(posterior::is_rvar), median)) |> View()
  
  
  # test: estimate decomposition rate for catotelm peat with clear dominance by S. fuscum
  res_macrofossils <- 
    dd_data_pmird_peat_cores_macrofossils |>
    dplyr::mutate(
      taxon_rank_value =
        dplyr::case_when(
          taxon_rank_value == "Erophorum vaginatum" ~ "Eriophorum vaginatum",
          TRUE ~ taxon_rank_value
        )
    ) |>
    dplyr::filter(macrofossil_volume_fraction > 90 & macrofossil_taxon_organ != "roots" & ! is.na(macrofossil_taxon_organ) & ! taxon_rank_value %in% c("Sphagnum sp."))
  
  res <- 
    dd_data_pmird_peat_cores_yhat |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat,# |>
        #dplyr::select(id_sample, core_label, sampling_date, sample_depth_upper, sample_depth_lower, bulk_density_1, bulk_density_1_in_pd, Ti),
      by = "id_sample"
    ) |>
    dplyr::select(-taxon_rank_value) |>
    dplyr::filter(core_label %in% c("eb1006_MK1", "eb1018_OD2")) |>
    dplyr::left_join(
      res_macrofossils |>
        dplyr::mutate(
          to_keep = TRUE
        ) |>
        dplyr::select(core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value, to_keep) |>
        dplyr::rename(
          sample_depth_upper_macrofossils = "sample_depth_upper",
          sample_depth_lower_macrofossils = "sample_depth_lower"
        ),
      by = dplyr::join_by(core_label, y$sample_depth_upper_macrofossils >= x$sample_depth_upper, y$sample_depth_lower_macrofossils <= x$sample_depth_lower)
    ) |>
    dplyr::filter(to_keep) |>
    dplyr::select(-to_keep)
  
  
  index <- 
    res |>
    dplyr::group_by(core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value) |>
    dplyr::summarise(
      sample_depth_upper_macrofossils = min(sample_depth_upper_macrofossils),
      sample_depth_lower_macrofossils = max(sample_depth_lower_macrofossils),
      .groups = "drop"
    ) |>
    dplyr::filter(sample_depth_upper == sample_depth_upper_macrofossils & sample_depth_lower == sample_depth_lower_macrofossils) |>
    dplyr::mutate(
      to_keep = TRUE
    ) |>
    dplyr::select(! dplyr::ends_with("_macrofossils"))
  
  res <- 
    res |>
    dplyr::left_join(
      index,
      by = c("core_label", "sample_depth_upper", "sample_depth_lower", "taxon_rank_value")
    ) |>
    dplyr::filter(to_keep & ! duplicated(paste0(core_label, "_", sample_depth_upper, "_", sample_depth_lower, "_", taxon_rank_value))) |>
    dplyr::group_split(core_label) |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          age =
            tibble::tibble(
              age_upper = suppressMessages(dd_get_age_for_depth(sample_depth_upper, core_label = unique(core_label), postbomb = 1, ssize = 8000, bacon_is_loaded = FALSE)),
              age_lower = suppressMessages(dd_get_age_for_depth(sample_depth_lower, core_label = unique(core_label), postbomb = 1, ssize = 8000, bacon_is_loaded = TRUE))
            ) |>
            dplyr::mutate(
              dplyr::across(dplyr::everything(), function(.x) (lubridate::year(sampling_date) - 1950) + .x)
            ) |>
            dplyr::mutate(
              wtd_ta = 
                dplyr::case_when(
                  core_label == "eb1006_MK1" ~ dd_data_eb1006_ta_wtd((sample_depth_upper + sample_depth_lower) / 2),
                  core_label == "eb1018_OD2" ~ dd_data_eb1018_ta_wtd((sample_depth_upper + sample_depth_lower) / 2)
                )
            )
        ) |>
        tidyr::unnest("age")
    }) |>
    dd_do_call("rbind")
    
  
  
  res |>
    dplyr::mutate(
      initial_mass_rate = (as.numeric(units::set_units(bulk_density_1, g/(m^2 * cm))) * (sample_depth_lower - sample_depth_upper)) /(1 - degree_of_decomposition) / (age_lower - age_upper)
    ) |>
    ggplot(aes(x = wtd_ta, y = median(initial_mass_rate), color = taxon_rank_value)) +
    ggdist::stat_pointinterval(aes(ydist = initial_mass_rate)) +
    geom_smooth(method = "lm") +
    facet_wrap(~ core_label, scales = "free_x")
    
    
    dplyr::mutate(
      age_upper = 
        dplyr::case_when(
          core_label == "eb1006_MK1" ~ dd_data_eb1006_age_depth_model(v = sample_depth_upper),
          core_label == "eb1018_OD2" ~ dd_data_eb1018_age_depth_model(v = sample_depth_upper)
        ),
      age_lower = 
        dplyr::case_when(
          core_label == "eb1006_MK1" ~ dd_data_eb1006_age_depth_model(v = sample_depth_lower),
          core_label == "eb1018_OD2" ~ dd_data_eb1018_age_depth_model(v = sample_depth_lower)
        ),
      duration = age_lower - age_upper,
      wtd_ta = dd_data_eb1006_ta_wtd((sample_depth_upper + sample_depth_lower) / 2)
    )
  
}



#### Simulations ####

#' Simulates the ratio of the true NPP to the NPP reconstructed naively
#' 
#' @export
dd_make_simulation_1 <- function() {
  
  tidyr::expand_grid(
    m_1_1 = 
      10^seq(log10(1/1000), log10(1000), length.out = 50),
    m_2_1 = 1,
    gamma_1 = c(seq(0.0, 0.90, length.out = 10), 0.999),
    gamma_2 = c(seq(0.0, 0.90, length.out = 10), 0.999)
  ) |>
    #dplyr::filter(! (m_1_0 == 0 & m_2_0 == 0)) |>
    dplyr::mutate(
      m_1_0 = m_1_1 / (1 - gamma_1),
      m_2_0 = m_2_1 / (1 - gamma_2),
      gamma = 
        dplyr::case_when(
          m_1_0 == 0 ~ (1 - m_2_1/m_2_0),
          m_2_0 == 0 ~ (1 - m_1_1/m_1_0),
          TRUE ~ 1/(m_1_1 + m_2_1) * (m_1_1 * (1 - m_1_1/m_1_0) + m_2_1 * (1 - m_2_1/m_2_0))
        ),
      m_0_classical = (m_1_1 + m_2_1) / (1 - gamma),
      m_0_true = m_1_0 + m_2_0,
      ratio_m0 = m_0_classical/m_0_true 
    )
  
}


#### NPP reconstructed for some peat smaples ####


#' Uses the naive approach to reconstruct NPP for some Sphagnum species
#' 
#' @param dd_rbacon Not actually used, just as dummy to allow targets to track 
#' changes because rbacon does not return a model object, but operates via files
#' written to disk.
#' 
#' @export
dd_make_reconstruction_initial_mass_rate_1 <- function(dd_data_pmird_peat_cores_yhat, dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_macrofossils, dd_rbacon, res_file) {
  
  res_macrofossils <- 
    dd_data_pmird_peat_cores_macrofossils |>
    dplyr::mutate(
      taxon_rank_value =
        dplyr::case_when(
          taxon_rank_value == "Erophorum vaginatum" ~ "Eriophorum vaginatum",
          TRUE ~ taxon_rank_value
        ),
      macrofossil_volume_fraction_threshold = 90
    ) |>
    dplyr::filter(macrofossil_volume_fraction > macrofossil_volume_fraction_threshold & macrofossil_taxon_organ != "roots" & ! is.na(macrofossil_taxon_organ) & ! taxon_rank_value %in% c("Sphagnum sp.") & stringr::str_detect(taxon_rank_value, "Sphagnum"))
  
  res <- 
    purrr::map(seq_along(dd_data_pmird_peat_cores_yhat), function(i) {
      dd_data_pmird_peat_cores_yhat[[i]] |>
        dplyr::mutate(
          id_model = i
        )
    }) |>
    dd_do_call("rbind") |>
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat,# |>
      #dplyr::select(id_sample, core_label, sampling_date, sample_depth_upper, sample_depth_lower, bulk_density_1, bulk_density_1_in_pd, Ti),
      by = "id_sample"
    ) |>
    dplyr::select(-taxon_rank_value) |>
    dplyr::filter(core_label %in% c("eb1005_MH1", "eb1006_MK1", "eb1018_OD2")) |>
    dplyr::left_join(
      res_macrofossils |>
        dplyr::mutate(
          to_keep = TRUE
        ) |>
        dplyr::select(core_label, id_sample, sample_depth_upper, sample_depth_lower, taxon_rank_value, macrofossil_volume_fraction_threshold, macrofossil_volume_fraction, to_keep) |>
        dplyr::rename(
          sample_depth_upper_macrofossils = "sample_depth_upper",
          sample_depth_lower_macrofossils = "sample_depth_lower",
          id_sample_macrofossils = "id_sample",
        ),
      by = dplyr::join_by(core_label, y$sample_depth_upper_macrofossils >= x$sample_depth_upper, y$sample_depth_lower_macrofossils <= x$sample_depth_lower)
    ) |>
    dplyr::filter(to_keep) |>
    dplyr::select(-to_keep)
  
  
  index <- 
    res |>
    dplyr::group_by(id_model, core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value) |>
    dplyr::summarise(
      sample_depth_upper_macrofossils = min(sample_depth_upper_macrofossils),
      sample_depth_lower_macrofossils = max(sample_depth_lower_macrofossils),
      id_sample_macrofossils = list(id_sample_macrofossils),
      .groups = "drop"
    ) |>
    dplyr::filter(sample_depth_upper == sample_depth_upper_macrofossils & sample_depth_lower == sample_depth_lower_macrofossils) |>
    dplyr::mutate(
      to_keep = TRUE
    ) |>
    dplyr::select(! dplyr::ends_with("_macrofossils") | dplyr::any_of(c("id_sample_macrofossils")))
  
  res <- 
    res |>
    dplyr::left_join(
      index,
      by = c("id_model", "core_label", "sample_depth_upper", "sample_depth_lower", "taxon_rank_value")
    ) |>
    dplyr::filter(to_keep & ! duplicated(paste0(id_model, "_", core_label, "_", sample_depth_upper, "_", sample_depth_lower, "_", taxon_rank_value))) |>
    dplyr::group_split(core_label) |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          age =
            tibble::tibble(
              age_upper = suppressMessages(dd_get_age_for_depth(sample_depth_upper, core_label = unique(core_label), postbomb = 1, ssize = 8000, bacon_is_loaded = FALSE)),
              age_lower = suppressMessages(dd_get_age_for_depth(sample_depth_lower, core_label = unique(core_label), postbomb = 1, ssize = 8000, bacon_is_loaded = TRUE))
            ) |>
            dplyr::mutate(
              dplyr::across(dplyr::everything(), function(.x) (lubridate::year(sampling_date) - 1950) + .x)
            )
        ) |>
        tidyr::unnest("age")
    }) |>
    dd_do_call("rbind")
  
  ## compute initial mass
  
  # first: Fix long tails in the marginal distribution of the degree of decomposition to avoid Inf initial masses
  res$degree_of_decomposition_plus05 <- res$degree_of_decomposition + 0.2
  for(i in seq_len(nrow(res))) {
    res$degree_of_decomposition[[i]][res$degree_of_decomposition[[i]] > 0.98] <- 0.98
    res$degree_of_decomposition_plus05[[i]][res$degree_of_decomposition_plus05[[i]] > 0.98] <- 0.98
  }
  
  res <- 
    res |>
    dplyr::mutate(
      initial_mass_rate = (as.numeric(units::set_units(bulk_density_1, g/(m^2 * cm))) * macrofossil_volume_fraction/100 * (sample_depth_lower - sample_depth_upper)) /(1 - degree_of_decomposition) / (age_lower - age_upper), #---todo: replace by measured bulk density
      initial_mass_rate_plus05 = (as.numeric(units::set_units(bulk_density_1, g/(m^2 * cm))) * macrofossil_volume_fraction/100 * (sample_depth_lower - sample_depth_upper)) /(1 - degree_of_decomposition_plus05) / (age_lower - age_upper) #---todo: replace by measured bulk density
    )
  
  saveRDS_rvars(res, res_file)
  
  res_file
 
}



#' Summarizes values in `dd_data_bona2018_moss_npp` for species groups
#' 
#' @export
dd_make_data_bona2018_moss_npp_summary_1 <- function(dd_data_bona2018_moss_npp) {
  
  # target_taxa <- 
  #   dd_reconstruction_initial_mass_rate_1 |>
  #   dplyr::pull(taxon_rank_value) |>
  #   unique() |>
  #   stringr::str_split(pattern = "/") |>
  #   unlist() |>
  #   unique()
  
  target_taxa <- c("Sphagnum fuscum", "Sphagnum medium", "Sphagnum divinum", "Sphagnum rubellum", "Sphagnum cuspidatum")
  
  dd_data_bona2018_moss_npp |>
    dplyr::filter(stringr::str_detect(moss_species_or_layer, paste0("(", paste(target_taxa, collapse = "|"), ")"))) |>
    dplyr::rename(
      taxon_rank_value = "moss_species_or_layer",
      npp = "moss_npp"
    ) |>
    dplyr::mutate(
      reference_id =
        dplyr::case_when(
          reference_id == 30 ~ "Moore.1989", 
          reference_id == 43 ~ "Reader.1972",
          reference_id == 47 ~ "Rochefort.1990",
          reference_id == 51 ~ "Szumigalski.1995",
          reference_id == 55 ~ "Thormann.1997",
          TRUE ~ as.character(reference_id)
        )
    ) |>
    dplyr::group_by(taxon_rank_value) |>
    dplyr::summarise(
      npp_mean = median(npp, na.rm = TRUE),
      npp_sd = sd(npp, na.rm = TRUE),
      n = length(npp),
      reference_id = list(unique(reference_id)),
      .groups = "drop"
    ) |>
    dplyr::filter(taxon_rank_value != "Sphagnum magellanicum Brid., Sphagnum fuscum (Schimp.) Klinggr.") |>
    dplyr::mutate(
      taxon_rank_value =
        taxon_rank_value |>
        stringr::str_remove_all(pattern = "(\\(|\\)|\\.|Ehrh|Hedw|Klinggr|Brid|Schimp)") |>
        stringr::str_replace_all(" +, +", ", ") |>
        stringr::str_remove_all(pattern = " +$") |>
        stringr::str_split(", ")
    ) |>
    dplyr::filter(purrr::map_lgl(taxon_rank_value, function(.x) length(.x) == 1)) |>
    tidyr::unnest("taxon_rank_value")
  
} 

#' Summarizes values in `dd_data_bengtsson2021_moss_npp`
#' 
#' @export
dd_make_data_bengtsson2021_moss_npp_summary_1 <- function(dd_data_bengtsson2021_moss_npp) {
  
  dd_data_bengtsson2021_moss_npp |>
    tidyr::pivot_longer(
      dplyr::all_of(c("prod13", "prod14")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::filter(species == "S.fuscum") |>
    dplyr::summarise(
      npp_mean = mean(value, na.rm = TRUE),
      npp_sd = sd(value, na.rm = TRUE)
    )
  
}











