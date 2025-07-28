#' Computes Pearson correlations between decomposition indicators
#' 
#' @export
dd_make_correlations_1 <- function(dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_macrofossils, dd_data_pmird_peat_cores_yhat) {
  
  dd_data_pmird_peat_cores_macrofossils <- 
    dd_data_pmird_peat_cores_macrofossils |>
    dplyr::filter(stringr::str_detect(taxon_rank_value, "Sphagnum") | stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum|Ericaceae)") | macrofossil_type %in% c("wood", "amorphous organic matter")) |>
    dplyr::select(core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value, macrofossil_type, macrofossil_volume_fraction) |>
    dplyr::mutate(
      taxon_rank_value = 
        dplyr::case_when(
          stringr::str_detect(taxon_rank_value, "Sphagnum") & stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_hummock",
          stringr::str_detect(taxon_rank_value, "Sphagnum") & ! stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_lawn_hollow",
          stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum)") ~ "sedge",
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
  
  res <- 
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat,
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
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
      macrofossil_volume_fraction_sphagnum =
        macrofossil_volume_fraction_sphagnum_lawn_hollow + macrofossil_volume_fraction_sphagnum_hummock,
      macrofossil_volume_fraction_sedge =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$sedge[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$sedge[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        }),
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
  

  
  res1 <- 
    res |> 
    dplyr::mutate(degree_of_decomposition = mean(degree_of_decomposition))
  
  # cor(
  #   x = res1$degree_of_decomposition, 
  #   y = 
  #     res1 |> 
  #     dplyr::select(195:221) |> 
  #     ir::ir_drop_spectra() |> 
  #     dplyr::select(! dplyr::ends_with("_in_pd") & ! dplyr::all_of(c("eac_1", "edc_1"))) 
  # )
  
  
  res_plot <- 
    res1 |> 
    dplyr::select(id_sample, core_label, 195:226) |> 
    ir::ir_drop_spectra() |> 
    dplyr::select(! dplyr::ends_with("_in_pd") & ! dplyr::all_of(c("eac_1", "edc_1", "macroporosity_1"))) |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), as.numeric)
    ) |>
    tidyr::pivot_longer(
      ! dplyr::all_of(c("id_sample", "core_label", "degree_of_decomposition", "macrofossil_volume_fraction_sphagnum", "macrofossil_volume_fraction_sphagnum_hummock", "macrofossil_volume_fraction_sphagnum_lawn_hollow", "macrofossil_volume_fraction_sedge", "macrofossil_volume_fraction_shrub", "macrofossil_volume_fraction_uom")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "bulk_density_1" ~ "BD (g cm<sup>-3</sup>)",
          variable == "C_to_N_1" ~ "C/N (g g<sup>-1</sup>)",
          variable == "H_to_C_1" ~ "H/C (g g<sup>-1</sup>)",
          variable == "O_to_C_1" ~ "O/C (g g<sup>-1</sup>)",
          variable == "d13C_1" ~ "&delta;<sup>13</sup>C (---todo)",
          variable == "d15N_1" ~ "&delta;<sup>15</sup>N (---todo)",
          variable == "dgf0_1" ~ "&Delta;G<sub>f</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
          variable == "nosc_1" ~ "NOSC (-)",
          variable == "hi_1630_1090" ~ "HI<sub>1630/1090</sub>",
          TRUE ~ variable
        )
    ) |>
    dplyr::arrange(macrofossil_volume_fraction_sphagnum_lawn_hollow) |>
    ggplot(aes(y = value, x = degree_of_decomposition)) +
    geom_point(aes(color = macrofossil_volume_fraction_sphagnum_lawn_hollow)) +
    facet_grid(variable ~ core_label, scales = "free_y") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      strip.text.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown()
    ) +
    labs(
      y = "Decomposition indicator value",
      x = "Degree of decomposition (%)"
    ) +
    guides(
      color = guide_colorbar(title = "Macrofossil abudance of <i>Sphagnum</i>")
    )
  
  
}