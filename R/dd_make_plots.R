#' Plots true and predicted degree of decomposition for samples while adding mineral spectra
#' 
#' @export
dd_make_plot_1 <- function(dd_model, dd_data_model, dd_data_pmird_mineral, config, prediction_domain, file_plot) {
  
  # selection
  d_selection <- 
    dd_data_model |>
    dplyr::filter(
      (id_dataset == 17 & id_sample_original == 4727 & taxon_rank_value == "Larix laricina") |
        (id_dataset == 14 & id_sample_original == 4211 & taxon_rank_value == "Sphagnum rubellum") |
        (id_dataset == 1 & id_sample_original == 207 & taxon_rank_value == "Sphagnum capillifolium") |
        (id_dataset == 1 & id_sample_original == 142 & taxon_rank_value == "Typha latifolia") |
        (id_dataset == 24 & id_sample_original == 9246 & taxon_rank_value == "Phragmites australis" & taxon_organ == "rhizomes") |
        (id_dataset == 24 & id_sample_original == 9228 & taxon_rank_value == "Phragmites australis" & taxon_organ == "leaves")
    ) |>
    ir::ir_interpolate() |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
    
  dd_data_pmird_mineral <- 
    dd_data_pmird_mineral |>
    ir::ir_interpolate() |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
  
  # add minerals
  d_selection <- 
    d_selection |>
    dplyr::group_split(id_sample) |>
    purrr::map(function(.x) {
      scale_factor <- seq(0, 2, length.out = 20)
      .x |>
        dplyr::slice(rep(1, length(scale_factor))) |>
        magrittr::add(
          dd_data_pmird_mineral |>
            dplyr::slice(rep(1, length(scale_factor))) |>
            magrittr::multiply_by(scale_factor)
          ) |>
        dplyr::mutate(
          mineral_scale_factor = scale_factor
        )
    })
  
  # make predictions
  d_selection <- 
    d_selection |>
    purrr::map(function(.x) {
      irp_degree_of_decomposition(
        x = .x, 
        model = dd_model, 
        config = config, 
        prediction_domain = prediction_domain
      )
    })
  
  # make plot
  res <- 
    d_selection |>
    purrr::map(function(.x) {
      .x |>
        dplyr::select(id_sample, taxon_rank_value, taxon_organ, mass_relative_mass, degree_of_decomposition, degree_of_decomposition_in_pd, mineral_scale_factor)
    })
  res <- do.call("rbind", res)
  
  res_plot =
    res |>
    dplyr::mutate(
      facet_label = paste0("<i>", taxon_rank_value, "</i>", ifelse(is.na(taxon_organ), "", paste0(" (", taxon_organ, ")")))
    ) |>
    ggplot(aes(x = mineral_scale_factor, group = id_sample)) +
    ggdist::stat_lineribbon(aes(ydist = degree_of_decomposition)) +
    geom_point(
      data = res |> 
        dplyr::filter(! duplicated(id_sample)) |>
        dplyr::mutate(
          facet_label = paste0("<i>", taxon_rank_value, "</i>", ifelse(is.na(taxon_organ), "", paste0(" (", taxon_organ, ")")))
        ),
      mapping = aes(y = 1 - mass_relative_mass, group = id_sample),
      shape = 21,
      fill = "white"
    ) +
    guides(fill = guide_legend(title = "Prediction interval level")) +
    scale_fill_brewer() +
    facet_wrap(~ facet_label) +
    scale_y_continuous(labels = function(x) format(100 * x,digits = 2)) +
    labs(
      y = "Degree of decomposition (%)",
      x = "Relative contribution of minerals (-)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      legend.position = "bottom"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 7, height = 4.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots depth profiles of the degree of decomposition for cores from the Venner Moor 
#' 
#' @export
dd_make_plot_2 <- function(dd_model, dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_yhat, file_plot) {
  
  res <- 
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat,
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    ) |>
    dplyr::filter(id_dataset == 8)
  
  res_plot <- 
    res |>
    dplyr::select(core_label, sample_depth_upper, degree_of_decomposition, degree_of_decomposition_in_pd, nitrogen_content_1) |>
    dplyr::mutate(
      core_label =
        dplyr::case_when(
          core_label == "1" ~ "Drained, extracted,<br>restored since the 1970's",
          core_label == "4" ~ "Drained, not restored"
        ),
      nitrogen_content_1 = as.numeric(nitrogen_content_1)
    ) |>
    ggplot(aes(x = sample_depth_upper, group = core_label)) +
    ggdist::stat_lineribbon(aes(ydist = degree_of_decomposition)) +
    geom_point(aes(y = median(degree_of_decomposition), shape = degree_of_decomposition_in_pd, color = nitrogen_content_1 <= 0.02)) +
    guides(fill = guide_legend(title = "Prediction interval level")) +
    scale_fill_brewer() +
    facet_wrap(~ core_label) +
    scale_y_continuous(labels = function(x) format(100 * x, digits = 2)) +
    labs(
      y = "Degree of decomposition (%)",
      x = "Depth below the peat surface (cm)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.position = "bottom",
      legend.box="vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 4.5, height = 7, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots depth profiles of the degree of decomposition for cores from Patagonia 
#' 
#' @export
dd_make_plot_3 <- function(dd_model, dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_yhat, file_plot) {
  
  res <- 
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat,
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    ) |>
    dplyr::filter(core_label %in% c("PBR 2", "Sky 1", "Sky 2"))
  
  d_ash_layers <- 
    dplyr::bind_rows(
      tibble::tibble(
        core_label = "PBR 2",
        ash_depth_upper = c(35, 110, 155),
        ash_depth_lower = c(50, 125, 195)
      ),
      tibble::tibble(
        core_label = "Sky 1",
        ash_depth_upper = c(95, 135, 185),
        ash_depth_lower = c(105, 145, 195)
      ),
      tibble::tibble(
        core_label = "Sky 2",
        ash_depth_upper = c(120, 195, 245),
        ash_depth_lower = c(135, 210, 275)
      )
    )
    
  
  res_plot <- 
    res |>
    dplyr::filter(! duplicated(paste0(core_label, "_", sample_depth_upper))) |>
    dplyr::select(core_label, sample_depth_upper, degree_of_decomposition, degree_of_decomposition_in_pd, nitrogen_content_1, carbon_content_1, hi_1630_1090) |>
    dplyr::arrange(core_label, sample_depth_upper) |>
    dplyr::mutate(
      nitrogen_content_1 = as.numeric(nitrogen_content_1)
    ) |>
    ggplot(aes(group = core_label)) +
    geom_rect(
      data = d_ash_layers,
      aes(xmin = ash_depth_upper, xmax = ash_depth_lower, ymin = -Inf, ymax = Inf) ,
      fill = "lightgrey"
    ) +
    ggdist::stat_lineribbon(aes(x = sample_depth_upper, ydist = degree_of_decomposition), alpha = 0.6) +
    geom_point(aes(x = sample_depth_upper, y = median(degree_of_decomposition), shape = degree_of_decomposition_in_pd, color = nitrogen_content_1 <= 0.015)) +
    geom_path(aes(x = sample_depth_upper, y = as.numeric(carbon_content_1)), color = "brown") +
    geom_ribbon(
      aes(
        x = sample_depth_upper, 
        ymin = as.numeric(carbon_content_1 - errors::errors(carbon_content_1)),
        ymax = as.numeric(carbon_content_1 + errors::errors(carbon_content_1))
      ), 
      fill = "brown", 
      alpha = 0.6
    ) +
    geom_path(aes(x = sample_depth_upper, y = hi_1630_1090), color = "darksalmon") +
    guides(
      fill = guide_legend(title = "Prediction interval level"),
      shape = guide_legend(title = "Spectrum is in prediction domain")
    ) +
    scale_fill_brewer() +
    facet_wrap(~ core_label) +
    scale_y_continuous(labels = function(x) format(100 * x, digits = 2)) +
    labs(
      y = "Degree of decomposition (%), C<sub>MIRS</sub> (%), or HI<sub>1630/1090</sub> (times 100)",
      x = "Depth below the peat surface (cm)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      axis.title.x = ggtext::element_markdown(),
      legend.position = "bottom",
      legend.box = "vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6.5, height = 7, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots depth profiles of the difference in degree of decomposition predicted by two models for cores from the Venner Moor 
#' 
#' @export
dd_make_plot_4 <- function(dd_model, dd_data_pmird_peat_cores, dd_data_pmird_peat_cores_yhat, file_plot) {
  
  dd_data_pmird_peat_cores_yhat <- 
    dd_data_pmird_peat_cores_yhat[[1]] |>
    dplyr::mutate(
      degree_of_decomposition = degree_of_decomposition - dd_data_pmird_peat_cores_yhat[[2]]$degree_of_decomposition
    )
  
  res <- 
    dplyr::left_join(
      dd_data_pmird_peat_cores,
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    ) |>
    dplyr::filter(id_dataset == 8)
  
  res_plot <- 
    res |>
    dplyr::select(core_label, sample_depth_upper, degree_of_decomposition, degree_of_decomposition_in_pd) |>
    dplyr::mutate(
      core_label =
        dplyr::case_when(
          core_label == "1" ~ "Drained, extracted,<br>restored since the 1970's",
          core_label == "4" ~ "Drained, not restored"
        )
    ) |>
    ggplot(aes(x = sample_depth_upper, group = core_label)) +
    ggdist::stat_lineribbon(aes(ydist = degree_of_decomposition)) +
    geom_point(aes(y = median(degree_of_decomposition), shape = degree_of_decomposition_in_pd), color = "lightgrey") +
    guides(fill = guide_legend(title = "Prediction interval level")) +
    scale_fill_brewer() +
    facet_wrap(~ core_label) +
    scale_y_continuous(labels = function(x) format(100 * x, digits = 2)) +
    labs(
      y = "&Delta; Degree of decomposition (%)",
      x = "Depth below the peat surface (cm)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.position = "bottom",
      legend.box="vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 4.5, height = 7, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots true and predicted degree of decomposition for samples while adding a decomposed peat spectrum
#' 
#' @export
dd_make_plot_5 <- function(dd_model, dd_data_model, dd_data_pmird_decomposed, dd_data_pmird_peat_cores, config, prediction_domain, file_plot) {
  
  # selection
  d_selection <- 
    dd_data_model |>
    dplyr::filter(
      (id_dataset == 17 & id_sample_original == 4727 & taxon_rank_value == "Larix laricina") |
        (id_dataset == 14 & id_sample_original == 4211 & taxon_rank_value == "Sphagnum rubellum") |
        (id_dataset == 1 & id_sample_original == 207 & taxon_rank_value == "Sphagnum capillifolium") |
        (id_dataset == 1 & id_sample_original == 142 & taxon_rank_value == "Typha latifolia") |
        (id_dataset == 24 & id_sample_original == 9246 & taxon_rank_value == "Phragmites australis" & taxon_organ == "rhizomes") |
        (id_dataset == 24 & id_sample_original == 9228 & taxon_rank_value == "Phragmites australis" & taxon_organ == "leaves")
    ) |>
    ir::ir_interpolate() |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
  
  dd_data_pmird_mineral <- 
    dd_data_pmird_decomposed |>
    ir::ir_interpolate() |>
    ir::ir_clip(range = data.frame(start = 600, end = 4000)) |>
    ir::ir_bc(method = "rubberband", do_impute = TRUE) |>
    ir::ir_normalise(method = "area")
  
  # add minerals
  d_selection <- 
    d_selection |>
    dplyr::group_split(id_sample) |>
    purrr::map(function(.x) {
      scale_factor <- seq(0, 2, length.out = 20)
      .x |>
        dplyr::slice(rep(1, length(scale_factor))) |>
        magrittr::add(
          dd_data_pmird_mineral |>
            dplyr::slice(rep(1, length(scale_factor))) |>
            magrittr::multiply_by(scale_factor)
        ) |>
        dplyr::mutate(
          mineral_scale_factor = scale_factor
        )
    })
  
  # make predictions
  d_selection <- 
    d_selection |>
    purrr::map(function(.x) {
      irp_degree_of_decomposition(
        x = .x, 
        model = dd_model, 
        config = config, 
        prediction_domain = prediction_domain
      )
    })
  
  # make plot
  res <- 
    d_selection |>
    purrr::map(function(.x) {
      .x |>
        dplyr::select(id_sample, taxon_rank_value, taxon_organ, mass_relative_mass, degree_of_decomposition, degree_of_decomposition_in_pd, mineral_scale_factor)
    })
  res <- do.call("rbind", res)
  
  res_plot =
    res |>
    dplyr::mutate(
      facet_label = paste0("<i>", taxon_rank_value, "</i>", ifelse(is.na(taxon_organ), "", paste0(" (", taxon_organ, ")")))
    ) |>
    ggplot(aes(x = mineral_scale_factor, group = id_sample)) +
    ggdist::stat_lineribbon(aes(ydist = degree_of_decomposition)) +
    geom_point(
      data = res |> 
        dplyr::filter(! duplicated(id_sample)) |>
        dplyr::mutate(
          facet_label = paste0("<i>", taxon_rank_value, "</i>", ifelse(is.na(taxon_organ), "", paste0(" (", taxon_organ, ")")))
        ),
      mapping = aes(y = 1 - mass_relative_mass, group = id_sample),
      shape = 21,
      fill = "white"
    ) +
    guides(fill = guide_legend(title = "Prediction interval level")) +
    scale_fill_brewer() +
    facet_wrap(~ facet_label) +
    scale_y_continuous(labels = function(x) format(100 * x,digits = 2)) +
    labs(
      y = "Degree of decomposition (%)",
      x = "Relative contribution of decomposed peat sample (-)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      legend.position = "bottom"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 7, height = 4.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots various decomposition indicators versus the suggested decomposition indicator
#' 
#' @export
dd_make_plot_6 <- function(dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_macrofossils, dd_data_pmird_peat_cores_yhat, variable_color, variable_color_legend_title = "Macrofossil abudance of <i>Sphagnum</i>", file_plot) {
  
  dd_data_pmird_peat_cores_macrofossils <- 
    dd_data_pmird_peat_cores_macrofossils |>
    dplyr::filter(stringr::str_detect(taxon_rank_value, "Sphagnum") | stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum|Ericaceae|Cyperaceae)") | macrofossil_type %in% c("wood", "amorphous organic matter")) |>
    dplyr::select(core_label, sample_depth_upper, sample_depth_lower, taxon_rank_value, macrofossil_type, macrofossil_volume_fraction) |>
    dplyr::mutate(
      taxon_rank_value = 
        dplyr::case_when(
          stringr::str_detect(taxon_rank_value, "Sphagnum") & stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_hummock",
          stringr::str_detect(taxon_rank_value, "Sphagnum") & ! stringr::str_detect(taxon_rank_value, "(fuscum|capillifolium)") ~ "sphagnum_lawn_hollow",
          stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum|Cyperaceae)") ~ "sedge",
          stringr::str_detect(taxon_rank_value, "Ericaceae")  ~ "shrub",
          macrofossil_type == "wood"  ~ "wood",
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
      dd_data_pmird_peat_cores_irpeat |>
        ir::ir_drop_spectra(),
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
        }),
      macrofossil_volume_fraction_wood =
        purrr::map_dbl(seq_along(id_sample), function(i) {
          index <- dd_data_pmird_peat_cores_macrofossils$core_label == core_label[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_upper >= sample_depth_upper[[i]] & dd_data_pmird_peat_cores_macrofossils$sample_depth_lower <= sample_depth_lower[[i]]
          if(sum(index) == 1L) {
            dd_data_pmird_peat_cores_macrofossils$wood[index]
          } else if(sum(index) > 1L) {
            mean(dd_data_pmird_peat_cores_macrofossils$wood[index], na.rm = TRUE)
          } else {
            NA_real_
          }
        })
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::ends_with("_in_pd"), as.logical)
    ) 
  
  res <- 
    res |> 
    dplyr::mutate(
      degree_of_decomposition_lower = posterior::quantile2(degree_of_decomposition, probs = 0.25),
      degree_of_decomposition_upper = posterior::quantile2(degree_of_decomposition, probs = 0.75),
      degree_of_decomposition = median(degree_of_decomposition)
    )
  
  res <- 
    res |> 
    dplyr::select(id_sample, core_label, C, N, d13C, d15N, bulk_density, 195:232) |> 
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), function(.x) {
        if(inherits(.x, "quantities")) {
          units::drop_units(.x)
        } else if (inherits(.x, "rvars")) {
          errors::set_errors(mean(.x), sd(.x))
        } else {
          .x
        }
      }),
      dplyr::across(
        dplyr::starts_with("macrofossil_volume_fraction_"),
        function(.x) .x * 100
      ),
      bulk_density_1 =
        dplyr::case_when(
          ! is.na(bulk_density) ~ bulk_density,
          TRUE ~ bulk_density_1
        ),
      carbon_content_1 =
        dplyr::case_when(
          ! is.na(C) ~ C,
          TRUE ~ carbon_content_1
        ),
      C_to_N_1 =
        dplyr::case_when(
          ! is.na(C/N) ~ C/N,
          TRUE ~ C_to_N_1
        ),
      d13C_1 =
        dplyr::case_when(
          ! is.na(d13C) ~ d13C,
          TRUE ~ d13C_1
        ),
      d15N_1 =
        dplyr::case_when(
          ! is.na(d15N) ~ d15N,
          TRUE ~ d15N_1
        ),
      core_label =
        dplyr::case_when(
          core_label == "1" ~ "Venner Moor<br>(rewetted)",
          core_label == "4" ~ "Venner Moor<br>(forested)",
          core_label == "eb1005_MH1" ~ "MH 1",
          core_label == "eb1006_MK1" ~ "MK 1",
          core_label == "eb1018_OD2" ~ "OD 2",
          TRUE ~ core_label
        ) |>
        factor(levels = c("Venner Moor<br>(rewetted)", "Venner Moor<br>(forested)", "OD 2", "MH 1", "MK 1", "PBR 2", "Sky 1", "Sky 2"))
    ) |>
    dplyr::select(! dplyr::ends_with("_in_pd") & ! dplyr::any_of(c("eac_1", "edc_1", "macroporosity_1", "d13C", "d15N", "C", "N", "bulk_density"))) |>
    tidyr::pivot_longer(
      ! dplyr::all_of(c("id_sample", "core_label", "macrofossil_volume_fraction_sphagnum", "macrofossil_volume_fraction_sphagnum_hummock", "macrofossil_volume_fraction_sphagnum_lawn_hollow", "macrofossil_volume_fraction_sedge", "macrofossil_volume_fraction_shrub", "macrofossil_volume_fraction_wood", "macrofossil_volume_fraction_uom", "carbon_content_1", "nitrogen_content_1")) & ! dplyr::starts_with("degree_of_decomposition"),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::left_join(
      res |>
        dplyr::select(id_sample, dplyr::ends_with("_in_pd"), C, N, d13C, d15N, bulk_density) |>
        dplyr::mutate(
          bulk_density_1_in_pd = 
            dplyr::case_when(
              ! is.na(bulk_density) ~ NA,
              TRUE ~ bulk_density_1_in_pd
            ),
          C_to_N_1_in_pd = 
            dplyr::case_when(
              ! is.na(C/N) ~ NA,
              TRUE ~ C_to_N_1_in_pd
            ),
          d13C_1_in_pd =
            dplyr::case_when(
              ! is.na(d13C) ~ NA,
              TRUE ~ d13C_1_in_pd
            ),
          d15N_1_in_pd =
            dplyr::case_when(
              ! is.na(d15N) ~ NA,
              TRUE ~ d15N_1_in_pd
            )
        ) |>
        tidyr::pivot_longer(
          cols = ! dplyr::any_of(c("id_sample", "C", "N", "d13C", "d15N", "bulk_density")),
          names_to = "variable",
          values_to = "is_in_pd"
        ) |>
        dplyr::select(id_sample, variable, is_in_pd) |>
        dplyr::mutate(
          variable =
            variable |>
            stringr::str_remove("_in_pd$"),
          is_in_pd =
            dplyr::case_when(
              is.na(is_in_pd) ~ "measured",
              is_in_pd ~ "in prediction domain",
              ! is_in_pd ~ "not in prediction domain"
            )
        ),
      by = c("id_sample", "variable")
    ) |>
    dplyr::mutate(
      variable_or = variable, 
      variable =
        dplyr::case_when(
          variable == "bulk_density_1" ~ "BD (g cm<sup>-3</sup>)",
          variable == "C_to_N_1" ~ "C/N (g g<sup>-1</sup>)",
          variable == "H_to_C_1" ~ "H/C (g g<sup>-1</sup>)",
          variable == "O_to_C_1" ~ "O/C (g g<sup>-1</sup>)",
          variable == "d13C_1" ~ "&delta;<sup>13</sup>C (&permil;)",
          variable == "d15N_1" ~ "&delta;<sup>15</sup>N (&permil;)",
          variable == "dgf0_1" ~ "&Delta;G<sub>f</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
          variable == "nosc_1" ~ "NOSC (-)",
          variable == "hi_1630_1090" ~ "HI<sub>1630/1090</sub>",
          TRUE ~ variable
        ),
      variable = 
        factor(variable, levels = unique(variable)),
      is_in_pd =
        dplyr::case_when(
          is.na(is_in_pd) ~ "measured",
          TRUE ~ is_in_pd
        ) |>
        factor(levels = c("measured", "in prediction domain", "not in prediction domain"))
    )
  
  res_correlation <- 
    res |>
    dplyr::group_by(variable_or) |>
    dplyr::summarise(
      correlation = cor(degree_of_decomposition, value),
      .groups = "drop"
    )
  
  res_plot <- 
    res |>
    dplyr::left_join(
      res_correlation,
      by = "variable_or"
    ) |>
    dplyr::arrange(variable) |>
    ggplot(aes(y = value, x = degree_of_decomposition, group = core_label)) +
    geom_errorbar(aes(ymin = value - errors::errors(value), ymax = value + errors::errors(value)), color = "lightgrey", width = 0) +
    geom_errorbarh(aes(xmin = degree_of_decomposition_lower, xmax = degree_of_decomposition_upper), color = "lightgrey") +
    geom_point(aes(color = as.numeric(eval(variable_color)), shape = is_in_pd)) +
    scale_x_continuous(labels = function(x) round(x * 100, 0)) +
    scale_color_binned(breaks = seq(0, 100, by = 20)) +
    facet_grid(paste0(variable, "<br>&rho; = ", round(correlation, digits = 2L)) |> factor(levels = unique(paste0(variable, "<br>&rho; = ", round(correlation, digits = 2L)))) ~ core_label, scales = "free_y") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 15),
      strip.text.y = ggtext::element_markdown(size = 12, lineheight = 1.5),
      legend.title = ggtext::element_markdown(size = 16),
      axis.title.x = ggtext::element_markdown(size = 18),
      axis.title.y = ggtext::element_markdown(size = 18),
      legend.text = ggtext::element_markdown(size = 15),
      legend.box = "horizontal",
      legend.direction = "horizontal"
    ) +
    labs(
      y = "Decomposition indicator value",
      x = "Degree of decomposition (%)"
    ) +
    guides(
      color = 
        guide_legend(
          title = variable_color_legend_title, 
          override.aes = list(size = 4)#,
          #theme = 
          #  theme(
          #    legend.key.width = unit(dev.size()[1] / 17, "inches"),
          #    legend.key.height = unit(1, "lines"),
          #    legend.box = "vertical",
          #    legend.direction = "horizontal"
          #)
        ),
      shape = guide_legend(title = "Reliability", override.aes = list(size = 4))
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 13.3, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots reconstructed initial C/N ratios versus depth for the cores
#' 
#' @export
dd_make_plot_7 <- function(dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_yhat, file_plot) {
  
  res <- 
    dplyr::left_join(
      dd_data_pmird_peat_cores_irpeat |>
        ir::ir_drop_spectra(),
      dd_data_pmird_peat_cores_yhat,
      by = "id_sample"
    )
  
  res |>
    dplyr::mutate(
      N =
        dplyr::case_when(
          is.na(N) ~ as.numeric(nitrogen_content_1),
          TRUE ~ N
        ),
      y = (0.4/(1 - degree_of_decomposition))/N
    ) |>
    dplyr::filter(! duplicated(paste0(core_label, "_", sample_depth_upper)) & mean(y) < 1e3) |>
    ggplot(aes(x = sample_depth_lower)) +
    ggdist::stat_lineribbon(aes(ydist = y, group = core_label)) +
    guides(fill = guide_legend(title = "Prediction interval level")) +
    scale_fill_brewer() +
    facet_wrap(~ core_label) +
    labs(
      y = "Degree of decomposition (%)",
      x = "Depth below the peat surface (cm)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      legend.position = "bottom",
      legend.box="vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
}



#' Plots model residuals versus N_MIRS conditional on the degree of decomposition for the training data
#' 
#' @export
dd_make_plot_8 <- function(dd_stan_1_fit, dd_data_model, file_plot) {
  
  dd_data_model <- 
    dd_data_model |>
    irpeat::irp_nitrogen_content_1(do_summary = TRUE) |>
    dplyr::mutate(
      is_in_pd =
        dplyr::case_when(
          ! is.na(N) ~ "measured",
          is.na(N) & nitrogen_content_1_in_pd ~ "in prediction domain",
          is.na(N) & ! nitrogen_content_1_in_pd ~ "not in prediction domain"
        ) |>
        factor(levels = c("measured", "in prediction domain", "not in prediction domain")),
      nitrogen_content_1 =
        dplyr::case_when(
          is.na(N) ~ nitrogen_content_1,
          TRUE ~ quantities::set_quantities(N, unit = "g/g", errors = N_err, mode = "standard")
        )
    )
  
  res <- 
    purrr::map2(dd_stan_1_fit, seq_along(dd_stan_1_fit), function(.x, i) {
      dd_data_model |>
        dplyr::select(dd_dataset_label, nitrogen_content_1, is_in_pd, mass_relative_mass, taxon_rank_value, taxon_organ) |>
        dplyr::mutate(
          yhat_residuals = residuals(.x)[, 1],
          id_model = i,
          degree_of_decomposition = 1 - mass_relative_mass
        )
    })
  res <- do.call("rbind", res)
  
  res_plot <- 
    res |>
    dplyr::mutate(
      taxon_rank_value =
        dplyr::case_when(
          taxon_rank_value %in% c("Sphagnum capillifolium", "Phragmites australis", "Typha latifolia") & dd_dataset_label != "pmird_undecomposed_litter" ~ taxon_rank_value,
          TRUE ~ "Others"
        ),
      plot_label =
        dplyr::case_when(
          taxon_rank_value == "Others" ~ taxon_rank_value,
          TRUE ~ paste0("<i>", taxon_rank_value, "</i>"," (", taxon_organ, ")")
        ) |>
        factor(levels = c("<i>Phragmites australis</i> (leaves)", "<i>Phragmites australis</i> (rhizomes)", "<i>Typha latifolia</i> (aboveground parts)", "<i>Sphagnum capillifolium</i> (whole plant)", "Others")),
      facet_label =
        dplyr::case_when(
          degree_of_decomposition == 0 ~ "Undecomposed",
          TRUE ~ cut_interval(degree_of_decomposition * 100, length = 20)
        ),
      facet_label =
        dplyr::case_when(
          facet_label %in% c("(60,80]", "(80,100]") ~ "(60,100]",
          TRUE ~ facet_label
        ),
      facet_label = factor(facet_label, levels = c("Undecomposed", unique(facet_label)[-1][order(unique(facet_label)[-1] |> stringr::str_extract(pattern = "\\d{1}\\.?\\d*,") |> stringr::str_remove(",$") |> as.numeric())]))
    ) |>
    ggplot(aes(y = yhat_residuals, x = as.numeric(nitrogen_content_1))) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_point(aes(color = plot_label, shape = is_in_pd)) +
    geom_smooth(method = "lm", color = "grey50", formula = 'y ~ x') +
    facet_grid(id_model ~ facet_label) +
    guides(
      color = guide_legend(title = "Species (organ)", nrow = 3L, byrow = TRUE, override.aes = list(size = 4), order = 1),
      shape = guide_legend(title = "Reliability", override.aes = list(size = 4, ncol = 1, order = 2))
    ) +
    scale_color_manual(values = c("lightsalmon3", "lightsalmon1", "lightblue", "darkseagreen", "grey70")) +
    scale_y_continuous(labels = function(x) x * 100) +
    labs(
      y = "Degree of decomposition residuals (%)",
      x = "N or N<sub>MIRS</sub> (g g<sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 11),
      strip.text.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(size = 11),
      axis.title.y = ggtext::element_markdown(size = 11),
      legend.title = ggtext::element_markdown(size = 11),
      legend.text = ggtext::element_markdown(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "vertical"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 6.3, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots model coefficients
#' 
#' @export
dd_make_plot_9 <- function(dd_stan_1_fit, dd_model_info, file_plot) {
  
  
  
  res <- 
    purrr::map(seq_along(dd_stan_1_fit), function(i) {
      .x <- dd_stan_1_fit[[i]]
      parnames <- brms::variables(.x)
      .x |>
        tidybayes::gather_rvars(!!!syms(parnames[stringr::str_detect(parnames, "^b_x")])) |>
        dplyr::rename(
          variable = ".variable",
          value = ".value"
        ) |>
        dplyr::mutate(
          id_model = i,
          x = 
            dd_model_info$x[[i]]$spectra[[1]]$x,
          y = 
            dd_model_info$x[[i]]$spectra[[1]]$y - min(dd_model_info$x[[i]]$spectra[[1]]$y) + 1
        )
    }) |>
    dplyr::bind_rows()
  
  res_plot <- 
    res |>
    ggplot(aes(x = x)) +
    geom_path(aes(y = y), color = "grey") +
    ggdist::stat_lineribbon(aes(ydist = value), linewidth = 0) +
    scale_fill_brewer() +
    geom_hline(yintercept = 0, color = "grey30") +
    facet_wrap(~ id_model, ncol = 1L) +
    guides(
      x = guide_axis(minor.ticks = TRUE),
      y = guide_axis(minor.ticks = TRUE)
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      strip.text.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown()
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Combines plot `dd_plot_2` and `dd_plot_3` into one plot for multiple models
#' 
#' @export
dd_make_plot_10 <- function(dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_yhat, dd_model_info, file_plot, ggsave_width = 11, ggsave_height = 8) {
                            
  # dd_plot_2
  res <- 
    purrr::map(seq_along(dd_data_pmird_peat_cores_yhat), function(i) {
      dplyr::left_join(
        dd_data_pmird_peat_cores_irpeat,
        dd_data_pmird_peat_cores_yhat[[i]],
        by = "id_sample"
      ) |>
        ir::ir_drop_spectra() |>
        dplyr::filter(id_dataset == 8) |>
        dplyr::mutate(
          id_model = dd_model_info$id_model[[i]],
          dplyr::across(dplyr::ends_with("_in_pd"), as.logical)
        )
    }) |>
    dplyr::bind_rows()
    
  res_plot1 <- 
    res |>
    dplyr::mutate(
      C = quantities::set_quantities(C, unit = "g/g", errors = C_err, mode = "standard")
    ) |>
    dplyr::select(id_model, core_label, sample_depth_upper, degree_of_decomposition, degree_of_decomposition_in_pd, nitrogen_content_1, carbon_content_1, hi_1630_1090, N, C) |>
    dplyr::arrange(core_label, sample_depth_upper) |>
    dplyr::mutate(
      nitrogen_content_1 = as.numeric(nitrogen_content_1),
      nitrogen_content_1 = 
        dplyr::case_when(
          is.na(N) ~ nitrogen_content_1,
          TRUE ~ N
        ),
      carbon_content_1 =
        dplyr::case_when(
          is.na(N) ~ carbon_content_1,
          TRUE ~ C
        ),
      point_fill = factor(nitrogen_content_1 < 0.02, levels = c(TRUE, FALSE))
    ) |>
    dplyr::mutate(
      core_label =
        dplyr::case_when(
          core_label == "1" ~ "Drained, extracted,<br>restored since the 1970's",
          core_label == "4" ~ "Drained, not restored"
        )
    ) |>
    ggplot(aes(x = sample_depth_upper, group = core_label)) +
    ggdist::stat_lineribbon(aes(ydist = degree_of_decomposition), alpha = 0.6, linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Prediction interval level", order = 1, override.aes = list(linewidth = 0, color = NA))
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(y = median(degree_of_decomposition), shape = degree_of_decomposition_in_pd, fill = point_fill)) +
    scale_fill_manual(values = c("white", "grey20")) +
    scale_shape_manual(values = c(21, 24)) +
    guides(
      fill = guide_legend(title = "N < 0.02 g g<sup>-1</sup>", order = 2, override.aes = list(size = 4, shape = 21)),
      shape = guide_legend(title = "Sample is in prediction domain", order = 3, override.aes = list(size = 4))
    ) +
    geom_path(aes(x = sample_depth_upper, y = as.numeric(carbon_content_1), color = "C or C<sub>MIRS</sub>")) +
    geom_ribbon(
      aes(
        x = sample_depth_upper, 
        ymin = as.numeric(carbon_content_1 - errors::errors(carbon_content_1)),
        ymax = as.numeric(carbon_content_1 + errors::errors(carbon_content_1))
      ), 
      fill = "brown", 
      alpha = 0.6
    ) +
    geom_path(aes(x = sample_depth_upper, y = hi_1630_1090, color = "HI<sub>1630/1090</sub>")) +
    scale_color_manual(values = c("brown", "darksalmon")) +
    guides(
      color = guide_legend(title = "Other peat properties", order = 4, override.aes = list(linewidth = 2))
    ) +
    facet_grid(id_model ~ core_label) +
    scale_y_continuous(labels = function(x) format(100 * x, digits = 2)) +
    labs(
      y = "&gamma;<sub>MIRS</sub> (%), C or C<sub>MIRS</sub> (%), or HI<sub>1630/1090</sub> (times 100)",
      x = "Depth below the peat surface (cm)",
      title = "Venner Moor"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 11),
      strip.text.y = ggtext::element_markdown(size = 11),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.text = ggtext::element_markdown(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
  
  ## dd_plot_3
  res <- 
    purrr::map(seq_along(dd_data_pmird_peat_cores_yhat), function(i) {
      dplyr::left_join(
        dd_data_pmird_peat_cores_irpeat,
        dd_data_pmird_peat_cores_yhat[[i]],
        by = "id_sample"
      ) |>
        ir::ir_drop_spectra() |>
        dplyr::filter(core_label %in% c("PBR 2", "Sky 1", "Sky 2")) |>
        dplyr::mutate(
          id_model = dd_model_info$id_model[[i]],
          dplyr::across(dplyr::ends_with("_in_pd"), as.logical)
        )
    }) |>
    dplyr::bind_rows()
  
  d_ash_layers <- 
    dplyr::bind_rows(
      tibble::tibble(
        core_label = "PBR 2",
        ash_depth_upper = c(35, 110, 155),
        ash_depth_lower = c(50, 125, 195)
      ),
      tibble::tibble(
        core_label = "Sky 1",
        ash_depth_upper = c(95, 135, 185),
        ash_depth_lower = c(105, 145, 195)
      ),
      tibble::tibble(
        core_label = "Sky 2",
        ash_depth_upper = c(120, 195, 245),
        ash_depth_lower = c(135, 210, 275)
      )
    )
  
  
  res_plot2 <- 
    res |>
    dplyr::filter(! duplicated(paste0(id_model, "_", core_label, "_", sample_depth_upper))) |>
    dplyr::select(id_model, core_label, sample_depth_upper, degree_of_decomposition, degree_of_decomposition_in_pd, nitrogen_content_1, carbon_content_1, N, C, hi_1630_1090) |>
    dplyr::arrange(core_label, sample_depth_upper) |>
    dplyr::mutate(
      nitrogen_content_1 = as.numeric(nitrogen_content_1),
      nitrogen_content_1 = 
        dplyr::case_when(
          is.na(N) ~ nitrogen_content_1,
          TRUE ~ N
        ),
      carbon_content_1 =
        dplyr::case_when(
          is.na(N) ~ carbon_content_1,
          TRUE ~ quantities::set_quantities(C, "g/g", 0, mode = "standard")
        ),
      point_fill = factor(nitrogen_content_1 < 0.02, levels = c(TRUE, FALSE))
    ) |>
    dplyr::bind_rows( #---note: dummy entry to make the legend for N content the same across all plots
      tibble::tibble(
        core_label = "Sky 1",
        sample_depth_upper = 0,
        nitrogen_content_1 = 0.03,
        point_fill = factor(nitrogen_content_1 < 0.02, levels = c(TRUE, FALSE)),
        degree_of_decomposition_in_pd = TRUE,
        id_model = unique(res$id_model)[[1]]
      )
    ) |>
    ggplot(aes(group = core_label)) +
    geom_rect(
      data = d_ash_layers,
      aes(xmin = ash_depth_upper, xmax = ash_depth_lower, ymin = -Inf, ymax = Inf) ,
      fill = "lightgrey"
    ) +
    ggdist::stat_lineribbon(aes(x = sample_depth_upper, ydist = degree_of_decomposition), alpha = 0.6, linewidth = 0.5) +
    scale_fill_brewer() +
    guides(
      fill = guide_legend(title = "Prediction interval level", order = 1, override.aes = list(linewidth = 0, color = NA))
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(x = sample_depth_upper, y = median(degree_of_decomposition), shape = degree_of_decomposition_in_pd, fill = point_fill)) +
    scale_fill_manual(values = c("white", "grey20")) +
    scale_shape_manual(values = c(21, 24)) +
    guides(
      fill = guide_legend(title = "N < 0.02 g g<sup>-1</sup>", order = 2, override.aes = list(size = 4, shape = 21)),
      shape = guide_legend(title = "Sample is in prediction domain", order = 3, override.aes = list(size = 4))
    ) +
    geom_path(aes(x = sample_depth_upper, y = as.numeric(carbon_content_1), color = "C or C<sub>MIRS</sub>")) +
    geom_ribbon(
      aes(
        x = sample_depth_upper, 
        ymin = as.numeric(carbon_content_1 - errors::errors(carbon_content_1)),
        ymax = as.numeric(carbon_content_1 + errors::errors(carbon_content_1))
      ), 
      fill = "brown", 
      alpha = 0.6
    ) +
    geom_path(aes(x = sample_depth_upper, y = hi_1630_1090, color = "HI<sub>1630/1090</sub>")) +
    scale_color_manual(values = c("brown", "darksalmon")) +
    guides(
      color = guide_legend(title = "Other peat properties", order = 4, override.aes = list(linewidth = 2))
    ) +
    facet_grid(id_model ~ core_label) +
    scale_y_continuous(labels = function(x) format(100 * x, digits = 2)) +
    labs(
      y = "&gamma;<sub>MIRS</sub> (%), C or C<sub>MIRS</sub> (%), or HI<sub>1630/1090</sub> (times 100)",
      x = "Depth below the peat surface (cm)",
      title = "Patagonian peatlands"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 11),
      strip.text.y = ggtext::element_markdown(size = 11),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.text = ggtext::element_markdown(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "vertical"
    ) +
    coord_flip() +
    scale_x_reverse()
  
  res_plot <- 
    patchwork::wrap_plots(list(res_plot1, res_plot2)) +
    patchwork::plot_layout(ncol = 2L, nrow = 1L, widths = c(1, 1.5), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      legend.position = "bottom"
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = ggsave_width, 
    height = ggsave_height, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Boxplot of the suggested decomposition indicator and HI_1630/1090 for undecomposed litter
#' 
#' @param kfold An list like `dd_stan_1_kfold_2`.
#' 
#' @export
dd_make_plot_11 <- function(dd_data_model_irpeat, kfold, dd_mir_preprocessing_config, dd_model_info, dd_data_pmird_peat_cores_irpeat, file_plot) {
  
  res <- 
    purrr::map(seq_len(nrow(dd_model_info)), function(i) {
      dd_data_model_irpeat |>
        dplyr::mutate(
          degree_of_decomposition = 
            brms::kfold_predict(kfold[[i]])$yrep |>
            posterior::as_draws_df() |>
            posterior::as_draws_rvars(),
          degree_of_decomposition = do.call("c", degree_of_decomposition)
        ) |>
        ir::ir_drop_spectra() |>
        dplyr::mutate(
          dplyr::across(dplyr::ends_with("_in_pd"), as.logical),
          id_model = dd_model_info$id_model[[i]],
          degree_of_decomposition_measured = 1 - mass_relative_mass
        ) |>
        dplyr::filter(mass_relative_mass == 1.0) |>
        dplyr::select(id_model, taxon_rank_value, taxon_organ, degree_of_decomposition_measured, degree_of_decomposition, nitrogen_content_1, hi_1630_1090)
    }) |>
    dplyr::bind_rows() |>
    tidyr::pivot_longer(
      dplyr::all_of(c("degree_of_decomposition_measured", "degree_of_decomposition", "hi_1630_1090")),
      names_to = "variable",
      values_to = "value"
    )
  
  res_peat <- 
    dd_data_pmird_peat_cores_irpeat |>
    dplyr::filter(C > 0.3) |>
    dplyr::summarise(
      hi_1630_1090_min = min(hi_1630_1090, na.rm = TRUE),
      hi_1630_1090_max = max(hi_1630_1090, na.rm = TRUE)
    )
  
  res <- 
    res |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "degree_of_decomposition_measured" ~ "Degree of decomposition (measured)",
          variable == "degree_of_decomposition" ~ "Degree of decomposition (predicted)",
          variable == "hi_1630_1090" ~ "HI<sub>1630/1090</sub>"
        ),
      taxon_rank_value =
        dplyr::case_when(
          is.na(taxon_organ) ~ paste0("<i>", taxon_rank_value, "</i>"),
          TRUE ~ paste0("<i>", taxon_rank_value, "</i> (", taxon_organ, ")")
        ),
      taxon_rank_value =
        factor(
          taxon_rank_value, 
          levels = 
            tibble::tibble(
              taxon_rank_value = taxon_rank_value,
              nitrogen_content_1 = nitrogen_content_1 
            ) |>
            dplyr::group_by(taxon_rank_value) |>
            dplyr::summarise(
              nitrogen_content_1 = median(nitrogen_content_1),
              .groups = "drop"
            ) |>
            dplyr::bind_rows(
              tibble::tibble(
                taxon_rank_value = "Peat",
                nitrogen_content_1 = quantities::set_quantities(200, g/g, 0)
              )
            ) |>
            dplyr::arrange(nitrogen_content_1) |>
            dplyr::pull(taxon_rank_value)
        )
    )
  
  res_plot <- 
    res |>
    ggplot(aes(x = taxon_rank_value)) +
    ggdist::stat_pointinterval(aes(ydist = value, color = variable, size = as.numeric(nitrogen_content_1)), linewidth = 0.5) +
    geom_errorbar(
      data = 
        res_peat |>
        dplyr::mutate(
          taxon_rank_value = factor("Peat", levels = levels(res$taxon_rank_value))
        ),
      aes(x = taxon_rank_value, ymin = hi_1630_1090_min, ymax = hi_1630_1090_max),
      color = "darksalmon"
    ) +
    facet_wrap(~ id_model) +
    labs(
      y = "Degree of decomposition (%) or HI<sub>1630/1090</sub> (times 100)",
      x = "Taxon (organ)"
    ) +
    guides(
      color = guide_legend(title = "Variable", ncol = 2L)
    ) +
    scale_color_manual(values = c("steelblue1", "steelblue3", "darksalmon")) +
    scale_y_continuous(labels = function(x) x * 100) +
    coord_flip() +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 12),
      strip.text.y = ggtext::element_markdown(size = 12),
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggtext::element_markdown(),
      axis.text.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.text = ggtext::element_markdown(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal"
    ) 
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 7.5, height = 7.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Scatterplot of other decomposition indicators versus the measured degree of decomposition for litterbag samples
#' 
#' @param dd_stan_1_fit One element of `dd_stan_1_fit`.
#' 
#' @export
dd_make_plot_12 <- function(dd_data_model_irpeat, dd_stan_1_fit, file_plot) {
  
  res <- 
    dd_data_model_irpeat |>
    dplyr::bind_cols(
      predict(dd_stan_1_fit) |>
        as.data.frame() |>
        dplyr::mutate(
          degree_of_decomposition_1 = quantities::set_quantities(Estimate * 100, unit = "1", errors = `Est.Error` * 100, mode = "standard")
        ) |>
        dplyr::select(degree_of_decomposition_1) |>
        dplyr::mutate(
          degree_of_decomposition_1_in_pd = TRUE
        )
    ) |>
    dplyr::filter(taxon_rank_value %in% c("Sphagnum capillifolium", "Typha latifolia", "Phragmites australis") & dd_dataset_label %in% c("reuter2019", "arsenault2024a")) |>
    dplyr::mutate(
      nitrogen_supply_growth =
        dplyr::case_when(
          dd_dataset_label != "reuter2019" ~ NA_character_,
          dd_dataset_label == "reuter2019" & taxon_organ == "leaves" ~ NA_character_,
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^K") ~ "high-N",
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^T") ~ "low-N",
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^S") ~ "medium-N"
        ),
      core_label = "Litter",
      degree_of_decomposition = 1 - mass_relative_mass,
      macrofossil_volume_fraction_sphagnum_hummock =
        dplyr::case_when(
          taxon_rank_value %in% c("Sphagnum fuscum", "Sphagnum capillifolium") ~ 1,
          TRUE ~ 0
        ),
      macrofossil_volume_fraction_sphagnum_lawn_hollow =
        dplyr::case_when(
          ! taxon_rank_value %in% c("Sphagnum fuscum", "Sphagnum capillifolium") & stringr::str_detect(taxon_rank_value, "^Sphagnum") ~ 1,
          TRUE ~ 0
        ),
      macrofossil_volume_fraction_sphagnum =
        macrofossil_volume_fraction_sphagnum_lawn_hollow + macrofossil_volume_fraction_sphagnum_hummock,
      macrofossil_volume_fraction_sedge =
        dplyr::case_when(
          stringr::str_detect(taxon_rank_value, "(Carex|Sedge|Eriophorum|Erophorum|Cyperaceae)") ~ 1,
          TRUE ~ 0
        ),
      macrofossil_volume_fraction_shrub =
        dplyr::case_when(
          taxon_rank_value %in% c("Empetrum rubrum", "Lepidothamnus fonkii", "Gaultheria antarctica", "Calluna vulgaris", "Vaccinium myrtillus", "Vaccinium uliginosum", "Empetrum nigrum", "Vaccinium oxycoccos", "Vaccinium vitis-idaea", "Chamedaphne calyculata", "Larix laricina", "Betula populifolia", "Vaccinium myrtilloides", "Rhododendron groenlandicum", "Kalmia angustifolia") ~ 1,
          TRUE ~ 0
        )
    ) |>
    dplyr::select(id_sample, core_label, taxon_rank_value, taxon_organ, nitrogen_supply_growth, degree_of_decomposition, degree_of_decomposition_1, dplyr::all_of(colnames(dd_data_model_irpeat)[which(stringr::str_detect(colnames(dd_data_model_irpeat), "_in_pd")) - 1L]), hi_1630_1090, dplyr::starts_with("macrofossil_volume_fraction_"), C, N, d13C, d15N, dplyr::ends_with("_in_pd")) |>
    dplyr::mutate(
      dplyr::across(dplyr::ends_with("_in_pd"), as.logical)
    )

  res <- 
    res |> 
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), function(.x) {
        if(inherits(.x, "quantities")) {
          units::drop_units(.x)
        } else if (inherits(.x, "rvars")) {
          errors::set_errors(mean(.x), sd(.x))
        } else {
          .x
        }
      }),
      carbon_content_1 =
        dplyr::case_when(
          ! is.na(C) ~ C,
          TRUE ~ carbon_content_1
        ),
      C_to_N_1 =
        dplyr::case_when(
          ! is.na(C/N) ~ C/N,
          TRUE ~ C_to_N_1
        ),
      d13C_1 =
        dplyr::case_when(
          ! is.na(d13C) ~ d13C,
          TRUE ~ NA_real_
        ),
      d15N_1 =
        dplyr::case_when(
          ! is.na(d15N) ~ d15N,
          TRUE ~ NA_real_
        )
    ) |>
    dplyr::select(! dplyr::ends_with("_in_pd") & ! dplyr::any_of(c("eac_1", "edc_1", "macroporosity_1", "d13C", "d15N", "C", "N", "bulk_density_1", "degree_of_decomposition_1"))) |>
    tidyr::pivot_longer(
      ! dplyr::all_of(c("id_sample", "core_label", "taxon_rank_value", "taxon_organ", "nitrogen_supply_growth", "macrofossil_volume_fraction_sphagnum", "macrofossil_volume_fraction_sphagnum_hummock", "macrofossil_volume_fraction_sphagnum_lawn_hollow", "macrofossil_volume_fraction_sedge", "macrofossil_volume_fraction_shrub", "carbon_content_1", "nitrogen_content_1", "degree_of_decomposition")),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::left_join(
      res |>
        dplyr::select(id_sample, dplyr::ends_with("_in_pd"), C, N, d13C, d15N) |>
        dplyr::mutate(
          C_to_N_1_in_pd = 
            dplyr::case_when(
              ! is.na(C/N) ~ NA,
              TRUE ~ C_to_N_1_in_pd
            ),
          d13C_1_in_pd =
            dplyr::case_when(
              ! is.na(d13C) ~ NA,
              TRUE ~ d13C_1_in_pd
            ),
          d15N_1_in_pd =
            dplyr::case_when(
              ! is.na(d15N) ~ NA,
              TRUE ~ d15N_1_in_pd
            )
        ) |>
        tidyr::pivot_longer(
          cols = ! dplyr::any_of(c("id_sample", "C", "N", "d13C", "d15N")),
          names_to = "variable",
          values_to = "is_in_pd"
        ) |>
        dplyr::select(id_sample, variable, is_in_pd) |>
        dplyr::mutate(
          variable =
            variable |>
            stringr::str_remove("_in_pd$"),
          is_in_pd =
            dplyr::case_when(
              is.na(is_in_pd) ~ "measured",
              is_in_pd ~ "in prediction domain",
              ! is_in_pd ~ "not in prediction domain"
            )
        ),
      by = c("id_sample", "variable")
    ) |>
    dplyr::mutate(
      variable_or = variable, 
      variable =
        dplyr::case_when(
          variable == "bulk_density_1" ~ "BD (g cm<sup>-3</sup>)",
          variable == "C_to_N_1" ~ "C/N (g g<sup>-1</sup>)",
          variable == "H_to_C_1" ~ "H/C (g g<sup>-1</sup>)",
          variable == "O_to_C_1" ~ "O/C (g g<sup>-1</sup>)",
          variable == "d13C_1" ~ "&delta;<sup>13</sup>C (&permil;)",
          variable == "d15N_1" ~ "&delta;<sup>15</sup>N (&permil;)",
          variable == "dgf0_1" ~ "&Delta;G<sub>f</sub><sup>0</sup> (kJ mol<sub>C</sub><sup>-1</sup>)",
          variable == "nosc_1" ~ "NOSC (-)",
          variable == "hi_1630_1090" ~ "HI<sub>1630/1090</sub>",
          variable == "degree_of_decomposition_1" ~ "Degree of decomposition (%)",
          TRUE ~ variable
        ),
      variable = 
        factor(variable, levels = unique(variable)),
      plot_label =
        paste0("<i>", taxon_rank_value, "</i>"," (", taxon_organ, ifelse(is.na(nitrogen_supply_growth), "", paste0(", ", nitrogen_supply_growth)), ")") |>
        factor(levels = c("<i>Phragmites australis</i> (leaves)", paste0("<i>Phragmites australis</i> (rhizomes, ", c("low", "medium", "high"), "-N)"), "<i>Typha latifolia</i> (aboveground parts)", "<i>Sphagnum capillifolium</i> (whole plant)")),
      is_in_pd =
        dplyr::case_when(
          is.na(is_in_pd) ~ "measured",
          TRUE ~ is_in_pd
        ) |>
        factor(levels = c("measured", "in prediction domain", "not in prediction domain"))
    )
  
  res_correlation <- 
    res |>
    dplyr::group_by(variable_or) |>
    dplyr::summarise(
      correlation = cor(degree_of_decomposition, value, use = "pairwise.complete.obs"),
      .groups = "drop"
    )
  
  res_colors <- c("thistle3", "lightsalmon1", "lightsalmon2", "lightsalmon3", "lightblue", "darkseagreen")
  
  res_plot <- 
    res |>
    dplyr::left_join(
      res_correlation,
      by = "variable_or"
    ) |>
    dplyr::arrange(variable) |>
    ggplot(aes(y = value, x = degree_of_decomposition * 100)) +
    geom_errorbar(aes(ymin = value - errors::errors(value), ymax = value + errors::errors(value)), color = "lightgrey", width = 0) +
    geom_smooth(aes(color = plot_label, group = plot_label, fill = plot_label), method = "lm") +
    geom_point(aes(color = plot_label, shape = is_in_pd)) +
    scale_color_manual(values = res_colors) +
    scale_fill_manual(values = res_colors) +
    facet_wrap(~ paste0(variable, ifelse(variable_or == "degree_of_decomposition_1", "<br>", ""), " &rho; = ", round(correlation, digits = 2L)) |> factor(levels = unique(paste0(variable, ifelse(variable_or == "degree_of_decomposition_1", "<br>", ""), " &rho; = ", round(correlation, digits = 2L)))), scales = "free_y", ncol = 4L) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      strip.text.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.text = ggtext::element_markdown(),
      legend.box = "horizontal",
      legend.direction = "vertical"
    ) +
    labs(
      y = "Decomposition indicator value",
      x = "Measured degree of decomposition (%)"
    ) +
    guides(
      color = guide_legend(title = "Species (organ, initial N content)", override.aes = list(size = 3), ncol = 2L),
      fill = guide_legend(title = "Species (organ, initial N content)"),
      shape = guide_legend(title = "Reliability", override.aes = list(size = 3))
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8.5, height = 5.2, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots results of `dd_simulation_1`
#' 
#' @export
dd_make_plot_13 <- function(dd_simulation_1, file_plot = "figures/dd_plot_13.pdf") {
  
  res_plot <- 
    dd_simulation_1 |>
    dplyr::filter(! is.na(m_0_classical)) |>
    ggplot(aes(y = ratio_m0, x = m_1_1/m_2_1)) +
    geom_path(aes(color = gamma_1 - gamma_2, group = paste0(gamma_1, "_", gamma_2))) +
    facet_grid( ~ gamma_1) +
    scale_color_gradient2(mid = "lightgrey") +
    scale_x_log10(labels = function(x) round(x, 1)) +
    guides(color = guide_colorbar(title = "<i>&gamma;<sub>1</sub>(t) - &gamma;<sub>2</sub>(t)</i>")) +
    labs(
      y = "<i>m(t<sub>0</sub>)</i><sub>naive</sub> / <i>m(t<sub>0</sub>)</i> (g g<sup>-1</sup>)",
      x = "log(<i>m<sub>1</sub>(t) / <i>m<sub>2</sub>(t))"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown()
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8.5, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots the prediction domains of the models along with spectra of the peat samples
#' 
#' @param dd_data_pmird_peat_cores_yhat for model 3.
#' 
#' @export
dd_make_plot_14 <- function(dd_data_pmird_peat_cores, dd_mir_preprocessing_config, dd_model_info, dd_data_pmird_peat_cores_yhat, dd_stan_1_fit, file_plot) {
  
  prediction_domain <- dd_model_info$prediction_domain
  
  # preprocess spectra
  res <- 
    purrr::map(seq_along(dd_mir_preprocessing_config), function(i) {
      irp_preprocess_eb1149(
        x = dd_data_pmird_peat_cores,
        config = dd_mir_preprocessing_config[[i]]
      ) |>
        dplyr::mutate(
          model_name = dd_model_names(i, dd_model_info)
        )
    })
  
  # model coefficients
  res_coefficients <- 
    purrr::map(seq_along(dd_stan_1_fit), function(i) {
      .x <- dd_stan_1_fit[[i]]
      parnames <- brms::variables(.x)
      .x |>
        tidybayes::gather_rvars(!!!syms(parnames[stringr::str_detect(parnames, "^b_x")])) |>
        dplyr::rename(
          variable = ".variable",
          value = ".value"
        ) |>
        dplyr::mutate(
          x = res[[1]]$spectra[[1]]$x,
          value = median(value),
          model_name = dd_model_names(i, dd_model_info)
        )
    }) |>
    dplyr::bind_rows()
  
  index <- 
    res[[1]] |>
    dplyr::mutate(
      id_sample = seq_along(spectra),
      degree_of_decomposition = median(dd_data_pmird_peat_cores_yhat$degree_of_decomposition)
    ) |>
    dplyr::arrange(degree_of_decomposition) |>
    dplyr::slice(seq(1, length(spectra), length.out = 5) |> as.integer())
  
  res_spectra <- 
    purrr::map(seq_along(dd_mir_preprocessing_config), function(i) {
      res[[i]] |>
        dplyr::slice(index$id_sample[seq_along(index$id_sample)]) |> #---note: strange construct to avoid vctrs error
        dplyr::select(model_name, id_sample, id_measurement, spectra) |>
        dplyr::mutate(
          intensity_3400 = index$degree_of_decomposition
        )
    }) |>
    dplyr::bind_rows()
  
  res_pd <- 
    purrr::map(seq_along(dd_mir_preprocessing_config), function(i) {
      res[[i]] |>
        irpeat::irp_as_irp_prediction_domain() |>
        dplyr::mutate(
          plot_label = "Peat data"
        ) |>
        dplyr::bind_rows(
          prediction_domain[[i]] |>
            dplyr::mutate(
              plot_label = "Training data"
            )
        ) |>
        dplyr::mutate(
          model_name = dd_model_names(i, dd_model_info)
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::left_join(
      res_coefficients |>
        dplyr::select(model_name, x, value),
      by = c("model_name", "x")
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("ymin", "ymax")),
        function(.x) .x * value
      )
    )
  
  res_plot <- 
    res_pd |>
    dplyr::mutate(
      plot_label =
        factor(plot_label, levels = rev(unique(plot_label)))
    ) |>
    as.data.frame() |>
    ggplot() +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, group = paste0(model_name, "_", plot_label), fill = plot_label), color = NA, alpha = 0.6) +
    geom_path(
      data = dd_data_pmird_peat_cores$spectra[[1]],
      aes(y = y, x = x)
    ) +
    geom_path(
      data = 
        res_spectra |>
        tidyr::unnest("spectra") |>
        dplyr::left_join(
          res_coefficients |>
            dplyr::select(model_name, x, value),
          by = c("model_name", "x")
        ) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::all_of(c("y")),
            function(.x) .x * value
          )
        ),
      aes(y = y, x = x, color = intensity_3400, group = id_sample)
    ) +
    facet_wrap(~ model_name, scales = "free_y", ncol = 1L) +
    scale_fill_manual(values = c("salmon", "steelblue")) +
    labs(
      y = "Intensity (-)",
      x = "Wavenumber (cm<sup>-1</sup>)"
    ) +
    scale_x_continuous(minor_breaks = seq(600 , 4000, 200)) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      panel.grid.major.x = element_line(colour = "grey", size = 0.5),
      panel.grid.minor.x = element_line(colour = "grey", size = 0.5)
    ) +
    geom_vline(xintercept = c(670, 3570, 3750))
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 10, height = 6, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}



#' Plots measured values versus predictions for samples held out during CV and for the training data
#' 
#' @param rmse_df A list with objects like `dd_stan_1_fit_rmse`.
#' 
#' @export
dd_make_plot_15 <- function(dd_data_model_evaluation_1, rmse_df, dd_model_info, file_plot = "figures/dd_plot_15.pdf") {
  
  rmse_df <- 
    purrr::map(seq_along(rmse_df), function(i) {
      readRDS_rvars(rmse_df[[i]]) |>
        dplyr::mutate(
          rmse = rmse * 100,
          rmse_mean = mean(rmse),
          rmse_lower = posterior::quantile2(rmse, probs = 0.025),
          rmse_upper = posterior::quantile2(rmse, probs = 0.975),
          label = paste0("RMSE = ", round(rmse_mean, 1), " (", round(rmse_lower, 1), ",", round(rmse_upper, 1), ")"),
          cv_type =
            switch(
              i,
              "1" = "Fitted values",
              "2" = "Stratified CV",
              "3" = "Grouped CV"
            )
        )
    }) |>
    dplyr::bind_rows() |>
    dplyr::select(-rmse)
  
  res_colors <- c("thistle3", "lightsalmon1", "lightsalmon2", "lightsalmon3", "lightblue", "darkseagreen", "grey")
  
  cv_type_levels <- c("Fitted values", "Stratified CV", "Grouped CV")
  
  res_plot <- 
    dd_data_model_evaluation_1 |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("y")& dplyr::where(is.numeric),
        function(.x) .x * 100
      ),
      cv_type =
        factor(cv_type, levels = cv_type_levels),
      id_model =
        dd_model_names(id_model, dd_model_info = dd_model_info),
      nitrogen_supply_growth =
        dplyr::case_when(
          dd_dataset_label != "reuter2019" ~ NA_character_,
          dd_dataset_label == "reuter2019" & taxon_organ == "leaves" ~ NA_character_,
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^K") ~ "high-N",
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^T") ~ "low-N",
          dd_dataset_label == "reuter2019" & stringr::str_detect(sample_label, "^S") ~ "medium-N"
        ),
      taxon_rank_value =
        dplyr::case_when(
          dd_dataset_label == "pmird_undecomposed_litter" ~ "Other",
          TRUE ~ taxon_rank_value
        ),
      plot_label =
        dplyr::case_when(
          dd_dataset_label == "pmird_undecomposed_litter" ~ "Other",
          TRUE ~ paste0("<i>", taxon_rank_value, "</i>"," (", taxon_organ, ifelse(is.na(nitrogen_supply_growth), "", paste0(", ", nitrogen_supply_growth)), ")")
        ) |>
        factor(levels = c("<i>Phragmites australis</i> (leaves)", paste0("<i>Phragmites australis</i> (rhizomes, ", c("low", "medium", "high"), "-N)"), "<i>Typha latifolia</i> (aboveground parts)", "<i>Sphagnum capillifolium</i> (whole plant)", "Other")),
    ) |>
    ggplot(aes(y = y, x = yhat_mean)) +
    geom_errorbarh(aes(xmin = yhat_lower, xmax = yhat_upper), color = "lightgrey", height = 0) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    #geom_smooth(aes(color = plot_label, fill = plot_label), method = "lm") +
    geom_point(aes(color = plot_label)) +
    geom_text(
      data = 
        rmse_df |>
        dplyr::mutate(
          cv_type = factor(cv_type, levels = cv_type_levels)
        ), 
      aes(x = 5, y = 90, label = label),
      size = 2, hjust = "left"
    ) +
    scale_color_manual(values = res_colors) +
    scale_fill_manual(values = res_colors) +
    facet_grid(id_model ~ cv_type) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      strip.text.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      legend.text = ggtext::element_markdown(),
      legend.box = "horizontal",
      legend.direction = "vertical"
    ) +
    coord_fixed() +
    labs(
      x = "Fitted or predicted degree of decomposition (%)",
      y = "Measured degree of decomposition (%)"
    ) +
    guides(
      color = guide_legend(title = "Species (organ, initial N content)", override.aes = list(size = 3), ncol = 2L),
      fill = guide_legend(title = "Species (organ, initial N content)")
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6, height = 6, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots the spectrum for `dd_data_pmird_mineral`.
#' 
#' @export
dd_make_plot_16 <- function(dd_data_pmird_mineral, file_plot = "figures/dd_plot_16.pdf") {
  
  res_plot <- 
    dd_data_pmird_mineral |>
    dd_preprocess_mir_for_plotting() |>
    plot() +
    labs(
      y = "Intensity (-)",
      x = "Wavenumber (cm<sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown()
    )
 
  ggsave(
    file_plot,
    plot = res_plot,
    width = 6, height = 3.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
   
}


#' Plot of `dd_data_pmird_peat_cores` spectra within and outside the prediction domain for the different models
#' 
#' @export
dd_make_plot_17 <- function(dd_data_pmird_peat_cores, dd_data_pmird_peat_cores_yhat, file_plot = "figures/dd_plot_17.pdf") {
  
  # import information on pd
  res_pd <- 
    purrr::map(seq_along(dd_data_pmird_peat_cores_yhat), function(i) {
      dd_data_pmird_peat_cores_yhat[[i]] |>
      readRDS_rvars() |>
        dplyr::mutate(
          id_model = i
        )
    }) |>
    dd_do_call("rbind")
  
  # preprocess spectra
  res <- 
    dd_data_pmird_peat_cores |>
    dd_preprocess_mir_for_plotting() |>
    dplyr::select(dd_dataset_label, id_dataset, id_sample, spectra)
  
  res <- 
    dplyr::left_join(
      res_pd, 
      res,
      by = "id_sample"
    ) |>
    ir::ir_as_ir() |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.logical), as.logical)
    )
  
  res_plot <- 
    res |>
    ir::ir_bin(width = 20) |>
    dplyr::mutate(
      degree_of_decomposition_in_pd =
        dplyr::case_when(
          degree_of_decomposition_in_pd ~ "Within prediction domain",
          TRUE ~ "Outside prediction domain"
        ) |>
        factor(levels = c("Within prediction domain", "Outside prediction domain"))
    ) |>
    tidyr::unnest("spectra") |>
    ggplot(aes(y = y, x = x, group = id_sample)) +
    geom_path() +
    facet_grid(id_model ~ degree_of_decomposition_in_pd) +
    labs(
      y = "Intensity (-)",
      x = "Wavenumber (cm<sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      strip.background.y = element_blank()
    )
  
  res |>
    ggplot(aes(y = C, x = degree_of_decomposition_in_pd, color = id_model)) +
    geom_boxplot()
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 8, height = 6, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plot of C content, HI_1630/1090 for peat samples within and outside the prediction domain for the different models
#' 
#' @export
dd_make_plot_18 <- function(dd_data_pmird_peat_cores_irpeat, dd_data_pmird_peat_cores_yhat, file_plot = "figures/dd_plot_18.pdf") {
  
  # import information on pd
  res_pd <- 
    purrr::map(seq_along(dd_data_pmird_peat_cores_yhat), function(i) {
      dd_data_pmird_peat_cores_yhat[[i]] |>
        readRDS_rvars() |>
        dplyr::mutate(
          id_model = i
        )
    }) |>
    dd_do_call("rbind")
  
  # preprocess spectra
  res <- 
    dd_data_pmird_peat_cores_irpeat
  
  res <- 
    dplyr::left_join(
      res_pd, 
      res |>
        ir::ir_drop_spectra(),
      by = "id_sample"
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.logical), as.logical)
    )
  
  res_plot <- 
    res |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), as.numeric),
      degree_of_decomposition_in_pd =
        dplyr::case_when(
          degree_of_decomposition_in_pd ~ "Within prediction domain",
          TRUE ~ "Outside prediction domain"
        ) |>
        factor(levels = c("Within prediction domain", "Outside prediction domain")),
      carbon_content_1 =
        dplyr::case_when(
          ! is.na(C) ~ C,
          TRUE ~ carbon_content_1
        ),
      hi_1630_1090 = hi_1630_1090
    ) |>
    dplyr::select(id_model, carbon_content_1, hi_1630_1090, degree_of_decomposition_in_pd) |>
    tidyr::pivot_longer(
      dplyr::all_of(c("carbon_content_1", "hi_1630_1090")),
      names_to = "variable",
      values_to ="value"
    ) |>
    dplyr::mutate(
      variable =
        dplyr::case_when(
          variable == "carbon_content_1" ~ "C",
          variable == "hi_1630_1090" ~ "HI<sub>1630/1090</sub>"
        )
    ) |>
    ggplot(aes(x = value)) +
    geom_density(aes(color = degree_of_decomposition_in_pd, fill = degree_of_decomposition_in_pd), alpha = 0.6) +
    facet_grid(id_model ~ variable, scales = "free") +
    labs(
      y = "Density (g g<sup>-1</sup>)",
      x = "C content (g g<sup>-1</sup>) or HI<sub>1630/1090</sub> (-)"
    ) +
    scale_color_manual(values = c("black", "grey")) +
    scale_fill_manual(values = c("black", "grey")) +
    guides(
      color = guide_legend(title = ""),
      fill = guide_legend(title = "")
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(),
      strip.text.y = ggtext::element_markdown()
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 4, height = 5.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots results of `dd_simulation_1`, but now for the conditions where we think NPP may be reconstructed
#' 
#' @export
dd_make_plot_19 <- function(dd_simulation_1, mass_fraction_threshold = 0.9, error_threshold = 0.2, file_plot = "figures/dd_plot_19.pdf") {
  
  res_plot <- 
    dd_simulation_1 |>
    dplyr::mutate(
      m_0_naive = m_1 * (m_1_1/m_1) / (1 - gamma_mirs), #---note: (m_1_1/m_1): assumes we know the mass of component 1
      ratio_m0 = m_0_naive / m_1_0
    ) |>
    dplyr::filter(! is.na(m_0_naive) & m_1_1/m_1 >= mass_fraction_threshold) |>
    ggplot(aes(y = ratio_m0, x = m_1_1/m_1)) +
    geom_path(aes(color = gamma_1 - gamma_2, group = paste0(gamma_1, "_", gamma_2))) +
    facet_grid( ~ gamma_1) +
    scale_color_gradient2(mid = "lightgrey") +
    scale_x_continuous(breaks = c(0.9, 0.95, 1.0)) +
    #scale_x_log10(labels = function(x) round(x, 1), breaks = c(10, 100, 900)) +
    guides(color = guide_colorbar(title = "<i>&gamma;<sub>1</sub>(t) - &gamma;<sub>2</sub>(t)</i>")) +
    labs(
      y = "<i>m(t<sub>0</sub>)</i><sub>naive</sub> / <i>m<sub>1</sub>(t<sub>0</sub>)</i> (g g<sup>-1</sup>)",
      x = "log(<i>m<sub>1</sub>(t) / <i>m<sub>2</sub>(t))"
    ) +
    geom_hline(yintercept = c(1.0), color = "grey50") +
    geom_hline(yintercept = c(1 - error_threshold, 1 + error_threshold), color = "grey50", linetype = 2) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      legend.title = ggtext::element_markdown(),
      panel.spacing = unit(0.2, "in")
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 10, height = 5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  list(
    file_plot = file_plot,
    mass_fraction_threshold = mass_fraction_threshold,
    error_threshold = error_threshold
  )
  
}



#' Plots reconstructed NPP versus layer age for different peat cores and models predicting the degree of decomposition used to reconstruct the NPP
#' 
#' @export
dd_make_plot_20 <- function(dd_reconstruction_initial_mass_rate_1, file_plot = "figures/dd_plot_20.pdf") {
  
  res <- 
    dd_reconstruction_initial_mass_rate_1 |>
    dplyr::mutate(
      layer_amar = (as.numeric(units::set_units(bulk_density_1, g/(m^2 * cm))) * (sample_depth_lower - sample_depth_upper)) / (age_lower - age_upper), #---todo: replace by measured bulk density
      core_label =
        dplyr::case_when(
          core_label == "eb1005_MH1" ~ "MH1",
          core_label == "eb1006_MK1" ~ "MK1",
          core_label == "eb1018_OD2" ~ "OD2"
        ) |>
        factor(levels = c("MH1", "MK1", "OD2")),
      taxon_rank_value =
        paste0("<i>", taxon_rank_value, "</i>")
    )
  
  res_plot <- 
    res |>
    ggplot(aes(x = mean((age_upper + age_lower)/2), color = taxon_rank_value)) +
    ggdist::stat_pointinterval(
      aes(ydist = layer_amar/1000), 
      .width = c(0.5, 0.9),
      point_size = 2
    ) +
    #geom_smooth(method = "lm") +
    facet_grid(id_model ~ core_label, scales = "free_x") +
    #coord_cartesian(ylim = c(0, 1.5)) +
    labs(
      y = "Estimated NPP (kg m<sup>-2</sup> yr<sup>-1</sup>)",
      x = "Average age of layer middle (yr)"
    ) +
    scale_color_manual(values = c("darkseagreen", "lightsalmon3", "lightsalmon1", "grey70")) +
    #scale_fill_manual(values = c("black", "grey")) +
    guides(
      color = guide_legend(title = "", override.aes = list(size = 3)),
      #fill = guide_legend(title = "")
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 14),
      strip.text.y = ggtext::element_markdown(size = 14),
      legend.text = ggtext::element_markdown(size = 12)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 4.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Dummy plot function for rbacon models
#' 
#' @export
dd_make_plot_21 <- function(dd_reconstruction_initial_mass_rate_1) {
  
  target_core_label <- 
    dd_reconstruction_initial_mass_rate_1 |>
    dplyr::pull(core_label) |>
    unique()
  
  list.files(paste0("data/derived_data/Bacon_runs/", target_core_label), pattern = "\\.pdf$", full.names = TRUE)
  
}


#' Plot for simulation 2
#' 
#' @export
dd_make_plot_22 <- function(dd_simulation_2, file_plot = "figures/dd_plot_22.pdf") {
  
  res_plot <- 
    dd_simulation_2 |>
    dplyr::filter(! is.na(bias)) |>
    dplyr::arrange(gamma_1, gamma_2, m_1_1/m_1) |>
    ggplot(aes(y = bias, x = m_1_1/m_1)) +
    geom_path(aes(color = gamma_1 - gamma_2, group = paste0(gamma_1, "_", gamma_2))) +
    facet_grid( ~ gamma_1) +
    scale_color_gradient2(
      mid = "lightgrey", 
      low = scales::muted("blue"), 
      high = scales::muted("red")
    ) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    #scale_x_log10(labels = function(x) round(x, 1)) +
    guides(color = guide_colorbar(title = "<i>&gamma;<sub>1</sub>(t) - &gamma;<sub>2</sub>(t)</i>")) +
    labs(
      y = "<i>&gamma;</i><sub>MIRS</sub>(<i>t</i>) - <i>&gamma;(t)</i> (g g<sup>-1</sup>)",
      x = "<i>m</i><sub>1</sub>(t) / (<i>m</i><sub>1</sub>(t) + <i>m</i><sub>2</sub>(t)) (g g<sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(), 
      panel.spacing = unit(1.2, "lines"),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      legend.title = ggtext::element_markdown(size = 14),
      strip.text.x = ggtext::element_markdown(size = 14)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 4, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plots for simulation 3
#' 
#' @export
dd_make_plot_23 <- function(dd_simulation_3, dd_simulation_3_species, variable_degree_of_decomposition, file_plot = "figures/dd_plot_23.pdf") {
  
  res_plot <- 
    dd_simulation_3  |>
    dplyr::mutate(
      shape_variable = 
        ifelse(gamma_1 > min(gamma_1[gamma_1 > 0]), ">0", "0") |>
        factor(levels = c("0", ">0"))
    ) |>
    #dplyr::filter(taxon_rank_value %in% c("Aulacomnium palustre", "Betula populifolia", "Calluna vulgaris", "Eriophorum vaginatum", "Phragmites australis", "Sphagnum capillifolium", "Sphagnum fallax", "Typha latifolia")) |>
    dplyr::filter(taxon_rank_value %in% dd_simulation_3_species) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("taxon_rank_value"),
        function(.x) {
          stringr::str_replace(.x, pattern = "[A-Z]{1}.+ ", replacement = paste0(stringr::str_sub(.x, start = 1L, end = 1L), ". "))
        }
      )
    ) |>
    ggplot(aes(y = bias_hat, x = as.numeric(eval(variable_degree_of_decomposition)) - gamma)) +
    geom_point(aes(color = gamma_1 - gamma_2, shape = shape_variable), size = 2) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_gradient2(
      mid = "lightgrey", 
      low = scales::muted("blue"), 
      high = scales::muted("red")
    ) +
    guides(
      color = guide_colorbar(title = "<i>&gamma;<sub>1</sub> - &gamma;<sub>2</sub></i>"),
      shape = guide_legend(title = "<i>&gamma;</i><sub>1</sub> (g g<sup>-1</sup>)", override.aes = list(size = 3))
    ) +
    labs(
      y = "<i>&gamma;</i><sub>MIRS</sub> - <i>&gamma;</i> (g g<sup>-1</sup>) with theoretical <i>&gamma;</i><sub>MIRS</sub>",
      x = "<i>&gamma;</i><sub>MIRS</sub> - <i>&gamma;</i> (g g<sup>-1</sup>) with predicted <i>&gamma;</i><sub>MIRS</sub>"
    ) +
    facet_grid(
      paste0("<i>", as.character(taxon_rank_value), "</i>") ~ paste0("<i>", as.character(taxon_rank_value2), "</i>")
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 12),
      axis.title.y = ggtext::element_markdown(size = 12),
      legend.title = ggtext::element_markdown(size = 12),
      #legend.text = ggtext::element_markdown(size = 10),
      strip.text.x = ggtext::element_markdown(size = 10),
      strip.text.y = ggtext::element_markdown(size = 10)
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 5.5, height = 7, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}


#' Plot for simulation 2: here: using the theoretical degree of decomposition measured on a bulk sample as estimate for the degree of decomposition of a dominant component
#' 
#' @export
dd_make_plot_24 <- function(dd_simulation_2, file_plot = "figures/dd_plot_24.pdf") {
  
  res_plot <- 
    dd_simulation_2 |>
    dplyr::filter(! is.na(bias)) |>
    dplyr::mutate(
      bias = abs(gamma_1 - gamma_mirs)
    ) |>
    dplyr::filter(bias <= 0.2) |>
    dplyr::arrange(gamma_1, gamma_2, m_1_1/m_1) |>
    ggplot(aes(y = gamma_1 - gamma_mirs, x = m_1_1/m_1)) +
    geom_path(aes(color = gamma_1 - gamma_2, group = paste0(gamma_1, "_", gamma_2))) +
    facet_grid( ~ gamma_1) +
    scale_color_gradient2(
      mid = "lightgrey", 
      low = scales::muted("blue"), 
      high = scales::muted("red")
    ) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    #scale_x_log10(labels = function(x) round(x, 1)) +
    guides(color = guide_colorbar(title = "<i>&gamma;<sub>1</sub>(t) - &gamma;<sub>2</sub>(t)</i>")) +
    labs(
      y = "<i>&gamma;</i><sub>MIRS</sub>(<i>t</i>) - <i>&gamma;<sub>1</sub>(t)</i> (g g<sup>-1</sup>)",
      x = "<i>m</i><sub>1</sub>(t) / (<i>m</i><sub>1</sub>(t) + <i>m</i><sub>2</sub>(t)) (g g<sup>-1</sup>)"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(), 
      panel.spacing = unit(1.2, "lines"),
      legend.position = "bottom",
      axis.title.x = ggtext::element_markdown(size = 14),
      axis.title.y = ggtext::element_markdown(size = 14),
      legend.title = ggtext::element_markdown(size = 14),
      strip.text.x = ggtext::element_markdown(size = 14)
    ) +
    geom_vline(xintercept = 0.9, color = "grey50", linetype = 2)
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 12, height = 4, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}

#' Plot for simulation 4
#' 
#' @export
dd_make_plot_25 <- function(dd_simulation_4, file_plot = "figures/dd_plot_25.pdf") {
  
  res_plot <- 
    dd_simulation_4 |>
    dplyr::arrange(k_0_1, alpha, k_0_2, decomposition_progress) |>
    ggplot(aes(y = gamma_difference, x = decomposition_progress, group = k_0_2)) +
    geom_path(aes(color = k_0_2)) +
    guides(
      color = guide_colorbar(title = "<i>k</i><sub>0,2</sub> (yr<sup>-1</sup>)")
    ) +
    facet_grid(k_0_1 ~ round(alpha, 0)) +
    labs(
      y = "&gamma;<sub>1</sub> - &gamma;<sub>s</sub> (g g<sup>-1</sup>)",
      x = "Decomposition progress"
    ) +
    scale_x_log10() +
    geom_hline(yintercept = c(-0.5, 0.5), color = "grey") +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      axis.title.x = ggtext::element_markdown(size = 13),
      axis.title.y = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 12),
      legend.position = "bottom",
      panel.spacing = unit(0.2, "in")
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 9, height = 7.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}

#' Plot for simulation 4: Relation between the degree of decomposition and the decomposition progress
#' 
#' @export
dd_make_plot_26 <- function(dd_simulation_4, file_plot = "figures/dd_plot_26.pdf") {
  
  res_plot <- 
    dd_simulation_4 |>
    dplyr::arrange(k_0_1, alpha, k_0_2, decomposition_progress) |>
    ggplot(aes(x = decomposition_progress, y = (gamma_1 + gamma_2)/2, group = k_0_2)) +
    geom_path(aes(color = k_0_2)) +
    guides(
      color = guide_colorbar(title = "<i>k</i><sub>0,2</sub> (yr<sup>-1</sup>)")
    ) +
    scale_x_log10() +
    facet_grid(k_0_1 ~ round(alpha, 0)) +
    labs(
      y = "&gamma; (g g<sup>-1</sup>)",
      x = "Decomposition progress"
    ) +
    theme_classic() +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = ggtext::element_markdown(size = 13),
      strip.text.y = ggtext::element_markdown(size = 13),
      axis.title.x = ggtext::element_markdown(size = 13),
      axis.title.y = ggtext::element_markdown(size = 13),
      legend.title = ggtext::element_markdown(size = 12),
      legend.position = "bottom",
      panel.spacing = unit(0.2, "in")
    )
  
  ggsave(
    file_plot,
    plot = res_plot,
    width = 9, height = 7.5, 
    dpi = 300,
    device = cairo_pdf
  )
  
  file_plot
  
}
