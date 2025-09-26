# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(quantities)
library(ir)
library(crew)
library(ggplot2)
library(posterior)
library(future)

crew_sequential <-
  crew::crew_controller_local(
    name = "crew_sequential",
    workers = 1L,
    tasks_max = 1L
  )

crew_parallel_1 <-
  crew::crew_controller_local(
    name = "crew_parallel_1",
    workers = 4L
  )

crew_parallel_2 <-
  crew::crew_controller_local(
    name = "crew_parallel_2",
    workers = 5L
  )

# Set target options:
tar_option_set(
  packages = c("tibble", "ir", "ggplot2", "posterior", "quantities", "errors", "brms", "rstan"),
  format = "rds",
  controller = crew_controller_group(crew_sequential, crew_parallel_1, crew_parallel_2)
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# to avoid deletion of aux files
options(tinytex.clean = FALSE)

# Replace the target list below with your own:
list(
  #### data import ####
  tar_target(
    dd_data_model,
    command = dd_make_data_model()
  ),
  tar_target(
    dd_data_pmird_peat_cores,
    dd_make_data_pmird_peat_cores()
  ),
  tar_target(
    dd_data_pmird_peat_cores_macrofossils,
    dd_make_data_pmird_peat_cores_macrofossils()
  ),
  tar_target(
    dd_data_pmird_peat_cores_accumulation,
    dd_make_data_pmird_peat_cores_accumulation()
  ),
  tar_target(
    dd_data_pmird_mineral,
    dd_make_data_pmird_mineral()
  ),
  tar_target(
    dd_data_pmird_decomposed,
    command = 
      dd_data_pmird_peat_cores |>
      dplyr::filter(id_dataset == 8 & sample_depth_lower == 240)
  ),
  tar_target(
    dd_data_pmird_peat_cores_irpeat,
    dd_make_data_pmird_add_irpeat(x = dd_data_pmird_peat_cores)
  ),
  tar_target(
    dd_data_pmird_peat_cores_14C,
    command = dd_make_data_pmird_peat_cores_14C()
  ),
  tar_target(
    dd_data_model_irpeat,
    dd_make_data_pmird_add_irpeat(x = dd_data_model)
  ),
  # tar_target(
  #   dd_data_age_depth_model,
  #   command = 
  #     dd_make_data_age_depth_model()
  # ),
  tar_target(
    dd_data_ta_wtd,
    command = 
      dd_make_data_ta_wtd()
  ),
  # tar_target(
  #  dd_mir_co2_h20_background_file,
  #  command = "data/derived_data/bg_spectrum_pure.rds",
  #  format = "file"
  #),
  #tar_target(
  #  dd_mir_co2_h20_background,
  #  command = readRDS("data/derived_data/bg_spectrum_pure.rds")
  #),
  tar_target(
    dd_data_bona2018_moss_npp,
    command = dd_make_data_bona2018_moss_npp()
  ),
  tar_target(
    dd_data_bengtsson2021_moss_npp,
    command = dd_make_data_bengtsson2021_moss_npp()
  ),
  #### simulations ####
  tar_target(
    dd_simulation_2,
    command = dd_make_simulation_2()
  ),
  tar_target(
    dd_simulation_3,
    command = dd_make_simulation_3(dd_data_model = dd_data_model)
  ),
  tar_target(
    dd_simulation_3_species,
    command = c("Aulacomnium palustre", "Eriophorum vaginatum", "Phragmites australis", "Sphagnum capillifolium", "Typha latifolia")
  ),
  tar_target(
    dd_simulation_4,
    command = dd_make_simulation_4()
  ),
  tar_target(
    dd_simulation_5,
    command = 
      dd_make_simulation_5(
        dd_data_model = dd_data_model
      )
  ),
  tar_target(
    dd_simulation_6,
    command = 
      dd_make_simulation_6(
        dd_data_model = dd_data_model
      )
  ),
  #### age-depth models (recomputed) ####
  tar_target(
    dd_rbacon_csv,
    command = dd_prepare_rbacon_csv(dd_data_pmird_peat_cores_14C = dd_data_pmird_peat_cores_14C)
  ),
  tar_target(
    dd_rbacon,
    pattern = map(dd_rbacon_csv),
    command = 
      dd_run_rbacon(
        dd_rbacon_csv = dd_rbacon_csv, 
        hiatus.depths = NA, 
        postbomb = 1,
        dd_stan_1_mcmc_settings = dd_stan_1_mcmc_settings
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  #### spectral preprocessing ####
  tar_target(
    dd_mir_preprocessing_settings,
    dd_get_mir_preprocessing_settings()
  ),
  tar_target(
    dd_data_model_preprocessed,
    dd_make_data_model_preprocessed(
      dd_data_model = dd_data_model,
      dd_mir_preprocessing_settings = dd_mir_preprocessing_settings
    )
  ),
  #### models ####
  tar_target(
    dd_model_info,
    dd_get_model_info(
      dd_data_model_preprocessed = dd_data_model_preprocessed, 
      dd_mir_preprocessing_settings = dd_mir_preprocessing_settings
    )
  ),
  tar_target(
    dd_id_model,
    command = dd_model_info$id_model
  ),
  tar_target(
    dd_brms_compiled_model,
    dd_make_brms_compiled_model(
      dd_model_info[1, ],
      backend = "cmdstanr",
      sig_figs = 12L
    )
  ),
  tar_target(
    dd_stan_1_mcmc_settings,
    dd_get_stan_1_mcmc_settings()
  ),
  tar_target(
    dd_stan_1_fit,
    pattern = map(dd_id_model),
    command = 
      update(
        dd_brms_compiled_model, 
        newdata = 
          dd_model_info |> 
          dplyr::filter(id_model == dd_id_model) |>
          dd_make_brms_data(),
        iter = dd_stan_1_mcmc_settings$iter,
        warmup = dd_stan_1_mcmc_settings$warmup,
        chains = dd_stan_1_mcmc_settings$chains,
        cores = dd_stan_1_mcmc_settings$chains,
        control = 
          list(
            max_treedepth = dd_stan_1_mcmc_settings$max_treedepth,
            adapt_delta = dd_stan_1_mcmc_settings$adapt_delta
          ),
        save_pars = brms::save_pars(all = TRUE),
        backend = "cmdstanr",
        save_warmup = TRUE,
        sig_figs = 14L
      ) |>
      list()
  ), 
  #### grouped cross-validation ####
  tar_target(
    dd_stan_1_fit_rmse,
    command = 
      dd_make_fit_rmse(
        fit = dd_stan_1_fit,
        dd_model_info = dd_model_info,
        file = "targets_rvars/dd_stan_1_fit_rmse.rds"
      ),
    format = "file"
  ),
  tar_target(
    dd_stan_1_kfold_cv_folds,
    command = dd_make_cv_groups(dd_data_model, K = 10L)
  ),
  tar_target(
    dd_stan_1_kfold,
    pattern = map(dd_id_model),
    command = 
      {
        options(future.globals.maxSize = 1300 * 1024^2)
        future::plan("sequential")
        res <- 
          brms::kfold(
            x = dd_stan_1_fit[[dd_id_model]],
            chains = 4L,
            cores = 1L,
            folds = dd_stan_1_kfold_cv_folds,
            joint = "obs",
            save_fits = TRUE,
            recompile = FALSE
          ) |>
          list()
        options(future.globals.maxSize = 500 * 1024^2)
        res
      },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ), 
  tar_target(
    dd_stan_1_kfold_comparison,
    command = brms::loo_compare(dd_stan_1_kfold)
  ),
  tar_target(
    dd_stan_1_kfold_rmse,
    command = 
      dd_make_kfold_rmse(
        kfold = dd_stan_1_kfold,
        dd_model_info = dd_model_info,
        file = "targets_rvars/dd_stan_1_kfold_rmse.rds"
      ),
    format = "file"
  ),
  #### stratified cross-validation ####
  tar_target(
    dd_stan_1_kfold_cv_folds_2,
    command = dd_make_cv_groups_2(dd_data_model, K = 10L)
  ),
  tar_target(
    dd_stan_1_kfold_2,
    pattern = map(dd_id_model),
    command = 
      {
        options(future.globals.maxSize = 1300 * 1024^2)
        future::plan("sequential")
        res <- 
          brms::kfold(
            x = dd_stan_1_fit[[dd_id_model]],
            chains = 4L,
            cores = 1L,
            folds = dd_stan_1_kfold_cv_folds_2,
            joint = "obs",
            save_fits = TRUE,
            recompile = FALSE
          ) |>
          list()
        options(future.globals.maxSize = 500 * 1024^2)
        res
      },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    dd_stan_1_kfold_2_rmse,
    command = 
      dd_make_kfold_rmse(
        kfold = dd_stan_1_kfold_2,
        dd_model_info = dd_model_info,
        file = "targets_rvars/dd_stan_1_kfold_2_rmse.rds"
      ),
    format = "file"
  ),
  tar_target(
    dd_stan_1_kfold_2_comparison,
    command = brms::loo_compare(dd_stan_1_kfold_2)
  ),
  #### model interface for prediction ####
  tar_target(
    dd_mir_preprocessing_config,
    command = 
      dd_make_mir_preprocessing_config(
        dd_mir_preprocessing_settings = dd_mir_preprocessing_settings, 
        dd_model_info = dd_model_info
      )
  ),
  # for irpeatmodels and irpeat
  tar_target(
    dd_output_1,
    command =
      dd_make_output_1(
        dd_model_info = dd_model_info,
        dd_stan_1_fit = dd_stan_1_fit,
        dd_mir_preprocessing_config = dd_mir_preprocessing_config
      ),
    format = "file"
  ),
  #### model evaluation ####
  tar_target(
    dd_stan_1_fit_mcse,
    command = 
      dd_get_mcse_1(
        dd_stan_1_fit = dd_stan_1_fit,
        dd_mir_preprocessing_config = dd_mir_preprocessing_config,
        dd_data_model = dd_data_model, 
        dd_data_pmird_peat_cores = dd_data_pmird_peat_cores, 
        dd_model_info = dd_model_info
      )
  ),
  tar_target(
    dd_stan_1_fit_rhat,
    command = 
      dd_get_rhat_1(
        dd_stan_1_fit = dd_stan_1_fit
      )
  ),
  tar_target(
    dd_data_model_evaluation_1,
    command = 
      dd_make_data_model_evaluation_1(
        dd_stan_1_fit = dd_stan_1_fit, 
        dd_stan_1_kfold = dd_stan_1_kfold, 
        dd_stan_1_kfold_2 = dd_stan_1_kfold_2, 
        dd_stan_1_kfold_cv_folds = dd_stan_1_kfold_cv_folds, 
        dd_stan_1_kfold_cv_folds_2 = dd_stan_1_kfold_cv_folds_2, 
        dd_data_model = dd_data_model,
        dd_model_info = dd_model_info
      )
  ),
  #### mixtures ####
  tar_target(
    dd_model_gamma_mirs_from_gamma_1,
    pattern = map(dd_id_model),
    command = 
      dd_make_model_gamma_mirs_from_gamma(
        dd_data_model = dd_data_model, 
        id_model = dd_id_model
      ) |>
      list(),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ), 
  tar_target(
    dd_mm2_file,
    command = "mixing_model/stan_2.stan",
    format = "file"
  ),
  tar_target(
    dd_mm2_cmd,
    command =
      cmdstanr::cmdstan_model(
        dd_mm2_file,
        cpp_options = list(stan_threads = TRUE, "CXXFLAGS+= -O3 -march=native -mtune=native"),
        # stanc_options = list("O1"),
        force_recompile = TRUE
      )
  ),
  tar_target(
    dd_data_preparation_for_mixing_model_1,
    command =
      dd_make_data_preparation_for_mixing_model_1(
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils, 
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat
      )
  ),
  tar_target(
    dd_stan_data_mixing_model_1,
    pattern = map(dd_id_model),
    command = 
      dd_make_stan_data_mixing_model_1(
        dd_data_preparation_for_mixing_model_1 = dd_data_preparation_for_mixing_model_1,
        id_model = dd_id_model
      ) |>
      list(),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      )
  ),
  tar_target(
    dd_fit_stan_mixing_model_1,
    pattern = map(dd_id_model),
    command = 
      dd_make_fit_stan_mixing_model_1(
        dd_mm2_cmd = dd_mm2_cmd, 
        dd_stan_data_mixing_model_1 = dd_stan_data_mixing_model_1[[dd_id_model]], 
        file = paste0("mixing_model/stan_2_fit_", dd_id_model, ".rds")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_res_reconstructions_1,
    command = 
      dd_make_res_reconstructions_1(
        dd_data_preparation_for_mixing_model_1 = dd_data_preparation_for_mixing_model_1, 
        dd_fit_stan_mixing_model_1 = purrr::map(dd_fit_stan_mixing_model_1, readRDS), 
        file = "mixing_model/dd_res_reconstructions.rds"
      ),
    format = "file"
  ),
  #### plots ####
  tar_target(
    dd_plot_1,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_1(
        dd_model = dd_stan_1_fit[[dd_id_model]], 
        dd_data_model = dd_data_model, 
        dd_data_pmird_mineral = dd_data_pmird_mineral, 
        config = dd_mir_preprocessing_config[[dd_id_model]], 
        prediction_domain = dd_model_info$prediction_domain[[dd_id_model]], 
        file_plot = paste0("figures/dd_plot_1_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_data_pmird_peat_cores_yhat,
    pattern = map(dd_id_model),
    command =
      {
        file_name <- paste0("targets_rvars/dd_data_pmird_peat_cores_yhat_", dd_id_model, ".rds")
        irp_degree_of_decomposition(
          x = dd_data_pmird_peat_cores, 
          model = dd_stan_1_fit[[dd_id_model]], 
          config = dd_mir_preprocessing_config[[dd_id_model]], 
          prediction_domain = dd_model_info$prediction_domain[[dd_id_model]]
        ) |>
          dplyr::select(dplyr::all_of("id_sample") | dplyr::starts_with("degree_of_decomposition")) |>
          saveRDS_rvars(file_name)
        file_name
      },
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_2,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_2(
        dd_model = dd_stan_1_fit[[dd_id_model]], 
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        file_plot = paste0("figures/dd_plot_2_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_3,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_3(
        dd_model = dd_stan_1_fit[[dd_id_model]], 
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        file_plot = paste0("figures/dd_plot_3_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_4,
    command =
      dd_make_plot_4(
        dd_model = dd_stan_1_fit[c(1, 3)], 
        dd_data_pmird_peat_cores = dd_data_pmird_peat_cores, 
        dd_data_pmird_peat_cores_yhat = 
          purrr::map(dd_data_pmird_peat_cores_yhat[c(1, 3)], readRDS_rvars), 
        file_plot = "figures/dd_plot_4.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_5,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_5(
        dd_model = dd_stan_1_fit[[dd_id_model]], 
        dd_data_model = dd_data_model,
        dd_data_pmird_decomposed = dd_data_pmird_decomposed,
        dd_data_pmird_peat_cores = dd_data_pmird_peat_cores, 
        config = dd_mir_preprocessing_config[[dd_id_model]], 
        prediction_domain = dd_model_info$prediction_domain[[dd_id_model]], 
        file_plot = paste0("figures/dd_plot_5_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_1,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = sym("macrofossil_volume_fraction_sphagnum_hummock"), 
        variable_color_legend_title = "Abundance of hummock<br><i>Sphagnum</i> (%)",
        file_plot = paste0("figures/dd_plot_6_1_", dd_id_model, ".pdf")
      ),
    resources =
     tar_resources(
       crew = tar_resources_crew(controller = "crew_parallel_1")
     ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_2,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat,
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = expression(carbon_content_1), 
        variable_color_legend_title = "C",
        file_plot = paste0("figures/dd_plot_6_2_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_3,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = expression(macrofossil_volume_fraction_wood), 
        variable_color_legend_title = "Abundance of woody remains (%)",
        file_plot = paste0("figures/dd_plot_6_3_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_4,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = expression(macrofossil_volume_fraction_shrub + macrofossil_volume_fraction_sedge + macrofossil_volume_fraction_sphagnum), 
        variable_color_legend_title = "Abundance of wood, sedge, and <i>Sphagnum</i> macrofossils (%)",
        file_plot = paste0("figures/dd_plot_6_4_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_5,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = expression(macrofossil_volume_fraction_shrub), 
        variable_color_legend_title = "Shrub macrofossil abundance (%)",
        file_plot = paste0("figures/dd_plot_6_5_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_6_6,
    pattern = map(dd_id_model),
    command =
      dd_make_plot_6(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils,
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[dd_id_model]]), 
        variable_color = expression(macrofossil_volume_fraction_sedge), 
        variable_color_legend_title = "Sedge macrofossil abundance (%)",
        file_plot = paste0("figures/dd_plot_6_6_", dd_id_model, ".pdf")
      ),
    resources =
      tar_resources(
        crew = tar_resources_crew(controller = "crew_parallel_1")
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_8,
    command = 
      dd_make_plot_8(
        dd_stan_1_fit = dd_stan_1_fit, 
        dd_data_model = dd_data_model, 
        file_plot = "figures/dd_plot_8.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_9,
    command = 
      dd_make_plot_9(
        dd_stan_1_fit = dd_stan_1_fit, 
        dd_model_info = dd_model_info, 
        file_plot = "figures/dd_plot_9.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_10_1,
    command =
      dd_make_plot_10(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat,
        dd_data_pmird_peat_cores_yhat = purrr::map(dd_data_pmird_peat_cores_yhat[-2], readRDS_rvars), 
        dd_model_info = dd_model_info[-2, ], 
        file_plot = "figures/dd_plot_10_1.pdf",
        ggsave_width = 11, 
        ggsave_height = 8
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_10_2,
    command =
      dd_make_plot_10(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat,
        dd_data_pmird_peat_cores_yhat = purrr::map(dd_data_pmird_peat_cores_yhat, readRDS_rvars), 
        dd_model_info = dd_model_info, 
        file_plot = "figures/dd_plot_10_2.pdf",
        ggsave_width = 11, 
        ggsave_height = 12
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_11,
    command = 
      dd_make_plot_11(
        dd_data_model_irpeat = dd_data_model_irpeat, 
        kfold = dd_stan_1_kfold_2, 
        dd_mir_preprocessing_config = dd_mir_preprocessing_config, 
        dd_model_info = dd_model_info, 
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat,
        file_plot = "figures/dd_plot_11.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_12,
    command = 
      dd_make_plot_12(
        dd_data_model_irpeat = dd_data_model_irpeat,
        dd_stan_1_fit = dd_stan_1_fit[[3]],
        file_plot = "figures/dd_plot_12.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_13,
    command = 
      dd_make_plot_13(
        dd_simulation_1 = dd_simulation_1,
        file_plot = "figures/dd_plot_13.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_14,
    command =
      dd_make_plot_14(
        dd_data_pmird_peat_cores = dd_data_pmird_peat_cores, 
        dd_mir_preprocessing_config = dd_mir_preprocessing_config,
        dd_model_info = dd_model_info, 
        dd_data_pmird_peat_cores_yhat = readRDS_rvars(dd_data_pmird_peat_cores_yhat[[3]]), 
        dd_stan_1_fit = dd_stan_1_fit, 
        file_plot = "figures/dd_plot_14.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_15,
    command = 
      dd_make_plot_15(
        dd_data_model_evaluation_1 = dd_data_model_evaluation_1,
        rmse_df = list(dd_stan_1_fit_rmse, dd_stan_1_kfold_2_rmse, dd_stan_1_kfold_rmse),
        dd_model_info,
        file_plot = "figures/dd_plot_15.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_16,
    command = 
      dd_make_plot_16(
        dd_data_pmird_mineral = dd_data_pmird_mineral, 
        file_plot = "figures/dd_plot_16.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_17,
    command = 
      dd_make_plot_17(
        dd_data_pmird_peat_cores = dd_data_pmird_peat_cores, 
        dd_data_pmird_peat_cores_yhat = dd_data_pmird_peat_cores_yhat,
        file_plot = "figures/dd_plot_17.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_18,
    command = 
      dd_make_plot_18(
        dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
        dd_data_pmird_peat_cores_yhat = dd_data_pmird_peat_cores_yhat,
        file_plot = "figures/dd_plot_18.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_19_1,
    command = 
      dd_make_plot_19(
        dd_simulation_1 = dd_simulation_2,
        mass_fraction_threshold = 0.9, 
        error_threshold = 0.2, 
        file_plot = "figures/dd_plot_19_1.pdf"
      )
  ),
  # tar_target(
  #   dd_plot_20,
  #   command =
  #     dd_make_plot_20(
  #       dd_reconstruction_initial_mass_rate_1 = readRDS_rvars(dd_reconstruction_initial_mass_rate_1), 
  #       file_plot = "figures/dd_plot_20.pdf"
  #     ),
  #   format = "file"
  # ),
  tar_target(
    dd_plot_21,
    command = 
      dd_make_plot_21(
        dd_rbacon_csv = dd_rbacon_csv
      )
  ),
  tar_target(
    dd_plot_22,
    command = 
      dd_make_plot_22(
        dd_simulation_2 = dd_simulation_2,
        file_plot = "figures/dd_plot_22.pdf"
      ),
    format = "file"
  ),
  # tar_target(
  #   dd_plot_23_1,
  #   command = 
  #     dd_make_plot_23(
  #       dd_simulation_3 = dd_simulation_3, 
  #       dd_simulation_3_species = dd_simulation_3_species,
  #       variable_degree_of_decomposition = sym("degree_of_decomposition_1"), 
  #       file_plot = "figures/dd_plot_23_1.pdf"
  #     ),
  #   format = "file"
  # ),
  # tar_target(
  #   dd_plot_23_2,
  #   command = 
  #     dd_make_plot_23(
  #       dd_simulation_3 = dd_simulation_3, 
  #       dd_simulation_3_species = dd_simulation_3_species,
  #       variable_degree_of_decomposition = sym("degree_of_decomposition_2"), 
  #       file_plot = "figures/dd_plot_23_2.pdf"
  #     ),
  #   format = "file"
  # ),
  # tar_target(
  #   dd_plot_23_3,
  #   command = 
  #     dd_make_plot_23(
  #       dd_simulation_3 = dd_simulation_3, 
  #       dd_simulation_3_species = dd_simulation_3_species,
  #       variable_degree_of_decomposition = sym("degree_of_decomposition_3"), 
  #       file_plot = "figures/dd_plot_23_3.pdf"
  #     ),
  #   format = "file"
  # ),
  tar_target(
    dd_plot_24,
    command = 
      dd_make_plot_24(
        dd_simulation_2 = dd_simulation_2,
        file_plot = "figures/dd_plot_24.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_25,
    command = 
      dd_make_plot_25(
        dd_simulation_4 = dd_simulation_4,
        file_plot = "figures/dd_plot_25.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_26,
    command = 
      dd_make_plot_26(
        dd_simulation_4 = dd_simulation_4,
        file_plot = "figures/dd_plot_26.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_27,
    command = 
      dd_make_plot_27(
        dd_simulation_3 = dd_simulation_3, 
        dd_model_gamma_mirs_from_gamma_1 = dd_model_gamma_mirs_from_gamma_1[[1]],
        file_plot = "figures/dd_plot_27.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_28,
    command = 
      dd_make_plot_28(
        dd_res_reconstructions_1 = readRDS_rvars(dd_res_reconstructions_1),
        dd_data_ta_wtd = dd_data_ta_wtd,
        file_plot = "figures/dd_plot_28.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_29,
    command =
      irp_make_plot_29(
        dd_res_reconstructions_1 = readRDS_rvars(dd_res_reconstructions_1), 
        file_plot = "figures/dd_plot_29.pdf"
      ),
    format = "file"
  ),
  tar_target(
    dd_plot_30,
    command =
      irp_make_plot_30(
        dd_res_reconstructions_1 = readRDS_rvars(dd_res_reconstructions_1), 
        file_plot = "figures/dd_plot_30.pdf"
      ),
    format = "file"
  ),
  #### reconstruction of NPP ####
  tar_target(
    dd_simulation_1,
    command =
      dd_make_simulation_1()
  ),
  # tar_target(
  #   dd_reconstruction_initial_mass_rate_1,
  #   command =
  #     dd_make_reconstruction_initial_mass_rate_1(
  #       dd_data_pmird_peat_cores_yhat = purrr::map(dd_data_pmird_peat_cores_yhat, readRDS_rvars), 
  #       dd_data_pmird_peat_cores_irpeat = dd_data_pmird_peat_cores_irpeat, 
  #       dd_data_pmird_peat_cores_macrofossils = dd_data_pmird_peat_cores_macrofossils, 
  #       dd_rbacon = dd_rbacon,
  #       res_file = "targets_rvars/dd_reconstruction_initial_mass_rate_1.rds"
  #     ),
  #   format = "file"
  # ),
  tar_target(
    dd_data_bona2018_moss_npp_summary_1,
    command = 
      dd_make_data_bona2018_moss_npp_summary_1(
        dd_data_bona2018_moss_npp = dd_data_bona2018_moss_npp#, 
        #dd_reconstruction_initial_mass_rate_1 = readRDS_rvars(dd_reconstruction_initial_mass_rate_1)
      )
  ),
  tar_target(
    dd_data_bengtsson2021_moss_npp_summary_1,
    command = 
      dd_make_data_bengtsson2021_moss_npp_summary_1(
        dd_data_bengtsson2021_moss_npp = dd_data_bengtsson2021_moss_npp
      )
  ),
  #### report ####
  tar_target(
    dd_citations_data,
    dd_make_citations_data(
      dd_data_model = dd_data_model,
      dd_data_pmird_peat_cores = dd_data_pmird_peat_cores
    )
  ),
  tar_target(
    dd_references_file,
    command = "references.bib",
    format = "file"
  ),
  tar_render(
   dd_supporting_info,
   path = "dd-supporting-info.Rmd"
  ),
  tar_render(
    dd_paper,
    path = "dd-paper.Rmd"
  )#,
  # tar_render(
  #   dd_bias,
  #   path = "dd-bias.Rmd"
  # ),
  # #### for projpred ####
  # tar_target(
  #   dd_model_info_2,
  #   dd_get_model_info_2(
  #     dd_data_model_preprocessed = dd_data_model_preprocessed, 
  #     dd_mir_preprocessing_settings = dd_mir_preprocessing_settings
  #   )
  # ),
  # tar_target(
  #   dd_brms_compiled_model_2,
  #   dd_make_brms_compiled_model_2(
  #     dd_model_info_2[1, ],
  #     backend = "cmdstanr",
  #     sig_figs = 12L
  #   )
  # ),
  # tar_target(
  #   dd_stan_2_fit,
  #   pattern = map(dd_id_model),
  #   command = 
  #     update(
  #       dd_brms_compiled_model_2, 
  #       newdata = 
  #         dd_model_info_2 |> 
  #         dplyr::filter(id_model == dd_id_model) |>
  #         dd_make_brms_data_2(),
  #       iter = dd_stan_1_mcmc_settings$iter,
  #       warmup = dd_stan_1_mcmc_settings$warmup,
  #       chains = dd_stan_1_mcmc_settings$chains,
  #       cores = dd_stan_1_mcmc_settings$chains,
  #       control = 
  #         list(
  #           max_treedepth = dd_stan_1_mcmc_settings$max_treedepth,
  #           adapt_delta = dd_stan_1_mcmc_settings$adapt_delta
  #         ),
  #       save_pars = brms::save_pars(all = TRUE),
  #       backend = "cmdstanr",
  #       save_warmup = TRUE,
  #       sig_figs = 14L
  #     ) |>
  #     list()
  # ),
  # tar_target(
  #   dd_stan_2_kfold_2,
  #   pattern = map(dd_id_model),
  #   command = 
  #     {
  #       if(dd_id_model == 1) {
  #         options(future.globals.maxSize = 1300 * 1024^2)
  #         future::plan("multicore")
  #         res <- 
  #           brms::kfold(
  #             x = dd_stan_2_fit[[dd_id_model]],
  #             chains = 4L,
  #             cores = 1L,
  #             folds = dd_stan_1_kfold_cv_folds_2,
  #             joint = "obs",
  #             save_fits = TRUE,
  #             recompile = FALSE
  #           ) |>
  #           list()
  #         options(future.globals.maxSize = 500 * 1024^2)
  #       } else {
  #         res <- NULL
  #       }
  #       res
  #     },
  #   resources =
  #     tar_resources(
  #       crew = tar_resources_crew(controller = "crew_parallel_2")
  #     )
  # )
)
