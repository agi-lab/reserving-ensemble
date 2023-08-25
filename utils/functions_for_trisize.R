# This file defines all the common functions to be used but have different implementations
#   for different triangle sizes

if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}
if (!tri.size %in% c(20, 40)){stop("Only triangles of size 20 and 40 are currently supported")}

if (tri.size == 40) {
    source('ensemble/defining_partitions_triangle_40.R')
    claims_df_train_val_test <- claims_df_train_val_test_40
    claims_df_in_out <- claims_df_in_out_40
    fit_all_component_models <- fit_all_component_models_40
    fit_all_partition_ensembles <- fit_all_partition_ensembles_40
    split_points <- split_points_40
    calc_crps_ensemble <- calc_crps_ensemble_40
    simulate_reserves <- simulate_reserves_40
    predict_component_partitions <- predict_component_partitions_40
    
} else if (tri.size == 20) {
    source('ensemble/defining_partitions_triangle_20.R')
    claims_df_train_val_test <- claims_df_train_val_test_20
    claims_df_in_out <- claims_df_in_out_20
    fit_all_component_models <- fit_all_component_models_20
    fit_all_partition_ensembles <- fit_all_partition_ensembles_20
    split_points <- split_points_20
    calc_crps_ensemble <- calc_crps_ensemble_20
    simulate_reserves <- simulate_reserves_20
    predict_component_partitions <- predict_component_partitions_20
    
} else {
    stop ("Error, parameter 'tri.size' is not defined correctly.")
}
