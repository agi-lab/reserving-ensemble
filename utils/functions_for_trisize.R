# This file defines all the common functions to be used but have different implementations
#   for different triangle sizes

if (!exists('tri.size')) {stop("Size of triangles 'tri.size' should be defined")}
if (!tri.size %in% c(10, 20, 40)){stop("Only triangles of size 10, 20 and 40 are currently supported")}

if (tri.size == 40) {
    source('ensemble/defining_partitions_triangle_40.R')
    claims_df_train_val_test <- claims_df_train_val_test_40
    claims_df_in_out <- claims_df_in_out_40
    fit_all_partition_ensembles <- fit_all_partition_ensembles_40
    split_points <- split_points_40
    
} else if (tri.size == 20) {
    source('ensemble/defining_partitions_triangle_20.R')
    claims_df_train_val_test <- claims_df_train_val_test_20
    claims_df_in_out <- claims_df_in_out_20
    fit_all_partition_ensembles <- fit_all_partition_ensembles_20
    split_points <- split_points_20
    
} else if (tri.size == 10) {
    source('ensemble/defining_partitions_triangle_10.R')
    claims_df_train_val_test <- claims_df_train_val_test_10
    claims_df_in_out <- claims_df_in_out_10
    fit_all_partition_ensembles <- fit_all_partition_ensembles_10
    split_points <- split_points_10
    
} else {
    stop ("Error, parameter 'tri.size' is not defined correcty.")
}