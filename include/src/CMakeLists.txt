project (relate)

add_library(relateStatic STATIC cxxopts filesystem plot fast_log collapsed_matrix data sample fast_painting anc mutations tree_builder branch_length_estimator anc_builder tree_comparer)
add_library(relateShared SHARED cxxopts filesystem plot fast_log collapsed_matrix data sample fast_painting anc mutations tree_builder branch_length_estimator anc_builder tree_comparer)
