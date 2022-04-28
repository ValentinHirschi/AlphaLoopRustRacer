import model aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':2}
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s

set_alphaLoop_option FORM_compile_cores 8
set_FORM_option cores 1
set_FORM_option FORM_parallel_cores 8

set_alphaLoop_option qgraf_model SM
set_FORM_option generate_integrated_UV_CTs False
set_FORM_option UV_min_dod_to_subtract 10
set_FORM_option on_shell_renormalisation True
set_FORM_option generate_renormalisation_graphs False
set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF

set_FORM_option optimize_c_output_mode_per_pair_of_factors "flat"
set_FORM_option max_n_lines_in_C_source 25000
set_FORM_option max_n_lines_in_C_function 1000

set_FORM_option reference_lmb {"SG_QG0":["pq2","pq4"]}
set_alphaLoop_option SG_name_list ["SG_QG0",]
qgraf_define j = d d~ g gh gh~
qgraf_generate a > d d~ / u s c b t QCD^2==1 QED^2==1 [QCD]
!rm -rf LU_a_ddx_NLO_SG_QG0_bare
output qgraf LU_a_ddx_NLO_SG_QG0_bare
