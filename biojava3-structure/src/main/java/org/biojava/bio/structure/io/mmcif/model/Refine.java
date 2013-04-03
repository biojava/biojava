package org.biojava.bio.structure.io.mmcif.model;

public class Refine {
	String entry_id;
	String ls_number_reflns_obs;
	String ls_number_reflns_all;
	String pdbx_ls_sigma_I;
	String pdbx_ls_sigma_F;
	String pdbx_data_cutoff_high_absF;
	String pdbx_data_cutoff_low_absF  ;
	String pdbx_data_cutoff_high_rms_absF;
	String ls_d_res_low  ;
	String ls_d_res_high  ;
	String ls_percent_reflns_obs;
	String ls_R_factor_obs ;
	String ls_R_factor_all  ;
	String ls_R_factor_R_work;
	String ls_R_factor_R_free ;
	String ls_R_factor_R_free_error;
	String ls_R_factor_R_free_error_details;
	String ls_percent_reflns_R_free;
	String ls_number_reflns_R_free;
	String ls_number_parameters;
	String ls_number_restraints;
	String occupancy_min;
	String occupancy_max;
	String B_iso_mean;
	String[][] aniso_B;
	String solvent_model_details ;
	String solvent_model_param_ksol;
	String solvent_model_param_bsol;
	String pdbx_ls_cross_valid_method;
	String details;
	String pdbx_starting_model;
	String pdbx_method_to_determine_struct;
	String pdbx_isotropic_thermal_model;
	String pdbx_stereochemistry_target_values;
	String pdbx_stereochem_target_val_spec_case;
	String pdbx_R_Free_selection_details;
	String pdbx_overall_ESU_R;
	String pdbx_overall_ESU_R_Free;
	String overall_SU_ML;
	String overall_SU_B;
	String ls_redundancy_reflns_obs;
	String pdbx_overall_phase_error   ;
	String B_iso_min;
	String B_iso_max;
	String correlation_coeff_Fo_to_Fc;
	String correlation_coeff_Fo_to_Fc_free;
	String pdbx_solvent_vdw_probe_radii;
	String pdbx_solvent_ion_probe_radii;
	String pdbx_solvent_shrinkage_radii;
	String overall_SU_R_Cruickshank_DPI;
	String overall_SU_R_free;
	String ls_wR_factor_R_free;
	String ls_wR_factor_R_work;
	String overall_FOM_free_R_set;
	String overall_FOM_work_R_set;
	String pdbx_refine_id;
	String pdbx_diffrn_id;
	String pdbx_TLS_residual_ADP_flag;
	String pdbx_overall_SU_R_free_Cruickshank_DPI;
	String pdbx_overall_SU_R_Blow_DPI;
	String pdbx_overall_SU_R_free_Blow_DPI;

	public Refine(){
		aniso_B = new String[3][3];
	}

	public String getEntry_id() {
		return entry_id;
	}

	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}

	public String getLs_number_reflns_obs() {
		return ls_number_reflns_obs;
	}

	public void setLs_number_reflns_obs(String ls_number_reflns_obs) {
		this.ls_number_reflns_obs = ls_number_reflns_obs;
	}

	public String getLs_number_reflns_all() {
		return ls_number_reflns_all;
	}

	public void setLs_number_reflns_all(String ls_number_reflns_all) {
		this.ls_number_reflns_all = ls_number_reflns_all;
	}

	public String getPdbx_ls_sigma_I() {
		return pdbx_ls_sigma_I;
	}

	public void setPdbx_ls_sigma_I(String pdbx_ls_sigma_I) {
		this.pdbx_ls_sigma_I = pdbx_ls_sigma_I;
	}

	public String getPdbx_ls_sigma_F() {
		return pdbx_ls_sigma_F;
	}

	public void setPdbx_ls_sigma_F(String pdbx_ls_sigma_F) {
		this.pdbx_ls_sigma_F = pdbx_ls_sigma_F;
	}

	public String getPdbx_data_cutoff_high_absF() {
		return pdbx_data_cutoff_high_absF;
	}

	public void setPdbx_data_cutoff_high_absF(String pdbx_data_cutoff_high_absF) {
		this.pdbx_data_cutoff_high_absF = pdbx_data_cutoff_high_absF;
	}

	public String getPdbx_data_cutoff_low_absF() {
		return pdbx_data_cutoff_low_absF;
	}

	public void setPdbx_data_cutoff_low_absF(String pdbx_data_cutoff_low_absF) {
		this.pdbx_data_cutoff_low_absF = pdbx_data_cutoff_low_absF;
	}

	public String getPdbx_data_cutoff_high_rms_absF() {
		return pdbx_data_cutoff_high_rms_absF;
	}

	public void setPdbx_data_cutoff_high_rms_absF(
			String pdbx_data_cutoff_high_rms_absF) {
		this.pdbx_data_cutoff_high_rms_absF = pdbx_data_cutoff_high_rms_absF;
	}

	public String getLs_d_res_low() {
		return ls_d_res_low;
	}

	public void setLs_d_res_low(String ls_d_res_low) {
		this.ls_d_res_low = ls_d_res_low;
	}

	public String getLs_d_res_high() {
		return ls_d_res_high;
	}

	public void setLs_d_res_high(String ls_d_res_high) {
		this.ls_d_res_high = ls_d_res_high;
	}

	public String getLs_percent_reflns_obs() {
		return ls_percent_reflns_obs;
	}

	public void setLs_percent_reflns_obs(String ls_percent_reflns_obs) {
		this.ls_percent_reflns_obs = ls_percent_reflns_obs;
	}

	public String getLs_R_factor_obs() {
		return ls_R_factor_obs;
	}

	public void setLs_R_factor_obs(String ls_R_factor_obs) {
		this.ls_R_factor_obs = ls_R_factor_obs;
	}

	public String getLs_R_factor_all() {
		return ls_R_factor_all;
	}

	public void setLs_R_factor_all(String ls_R_factor_all) {
		this.ls_R_factor_all = ls_R_factor_all;
	}

	public String getLs_R_factor_R_work() {
		return ls_R_factor_R_work;
	}

	public void setLs_R_factor_R_work(String ls_R_factor_R_work) {
		this.ls_R_factor_R_work = ls_R_factor_R_work;
	}

	public String getLs_R_factor_R_free() {
		return ls_R_factor_R_free;
	}

	public void setLs_R_factor_R_free(String ls_R_factor_R_free) {
		this.ls_R_factor_R_free = ls_R_factor_R_free;
	}

	public String getLs_R_factor_R_free_error() {
		return ls_R_factor_R_free_error;
	}

	public void setLs_R_factor_R_free_error(String ls_R_factor_R_free_error) {
		this.ls_R_factor_R_free_error = ls_R_factor_R_free_error;
	}

	public String getLs_R_factor_R_free_error_details() {
		return ls_R_factor_R_free_error_details;
	}

	public void setLs_R_factor_R_free_error_details(
			String ls_R_factor_R_free_error_details) {
		this.ls_R_factor_R_free_error_details = ls_R_factor_R_free_error_details;
	}

	public String getLs_percent_reflns_R_free() {
		return ls_percent_reflns_R_free;
	}

	public void setLs_percent_reflns_R_free(String ls_percent_reflns_R_free) {
		this.ls_percent_reflns_R_free = ls_percent_reflns_R_free;
	}

	public String getLs_number_reflns_R_free() {
		return ls_number_reflns_R_free;
	}

	public void setLs_number_reflns_R_free(String ls_number_reflns_R_free) {
		this.ls_number_reflns_R_free = ls_number_reflns_R_free;
	}

	public String getLs_number_parameters() {
		return ls_number_parameters;
	}

	public void setLs_number_parameters(String ls_number_parameters) {
		this.ls_number_parameters = ls_number_parameters;
	}

	public String getLs_number_restraints() {
		return ls_number_restraints;
	}

	public void setLs_number_restraints(String ls_number_restraints) {
		this.ls_number_restraints = ls_number_restraints;
	}

	public String getOccupancy_min() {
		return occupancy_min;
	}

	public void setOccupancy_min(String occupancy_min) {
		this.occupancy_min = occupancy_min;
	}

	public String getOccupancy_max() {
		return occupancy_max;
	}

	public void setOccupancy_max(String occupancy_max) {
		this.occupancy_max = occupancy_max;
	}

	public String getB_iso_mean() {
		return B_iso_mean;
	}

	public void setB_iso_mean(String b_iso_mean) {
		B_iso_mean = b_iso_mean;
	}

	public String[][] getAniso_B() {
		return aniso_B;
	}

	public void setAniso_B(String[][] aniso_B) {
		this.aniso_B = aniso_B;
	}

	public String getSolvent_model_details() {
		return solvent_model_details;
	}

	public void setSolvent_model_details(String solvent_model_details) {
		this.solvent_model_details = solvent_model_details;
	}

	public String getSolvent_model_param_ksol() {
		return solvent_model_param_ksol;
	}

	public void setSolvent_model_param_ksol(String solvent_model_param_ksol) {
		this.solvent_model_param_ksol = solvent_model_param_ksol;
	}

	public String getSolvent_model_param_bsol() {
		return solvent_model_param_bsol;
	}

	public void setSolvent_model_param_bsol(String solvent_model_param_bsol) {
		this.solvent_model_param_bsol = solvent_model_param_bsol;
	}

	public String getPdbx_ls_cross_valid_method() {
		return pdbx_ls_cross_valid_method;
	}

	public void setPdbx_ls_cross_valid_method(String pdbx_ls_cross_valid_method) {
		this.pdbx_ls_cross_valid_method = pdbx_ls_cross_valid_method;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getPdbx_starting_model() {
		return pdbx_starting_model;
	}

	public void setPdbx_starting_model(String pdbx_starting_model) {
		this.pdbx_starting_model = pdbx_starting_model;
	}

	public String getPdbx_method_to_determine_struct() {
		return pdbx_method_to_determine_struct;
	}

	public void setPdbx_method_to_determine_struct(
			String pdbx_method_to_determine_struct) {
		this.pdbx_method_to_determine_struct = pdbx_method_to_determine_struct;
	}

	public String getPdbx_isotropic_thermal_model() {
		return pdbx_isotropic_thermal_model;
	}

	public void setPdbx_isotropic_thermal_model(String pdbx_isotropic_thermal_model) {
		this.pdbx_isotropic_thermal_model = pdbx_isotropic_thermal_model;
	}

	public String getPdbx_stereochemistry_target_values() {
		return pdbx_stereochemistry_target_values;
	}

	public void setPdbx_stereochemistry_target_values(
			String pdbx_stereochemistry_target_values) {
		this.pdbx_stereochemistry_target_values = pdbx_stereochemistry_target_values;
	}

	public String getPdbx_stereochem_target_val_spec_case() {
		return pdbx_stereochem_target_val_spec_case;
	}

	public void setPdbx_stereochem_target_val_spec_case(
			String pdbx_stereochem_target_val_spec_case) {
		this.pdbx_stereochem_target_val_spec_case = pdbx_stereochem_target_val_spec_case;
	}

	public String getPdbx_R_Free_selection_details() {
		return pdbx_R_Free_selection_details;
	}

	public void setPdbx_R_Free_selection_details(
			String pdbx_R_Free_selection_details) {
		this.pdbx_R_Free_selection_details = pdbx_R_Free_selection_details;
	}

	public String getPdbx_overall_ESU_R() {
		return pdbx_overall_ESU_R;
	}

	public void setPdbx_overall_ESU_R(String pdbx_overall_ESU_R) {
		this.pdbx_overall_ESU_R = pdbx_overall_ESU_R;
	}

	public String getPdbx_overall_ESU_R_Free() {
		return pdbx_overall_ESU_R_Free;
	}

	public void setPdbx_overall_ESU_R_Free(String pdbx_overall_ESU_R_Free) {
		this.pdbx_overall_ESU_R_Free = pdbx_overall_ESU_R_Free;
	}

	public String getOverall_SU_ML() {
		return overall_SU_ML;
	}

	public void setOverall_SU_ML(String overall_SU_ML) {
		this.overall_SU_ML = overall_SU_ML;
	}

	public String getOverall_SU_B() {
		return overall_SU_B;
	}

	public void setOverall_SU_B(String overall_SU_B) {
		this.overall_SU_B = overall_SU_B;
	}

	public String getPdbx_refine_id() {
		return pdbx_refine_id;
	}

	public void setPdbx_refine_id(String pdbx_refine_id) {
		this.pdbx_refine_id = pdbx_refine_id;
	}

	public String getLs_redundancy_reflns_obs() {
		return ls_redundancy_reflns_obs;
	}

	public void setLs_redundancy_reflns_obs(String ls_redundancy_reflns_obs) {
		this.ls_redundancy_reflns_obs = ls_redundancy_reflns_obs;
	}

	public String getPdbx_overall_phase_error() {
		return pdbx_overall_phase_error;
	}

	public void setPdbx_overall_phase_error(String pdbx_overall_phase_error) {
		this.pdbx_overall_phase_error = pdbx_overall_phase_error;
	}

	public String getB_iso_min() {
		return B_iso_min;
	}

	public void setB_iso_min(String b_iso_min) {
		B_iso_min = b_iso_min;
	}

	public String getB_iso_max() {
		return B_iso_max;
	}

	public void setB_iso_max(String b_iso_max) {
		B_iso_max = b_iso_max;
	}

	public String getCorrelation_coeff_Fo_to_Fc() {
		return correlation_coeff_Fo_to_Fc;
	}

	public void setCorrelation_coeff_Fo_to_Fc(String correlation_coeff_Fo_to_Fc) {
		this.correlation_coeff_Fo_to_Fc = correlation_coeff_Fo_to_Fc;
	}

	public String getCorrelation_coeff_Fo_to_Fc_free() {
		return correlation_coeff_Fo_to_Fc_free;
	}

	public void setCorrelation_coeff_Fo_to_Fc_free(
			String correlation_coeff_Fo_to_Fc_free) {
		this.correlation_coeff_Fo_to_Fc_free = correlation_coeff_Fo_to_Fc_free;
	}

	public String getPdbx_solvent_vdw_probe_radii() {
		return pdbx_solvent_vdw_probe_radii;
	}

	public void setPdbx_solvent_vdw_probe_radii(String pdbx_solvent_vdw_probe_radii) {
		this.pdbx_solvent_vdw_probe_radii = pdbx_solvent_vdw_probe_radii;
	}

	public String getPdbx_solvent_ion_probe_radii() {
		return pdbx_solvent_ion_probe_radii;
	}

	public void setPdbx_solvent_ion_probe_radii(String pdbx_solvent_ion_probe_radii) {
		this.pdbx_solvent_ion_probe_radii = pdbx_solvent_ion_probe_radii;
	}

	public String getPdbx_solvent_shrinkage_radii() {
		return pdbx_solvent_shrinkage_radii;
	}

	public void setPdbx_solvent_shrinkage_radii(String pdbx_solvent_shrinkage_radii) {
		this.pdbx_solvent_shrinkage_radii = pdbx_solvent_shrinkage_radii;
	}

	public String getOverall_SU_R_Cruickshank_DPI() {
		return overall_SU_R_Cruickshank_DPI;
	}

	public void setOverall_SU_R_Cruickshank_DPI(String overall_SU_R_Cruickshank_DPI) {
		this.overall_SU_R_Cruickshank_DPI = overall_SU_R_Cruickshank_DPI;
	}

	public String getOverall_SU_R_free() {
		return overall_SU_R_free;
	}

	public void setOverall_SU_R_free(String overall_SU_R_free) {
		this.overall_SU_R_free = overall_SU_R_free;
	}

	public String getLs_wR_factor_R_free() {
		return ls_wR_factor_R_free;
	}

	public void setLs_wR_factor_R_free(String ls_wR_factor_R_free) {
		this.ls_wR_factor_R_free = ls_wR_factor_R_free;
	}

	public String getLs_wR_factor_R_work() {
		return ls_wR_factor_R_work;
	}

	public void setLs_wR_factor_R_work(String ls_wR_factor_R_work) {
		this.ls_wR_factor_R_work = ls_wR_factor_R_work;
	}

	public String getOverall_FOM_free_R_set() {
		return overall_FOM_free_R_set;
	}

	public void setOverall_FOM_free_R_set(String overall_FOM_free_R_set) {
		this.overall_FOM_free_R_set = overall_FOM_free_R_set;
	}

	public String getOverall_FOM_work_R_set() {
		return overall_FOM_work_R_set;
	}

	public void setOverall_FOM_work_R_set(String overall_FOM_work_R_set) {
		this.overall_FOM_work_R_set = overall_FOM_work_R_set;
	}

	public String getPdbx_diffrn_id() {
		return pdbx_diffrn_id;
	}

	public void setPdbx_diffrn_id(String pdbx_diffrn_id) {
		this.pdbx_diffrn_id = pdbx_diffrn_id;
	}

	public String getPdbx_TLS_residual_ADP_flag() {
		return pdbx_TLS_residual_ADP_flag;
	}

	public void setPdbx_TLS_residual_ADP_flag(String pdbx_TLS_residual_ADP_flag) {
		this.pdbx_TLS_residual_ADP_flag = pdbx_TLS_residual_ADP_flag;
	}

	public String getPdbx_overall_SU_R_free_Cruickshank_DPI() {
		return pdbx_overall_SU_R_free_Cruickshank_DPI;
	}

	public void setPdbx_overall_SU_R_free_Cruickshank_DPI(
			String pdbx_overall_SU_R_free_Cruickshank_DPI) {
		this.pdbx_overall_SU_R_free_Cruickshank_DPI = pdbx_overall_SU_R_free_Cruickshank_DPI;
	}

	public String getPdbx_overall_SU_R_Blow_DPI() {
		return pdbx_overall_SU_R_Blow_DPI;
	}

	public void setPdbx_overall_SU_R_Blow_DPI(String pdbx_overall_SU_R_Blow_DPI) {
		this.pdbx_overall_SU_R_Blow_DPI = pdbx_overall_SU_R_Blow_DPI;
	}

	public String getPdbx_overall_SU_R_free_Blow_DPI() {
		return pdbx_overall_SU_R_free_Blow_DPI;
	}

	public void setPdbx_overall_SU_R_free_Blow_DPI(
			String pdbx_overall_SU_R_free_Blow_DPI) {
		this.pdbx_overall_SU_R_free_Blow_DPI = pdbx_overall_SU_R_free_Blow_DPI;
	}
	

}
