package org.biojava.bio.structure.io.mmcif.model;

import org.biojava.bio.structure.io.mmcif.chem.ChemCompTools;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;

public class ChemComp {
	String id ;
	String name;
	String type;
	String pdbx_type;
	String formula;
	String mon_nstd_parent_comp_id;
	String pdbx_synonyms;
	String pdbx_formal_charge;
	String pdbx_initial_date ;
	String pdbx_modified_date;
	String pdbx_ambiguous_flag;
	String pdbx_release_status ;
	String pdbx_replaced_by;
	String pdbx_replaces;
	String formula_weight;
	String one_letter_code;
	String three_letter_code;
	String pdbx_model_coordinates_details;
	String pdbx_model_coordinates_missing_flag;
	String pdbx_ideal_coordinates_details;
	String pdbx_ideal_coordinates_missing_flag;
	String pdbx_model_coordinates_db_code;
	String pdbx_subcomponent_list;
	String pdbx_processing_site;
	String mon_nstd_flag;

	// and some derived data for easier processing...
	ResidueType residueType;
	PolymerType polymerType;
	boolean standard;

	public String toString(){
		StringBuffer buf = new StringBuffer("ChemComp ");
		buf.append(id);
		buf.append(" ");
		buf.append(one_letter_code);
		buf.append(" ");
		buf.append(three_letter_code);
		buf.append(" poly:");
		buf.append(getPolymerType());
		buf.append(" resi:");
		buf.append(getResidueType());
		if (isStandard())
			buf.append(" standard");
		else
			buf.append(" modified");
		buf.append(" ");

		buf.append(name);
		buf.append(" ");
		buf.append(pdbx_type);
		buf.append(" ");
		buf.append(formula);
		buf.append(" parent:");
		buf.append(mon_nstd_parent_comp_id);
		return buf.toString();
	}

	public boolean hasParent(){
		String pid = mon_nstd_parent_comp_id;
		if ((pid != null ) && (! pid.equals("?"))){
			return true;
		}
		return false;
	}

	public boolean isStandard(){
		return standard;
	}

	private void setStandardFlag(){
		standard = ChemCompTools.isStandardChemComp(this);
	}



	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;

		residueType = ResidueType.getResidueTypeFromString(type);
		if ( residueType != null){
			polymerType = residueType.polymerType;
		}

	}


	public ResidueType getResidueType() {
		return residueType;
	}

	public void setResidueType(ResidueType residueType) {
		this.residueType = residueType;
	}

	public PolymerType getPolymerType() {
		return polymerType;
	}

	public void setPolymerType(PolymerType polymerType) {
		this.polymerType = polymerType;
	}

	public String getPdbx_type() {
		return pdbx_type;
	}
	public void setPdbx_type(String pdbx_type) {
		this.pdbx_type = pdbx_type;
	}
	public String getFormula() {
		return formula;
	}
	public void setFormula(String formula) {
		this.formula = formula;
	}
	public String getMon_nstd_parent_comp_id() {
		return mon_nstd_parent_comp_id;
	}
	public void setMon_nstd_parent_comp_id(String mon_nstd_parent_comp_id) {
		this.mon_nstd_parent_comp_id = mon_nstd_parent_comp_id;
		setStandardFlag();
	}
	public String getPdbx_synonyms() {
		return pdbx_synonyms;
	}
	public void setPdbx_synonyms(String pdbx_synonyms) {
		this.pdbx_synonyms = pdbx_synonyms;
	}
	public String getPdbx_formal_charge() {
		return pdbx_formal_charge;
	}
	public void setPdbx_formal_charge(String pdbx_formal_charge) {
		this.pdbx_formal_charge = pdbx_formal_charge;
	}
	public String getPdbx_initial_date() {
		return pdbx_initial_date;
	}
	public void setPdbx_initial_date(String pdbx_initial_date) {
		this.pdbx_initial_date = pdbx_initial_date;
	}
	public String getPdbx_modified_date() {
		return pdbx_modified_date;
	}
	public void setPdbx_modified_date(String pdbx_modified_date) {
		this.pdbx_modified_date = pdbx_modified_date;
	}
	public String getPdbx_ambiguous_flag() {
		return pdbx_ambiguous_flag;
	}
	public void setPdbx_ambiguous_flag(String pdbx_ambiguous_flag) {
		this.pdbx_ambiguous_flag = pdbx_ambiguous_flag;
	}
	public String getPdbx_release_status() {
		return pdbx_release_status;
	}
	public void setPdbx_release_status(String pdbx_release_status) {
		this.pdbx_release_status = pdbx_release_status;
	}
	public String getPdbx_replaced_by() {
		return pdbx_replaced_by;
	}
	public void setPdbx_replaced_by(String pdbx_replaced_by) {
		this.pdbx_replaced_by = pdbx_replaced_by;
	}
	public String getPdbx_replaces() {
		return pdbx_replaces;
	}
	public void setPdbx_replaces(String pdbx_replaces) {
		this.pdbx_replaces = pdbx_replaces;
	}
	public String getFormula_weight() {
		return formula_weight;
	}
	public void setFormula_weight(String formula_weight) {
		this.formula_weight = formula_weight;
	}
	public String getOne_letter_code() {
		return one_letter_code;
	}
	public void setOne_letter_code(String one_letter_code) {
		this.one_letter_code = one_letter_code;
		setStandardFlag();
	}
	public String getThree_letter_code() {
		return three_letter_code;
	}
	public void setThree_letter_code(String three_letter_code) {
		this.three_letter_code = three_letter_code;
	}
	public String getPdbx_model_coordinates_details() {
		return pdbx_model_coordinates_details;
	}
	public void setPdbx_model_coordinates_details(
			String pdbx_model_coordinates_details) {
		this.pdbx_model_coordinates_details = pdbx_model_coordinates_details;
	}
	public String getPdbx_model_coordinates_missing_flag() {
		return pdbx_model_coordinates_missing_flag;
	}
	public void setPdbx_model_coordinates_missing_flag(
			String pdbx_model_coordinates_missing_flag) {
		this.pdbx_model_coordinates_missing_flag = pdbx_model_coordinates_missing_flag;
	}
	public String getPdbx_ideal_coordinates_details() {
		return pdbx_ideal_coordinates_details;
	}
	public void setPdbx_ideal_coordinates_details(
			String pdbx_ideal_coordinates_details) {
		this.pdbx_ideal_coordinates_details = pdbx_ideal_coordinates_details;
	}
	public String getPdbx_ideal_coordinates_missing_flag() {
		return pdbx_ideal_coordinates_missing_flag;
	}
	public void setPdbx_ideal_coordinates_missing_flag(
			String pdbx_ideal_coordinates_missing_flag) {
		this.pdbx_ideal_coordinates_missing_flag = pdbx_ideal_coordinates_missing_flag;
	}
	public String getPdbx_model_coordinates_db_code() {
		return pdbx_model_coordinates_db_code;
	}
	public void setPdbx_model_coordinates_db_code(
			String pdbx_model_coordinates_db_code) {
		this.pdbx_model_coordinates_db_code = pdbx_model_coordinates_db_code;
	}
	public String getPdbx_subcomponent_list() {
		return pdbx_subcomponent_list;
	}
	public void setPdbx_subcomponent_list(String pdbx_subcomponent_list) {
		this.pdbx_subcomponent_list = pdbx_subcomponent_list;
	}
	public String getPdbx_processing_site() {
		return pdbx_processing_site;
	}
	public void setPdbx_processing_site(String pdbx_processing_site) {
		this.pdbx_processing_site = pdbx_processing_site;
	}

	public void setStandard(boolean standard) {
		this.standard = standard;
	}

   public String getMon_nstd_flag()
   {
      return mon_nstd_flag;
   }

   public void setMon_nstd_flag(String mon_nstd_flag)
   {
      this.mon_nstd_flag = mon_nstd_flag;
   }


}
