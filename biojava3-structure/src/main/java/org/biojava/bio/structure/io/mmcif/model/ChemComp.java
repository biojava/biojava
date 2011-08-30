package org.biojava.bio.structure.io.mmcif.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.io.mmcif.chem.ChemCompTools;
import org.biojava.bio.structure.io.mmcif.chem.PolymerType;
import org.biojava.bio.structure.io.mmcif.chem.ResidueType;

/** A definition for a Chemical Component, as maintained by the wwPDB. For access to all definitions,
 * please download the components.cif.gz file from the wwPDB website.
 * 
 * @author Andreas Prlic
 *
 */
public class ChemComp implements Serializable, Comparable<ChemComp>{
	/**
	 * 
	 */
	private static final long serialVersionUID = -4736341142030215915L;

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

	List<ChemCompDescriptor> descriptors = new ArrayList<ChemCompDescriptor>();

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

	public List<ChemCompDescriptor> getDescriptors() {
		return descriptors;
	}

	public void setDescriptors(List<ChemCompDescriptor> descriptors) {
		this.descriptors = descriptors;
	}

	
	public int compareTo(ChemComp arg0) {
		if ( this.equals(arg0))
			return 0;
		return this.getId().compareTo(arg0.getId());
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((descriptors == null) ? 0 : descriptors.hashCode());
		result = prime * result + ((formula == null) ? 0 : formula.hashCode());
		result = prime * result
				+ ((formula_weight == null) ? 0 : formula_weight.hashCode());
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		result = prime * result
				+ ((mon_nstd_flag == null) ? 0 : mon_nstd_flag.hashCode());
		result = prime
				* result
				+ ((mon_nstd_parent_comp_id == null) ? 0
						: mon_nstd_parent_comp_id.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result
				+ ((one_letter_code == null) ? 0 : one_letter_code.hashCode());
		result = prime
				* result
				+ ((pdbx_ambiguous_flag == null) ? 0 : pdbx_ambiguous_flag
						.hashCode());
		result = prime
				* result
				+ ((pdbx_formal_charge == null) ? 0 : pdbx_formal_charge
						.hashCode());
		result = prime
				* result
				+ ((pdbx_ideal_coordinates_details == null) ? 0
						: pdbx_ideal_coordinates_details.hashCode());
		result = prime
				* result
				+ ((pdbx_ideal_coordinates_missing_flag == null) ? 0
						: pdbx_ideal_coordinates_missing_flag.hashCode());
		result = prime
				* result
				+ ((pdbx_initial_date == null) ? 0 : pdbx_initial_date
						.hashCode());
		result = prime
				* result
				+ ((pdbx_model_coordinates_db_code == null) ? 0
						: pdbx_model_coordinates_db_code.hashCode());
		result = prime
				* result
				+ ((pdbx_model_coordinates_details == null) ? 0
						: pdbx_model_coordinates_details.hashCode());
		result = prime
				* result
				+ ((pdbx_model_coordinates_missing_flag == null) ? 0
						: pdbx_model_coordinates_missing_flag.hashCode());
		result = prime
				* result
				+ ((pdbx_modified_date == null) ? 0 : pdbx_modified_date
						.hashCode());
		result = prime
				* result
				+ ((pdbx_processing_site == null) ? 0 : pdbx_processing_site
						.hashCode());
		result = prime
				* result
				+ ((pdbx_release_status == null) ? 0 : pdbx_release_status
						.hashCode());
		result = prime
				* result
				+ ((pdbx_replaced_by == null) ? 0 : pdbx_replaced_by.hashCode());
		result = prime * result
				+ ((pdbx_replaces == null) ? 0 : pdbx_replaces.hashCode());
		result = prime
				* result
				+ ((pdbx_subcomponent_list == null) ? 0
						: pdbx_subcomponent_list.hashCode());
		result = prime * result
				+ ((pdbx_synonyms == null) ? 0 : pdbx_synonyms.hashCode());
		result = prime * result
				+ ((pdbx_type == null) ? 0 : pdbx_type.hashCode());
		result = prime * result
				+ ((polymerType == null) ? 0 : polymerType.hashCode());
		result = prime * result
				+ ((residueType == null) ? 0 : residueType.hashCode());
		result = prime * result + (standard ? 1231 : 1237);
		result = prime
				* result
				+ ((three_letter_code == null) ? 0 : three_letter_code
						.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ChemComp other = (ChemComp) obj;
		if (descriptors == null) {
			if (other.descriptors != null)
				return false;
		} else if (!descriptors.equals(other.descriptors))
			return false;
		if (formula == null) {
			if (other.formula != null)
				return false;
		} else if (!formula.equals(other.formula))
			return false;
		if (formula_weight == null) {
			if (other.formula_weight != null)
				return false;
		} else if (!formula_weight.equals(other.formula_weight))
			return false;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		if (mon_nstd_flag == null) {
			if (other.mon_nstd_flag != null)
				return false;
		} else if (!mon_nstd_flag.equals(other.mon_nstd_flag))
			return false;
		if (mon_nstd_parent_comp_id == null) {
			if (other.mon_nstd_parent_comp_id != null)
				return false;
		} else if (!mon_nstd_parent_comp_id
				.equals(other.mon_nstd_parent_comp_id))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (one_letter_code == null) {
			if (other.one_letter_code != null)
				return false;
		} else if (!one_letter_code.equals(other.one_letter_code))
			return false;
		if (pdbx_ambiguous_flag == null) {
			if (other.pdbx_ambiguous_flag != null)
				return false;
		} else if (!pdbx_ambiguous_flag.equals(other.pdbx_ambiguous_flag))
			return false;
		if (pdbx_formal_charge == null) {
			if (other.pdbx_formal_charge != null)
				return false;
		} else if (!pdbx_formal_charge.equals(other.pdbx_formal_charge))
			return false;
		if (pdbx_ideal_coordinates_details == null) {
			if (other.pdbx_ideal_coordinates_details != null)
				return false;
		} else if (!pdbx_ideal_coordinates_details
				.equals(other.pdbx_ideal_coordinates_details))
			return false;
		if (pdbx_ideal_coordinates_missing_flag == null) {
			if (other.pdbx_ideal_coordinates_missing_flag != null)
				return false;
		} else if (!pdbx_ideal_coordinates_missing_flag
				.equals(other.pdbx_ideal_coordinates_missing_flag))
			return false;
		if (pdbx_initial_date == null) {
			if (other.pdbx_initial_date != null)
				return false;
		} else if (!pdbx_initial_date.equals(other.pdbx_initial_date))
			return false;
		if (pdbx_model_coordinates_db_code == null) {
			if (other.pdbx_model_coordinates_db_code != null)
				return false;
		} else if (!pdbx_model_coordinates_db_code
				.equals(other.pdbx_model_coordinates_db_code))
			return false;
		if (pdbx_model_coordinates_details == null) {
			if (other.pdbx_model_coordinates_details != null)
				return false;
		} else if (!pdbx_model_coordinates_details
				.equals(other.pdbx_model_coordinates_details))
			return false;
		if (pdbx_model_coordinates_missing_flag == null) {
			if (other.pdbx_model_coordinates_missing_flag != null)
				return false;
		} else if (!pdbx_model_coordinates_missing_flag
				.equals(other.pdbx_model_coordinates_missing_flag))
			return false;
		if (pdbx_modified_date == null) {
			if (other.pdbx_modified_date != null)
				return false;
		} else if (!pdbx_modified_date.equals(other.pdbx_modified_date))
			return false;
		if (pdbx_processing_site == null) {
			if (other.pdbx_processing_site != null)
				return false;
		} else if (!pdbx_processing_site.equals(other.pdbx_processing_site))
			return false;
		if (pdbx_release_status == null) {
			if (other.pdbx_release_status != null)
				return false;
		} else if (!pdbx_release_status.equals(other.pdbx_release_status))
			return false;
		if (pdbx_replaced_by == null) {
			if (other.pdbx_replaced_by != null)
				return false;
		} else if (!pdbx_replaced_by.equals(other.pdbx_replaced_by))
			return false;
		if (pdbx_replaces == null) {
			if (other.pdbx_replaces != null)
				return false;
		} else if (!pdbx_replaces.equals(other.pdbx_replaces))
			return false;
		if (pdbx_subcomponent_list == null) {
			if (other.pdbx_subcomponent_list != null)
				return false;
		} else if (!pdbx_subcomponent_list.equals(other.pdbx_subcomponent_list))
			return false;
		if (pdbx_synonyms == null) {
			if (other.pdbx_synonyms != null)
				return false;
		} else if (!pdbx_synonyms.equals(other.pdbx_synonyms))
			return false;
		if (pdbx_type == null) {
			if (other.pdbx_type != null)
				return false;
		} else if (!pdbx_type.equals(other.pdbx_type))
			return false;
		if (polymerType != other.polymerType)
			return false;
		if (residueType != other.residueType)
			return false;
		if (standard != other.standard)
			return false;
		if (three_letter_code == null) {
			if (other.three_letter_code != null)
				return false;
		} else if (!three_letter_code.equals(other.three_letter_code))
			return false;
		if (type == null) {
			if (other.type != null)
				return false;
		} else if (!type.equals(other.type))
			return false;
		return true;
	}


}
