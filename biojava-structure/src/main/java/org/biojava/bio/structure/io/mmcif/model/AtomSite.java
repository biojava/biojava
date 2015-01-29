/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * created at Apr 26, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

public class AtomSite extends AbstractBean{
	String group_PDB;
	String id;
	String type_symbol;
	String label_atom_id;
	String label_alt_id;
	String label_comp_id;
	String label_asym_id;
	String label_entity_id;
	String label_seq_id;
	String pdbx_PDB_ins_code;
	String Cartn_x;
	String Cartn_y;
	String Cartn_z;
	String occupancy;
	String B_iso_or_equiv;
	String Cartn_x_esd;
	String Cartn_y_esd;
	String Cartn_z_esd;
	String auth_seq_id;
	String auth_comp_id;
	String auth_asym_id;
	String auth_atom_id;
	String pdbx_PDB_model_num;
	String occupancy_esd;
	String B_iso_or_equiv_esd;
	String pdbx_formal_charge;

	public String getGroup_PDB() {
		return group_PDB;
	}
	public void setGroup_PDB(String group_PDB) {
		this.group_PDB = group_PDB;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getType_symbol() {
		return type_symbol;
	}
	public void setType_symbol(String type_symbol) {
		this.type_symbol = type_symbol;
	}

	public String getLabel_alt_id() {
		return label_alt_id;
	}
	public void setLabel_alt_id(String label_alt_id) {
		this.label_alt_id = label_alt_id;
	}
	public String getLabel_comp_id() {
		return label_comp_id;
	}
	public void setLabel_comp_id(String label_comp_id) {
		this.label_comp_id = label_comp_id;
	}
	public String getLabel_entity_id() {
		return label_entity_id;
	}
	public void setLabel_entity_id(String label_entity_id) {
		this.label_entity_id = label_entity_id;
	}
	public String getLabel_seq_id() {
		return label_seq_id;
	}
	public void setLabel_seq_id(String label_seq_id) {
		this.label_seq_id = label_seq_id;
	}
	public String getPdbx_PDB_ins_code() {
		return pdbx_PDB_ins_code;
	}
	public void setPdbx_PDB_ins_code(String pdbx_PDB_ins_code) {
		this.pdbx_PDB_ins_code = pdbx_PDB_ins_code;
	}
	public String getCartn_x() {
		return Cartn_x;
	}
	public void setCartn_x(String cartn_x) {
		Cartn_x = cartn_x;
	}
	public String getCartn_y() {
		return Cartn_y;
	}
	public void setCartn_y(String cartn_y) {
		Cartn_y = cartn_y;
	}
	public String getCartn_z() {
		return Cartn_z;
	}
	public void setCartn_z(String cartn_z) {
		Cartn_z = cartn_z;
	}
	public String getOccupancy() {
		return occupancy;
	}
	public void setOccupancy(String occupancy) {
		this.occupancy = occupancy;
	}
	public String getB_iso_or_equiv() {
		return B_iso_or_equiv;
	}
	public void setB_iso_or_equiv(String b_iso_or_equiv) {
		B_iso_or_equiv = b_iso_or_equiv;
	}
	public String getCartn_x_esd() {
		return Cartn_x_esd;
	}
	public void setCartn_x_esd(String cartn_x_esd) {
		Cartn_x_esd = cartn_x_esd;
	}
	public String getCartn_y_esd() {
		return Cartn_y_esd;
	}
	public void setCartn_y_esd(String cartn_y_esd) {
		Cartn_y_esd = cartn_y_esd;
	}
	public String getCartn_z_esd() {
		return Cartn_z_esd;
	}
	public void setCartn_z_esd(String cartn_z_esd) {
		Cartn_z_esd = cartn_z_esd;
	}
	public String getAuth_seq_id() {
		return auth_seq_id;
	}
	public void setAuth_seq_id(String auth_seq_id) {
		this.auth_seq_id = auth_seq_id;
	}
	public String getAuth_comp_id() {
		return auth_comp_id;
	}
	public void setAuth_comp_id(String auth_comp_id) {
		this.auth_comp_id = auth_comp_id;
	}
	public String getAuth_asym_id() {
		return auth_asym_id;
	}
	public void setAuth_asym_id(String auth_asym_id) {
		this.auth_asym_id = auth_asym_id;
	}
	public String getAuth_atom_id() {
		return auth_atom_id;
	}
	public void setAuth_atom_id(String auth_atom_id) {
		this.auth_atom_id = auth_atom_id;
	}
	public String getPdbx_PDB_model_num() {
		return pdbx_PDB_model_num;
	}
	public void setPdbx_PDB_model_num(String pdbx_PDB_model_num) {
		this.pdbx_PDB_model_num = pdbx_PDB_model_num;
	}
	public String getLabel_atom_id() {
		return label_atom_id;
	}
	public void setLabel_atom_id(String label_atom_id) {
		this.label_atom_id = label_atom_id;
	}
	public String getLabel_asym_id() {
		return label_asym_id;
	}
	public void setLabel_asym_id(String label_asym_id) {
		this.label_asym_id = label_asym_id;
	}
	public String getOccupancy_esd() {
		return occupancy_esd;
	}
	public void setOccupancy_esd(String occupancy_esd) {
		this.occupancy_esd = occupancy_esd;
	}
	public String getB_iso_or_equiv_esd() {
		return B_iso_or_equiv_esd;
	}
	public void setB_iso_or_equiv_esd(String b_iso_or_equiv_esd) {
		B_iso_or_equiv_esd = b_iso_or_equiv_esd;
	}
	public String getPdbx_formal_charge() {
		return pdbx_formal_charge;
	}
	public void setPdbx_formal_charge(String pdbx_formal_charge) {
		this.pdbx_formal_charge = pdbx_formal_charge;
	}




}
