/**
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
 * Created on Feb 5, 2013
 * Created by Andreas Prlic
 *
 * @since 3.0.6
 */
package org.biojava.bio.structure.io.mmcif.model;

import java.io.Serializable;

/** stores these fields:
 * 
 * _chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_align
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_component_comp_id
_chem_comp_atom.pdbx_residue_numbering
_chem_comp_atom.pdbx_component_atom_id
_chem_comp_atom.pdbx_polymer_type
_chem_comp_atom.pdbx_ref_id
_chem_comp_atom.pdbx_component_id
_chem_comp_atom.pdbx_ordinal
 * 
 * @author Andreas Prlic
 *
 */
public class ChemCompAtom implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 4070599340294758941L;
	String comp_id;
	String atom_id;
	String alt_atom_id;
	String type_symbol;
	String charge;
	String pdbx_align;
	String pdbx_aromatic_flag;
	String pdbx_leaving_atom_flag;
	String pdbx_stereo_config;
	String model_Cartn_x;
	String model_Cartn_y;
	String model_Cartn_z;
	String pdbx_model_Cartn_x_ideal;
	String pdbx_model_Cartn_y_ideal;
	String pdbx_model_Cartn_z_ideal;
	String pdbx_component_comp_id;
	String pdbx_residue_numbering;
	String pdbx_component_atom_id;
	String pdbx_polymer_type;
	String pdbx_ref_id;
	String pdbx_component_id;
	String pdbx_ordinal;
	public String getComp_id() {
		return comp_id;
	}
	public void setComp_id(String comp_id) {
		this.comp_id = comp_id;
	}
	public String getAtom_id() {
		return atom_id;
	}
	public void setAtom_id(String atom_id) {
		this.atom_id = atom_id;
	}
	public String getAlt_atom_id() {
		return alt_atom_id;
	}
	public void setAlt_atom_id(String alt_atom_id) {
		this.alt_atom_id = alt_atom_id;
	}
	public String getType_symbol() {
		return type_symbol;
	}
	public void setType_symbol(String type_symbol) {
		this.type_symbol = type_symbol;
	}
	public String getCharge() {
		return charge;
	}
	public void setCharge(String charge) {
		this.charge = charge;
	}
	public String getPdbx_align() {
		return pdbx_align;
	}
	public void setPdbx_align(String pdbx_align) {
		this.pdbx_align = pdbx_align;
	}
	public String getPdbx_aromatic_flag() {
		return pdbx_aromatic_flag;
	}
	public void setPdbx_aromatic_flag(String pdbx_aromatic_flag) {
		this.pdbx_aromatic_flag = pdbx_aromatic_flag;
	}
	public String getPdbx_leaving_atom_flag() {
		return pdbx_leaving_atom_flag;
	}
	public void setPdbx_leaving_atom_flag(String pdbx_leaving_atom_flag) {
		this.pdbx_leaving_atom_flag = pdbx_leaving_atom_flag;
	}
	public String getPdbx_stereo_config() {
		return pdbx_stereo_config;
	}
	public void setPdbx_stereo_config(String pdbx_stereo_config) {
		this.pdbx_stereo_config = pdbx_stereo_config;
	}
	public String getModel_Cartn_x() {
		return model_Cartn_x;
	}
	public void setModel_Cartn_x(String model_Cartn_x) {
		this.model_Cartn_x = model_Cartn_x;
	}
	public String getModel_Cartn_y() {
		return model_Cartn_y;
	}
	public void setModel_Cartn_y(String model_Cartn_y) {
		this.model_Cartn_y = model_Cartn_y;
	}
	public String getModel_Cartn_z() {
		return model_Cartn_z;
	}
	public void setModel_Cartn_z(String model_Cartn_z) {
		this.model_Cartn_z = model_Cartn_z;
	}
	public String getPdbx_model_Cartn_x_ideal() {
		return pdbx_model_Cartn_x_ideal;
	}
	public void setPdbx_model_Cartn_x_ideal(String pdbx_model_Cartn_x_ideal) {
		this.pdbx_model_Cartn_x_ideal = pdbx_model_Cartn_x_ideal;
	}
	public String getPdbx_model_Cartn_y_ideal() {
		return pdbx_model_Cartn_y_ideal;
	}
	public void setPdbx_model_Cartn_y_ideal(String pdbx_model_Cartn_y_ideal) {
		this.pdbx_model_Cartn_y_ideal = pdbx_model_Cartn_y_ideal;
	}
	public String getPdbx_model_Cartn_z_ideal() {
		return pdbx_model_Cartn_z_ideal;
	}
	public void setPdbx_model_Cartn_z_ideal(String pdbx_model_Cartn_z_ideal) {
		this.pdbx_model_Cartn_z_ideal = pdbx_model_Cartn_z_ideal;
	}
	public String getPdbx_component_comp_id() {
		return pdbx_component_comp_id;
	}
	public void setPdbx_component_comp_id(String pdbx_component_comp_id) {
		this.pdbx_component_comp_id = pdbx_component_comp_id;
	}
	public String getPdbx_residue_numbering() {
		return pdbx_residue_numbering;
	}
	public void setPdbx_residue_numbering(String pdbx_residue_numbering) {
		this.pdbx_residue_numbering = pdbx_residue_numbering;
	}
	public String getPdbx_component_atom_id() {
		return pdbx_component_atom_id;
	}
	public void setPdbx_component_atom_id(String pdbx_component_atom_id) {
		this.pdbx_component_atom_id = pdbx_component_atom_id;
	}
	public String getPdbx_polymer_type() {
		return pdbx_polymer_type;
	}
	public void setPdbx_polymer_type(String pdbx_polymer_type) {
		this.pdbx_polymer_type = pdbx_polymer_type;
	}
	public String getPdbx_ref_id() {
		return pdbx_ref_id;
	}
	public void setPdbx_ref_id(String pdbx_ref_id) {
		this.pdbx_ref_id = pdbx_ref_id;
	}
	public String getPdbx_component_id() {
		return pdbx_component_id;
	}
	public void setPdbx_component_id(String pdbx_component_id) {
		this.pdbx_component_id = pdbx_component_id;
	}
	public String getPdbx_ordinal() {
		return pdbx_ordinal;
	}
	public void setPdbx_ordinal(String pdbx_ordinal) {
		this.pdbx_ordinal = pdbx_ordinal;
	}
	
	
	
}
