/*
 *                    PDB web development code
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
 *
 * Created on Mar 05, 2014
 * Created by Peter Rose
 *
 */

package org.biojava.nbio.structure.io.mmcif.model;
/**
 * A bean that stores data from the mmcif category _struct_conn
 * @author Peter Rose
 *
 */
public class StructConn extends AbstractBean
{
	private String id;
	private String conn_type_id;
	private String pdbx_PDB_id;
	private String ptnr1_label_asym_id;
	private String ptnr1_label_comp_id;
	private String ptnr1_label_seq_id;
	private String ptnr1_label_atom_id;
	private String pdbx_ptnr1_label_alt_id;
	private String pdbx_ptnr1_PDB_ins_code;
	private String pdbx_ptnr1_standard_comp_id;
	private String ptnr1_symmetry;
	private String ptnr2_label_asym_id;
	private String ptnr2_label_comp_id;
	private String ptnr2_label_seq_id;
	private String ptnr2_label_atom_id;
	private String pdbx_ptnr2_label_alt_id;
	private String pdbx_ptnr2_PDB_ins_code;
	private String ptnr1_auth_asym_id;
	private String ptnr1_auth_comp_id;
	private String ptnr1_auth_seq_id;
	private String ptnr2_auth_asym_id;
	private String ptnr2_auth_comp_id;
	private String ptnr2_auth_seq_id;
	private String ptnr2_symmetry;
	private String pdbx_ptnr3_label_atom_id;
	private String pdbx_ptnr3_label_seq_id;
	private String pdbx_ptnr3_label_comp_id;
	private String pdbx_ptnr3_label_asym_id;
	private String pdbx_ptnr3_label_alt_id;
	private String pdbx_ptnr3_PDB_ins_code;
	private String details;
	private String pdbx_dist_value;
	private String pdbx_value_order;
	private String pdbx_leaving_atom_flag;
	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}
	/**
	 * @param id the id to set
	 */
	public void setId(String id) {
		this.id = id;
	}
	/**
	 * @return the conn_type_id
	 */
	public String getConn_type_id() {
		return conn_type_id;
	}
	/**
	 * @param conn_type_id the conn_type_id to set
	 */
	public void setConn_type_id(String conn_type_id) {
		this.conn_type_id = conn_type_id;
	}
	/**
	 * @return the pdbx_PDB_id
	 */
	public String getPdbx_PDB_id() {
		return pdbx_PDB_id;
	}
	/**
	 * @param pdbx_PDB_id the pdbx_PDB_id to set
	 */
	public void setPdbx_PDB_id(String pdbx_PDB_id) {
		this.pdbx_PDB_id = pdbx_PDB_id;
	}
	/**
	 * @return the ptnr1_label_asym_id
	 */
	public String getPtnr1_label_asym_id() {
		return ptnr1_label_asym_id;
	}
	/**
	 * @param ptnr1_label_asym_id the ptnr1_label_asym_id to set
	 */
	public void setPtnr1_label_asym_id(String ptnr1_label_asym_id) {
		this.ptnr1_label_asym_id = ptnr1_label_asym_id;
	}
	/**
	 * @return the ptnr1_label_comp_id
	 */
	public String getPtnr1_label_comp_id() {
		return ptnr1_label_comp_id;
	}
	/**
	 * @param ptnr1_label_comp_id the ptnr1_label_comp_id to set
	 */
	public void setPtnr1_label_comp_id(String ptnr1_label_comp_id) {
		this.ptnr1_label_comp_id = ptnr1_label_comp_id;
	}
	/**
	 * @return the ptnr1_label_seq_id
	 */
	public String getPtnr1_label_seq_id() {
		return ptnr1_label_seq_id;
	}
	/**
	 * @param ptnr1_label_seq_id the ptnr1_label_seq_id to set
	 */
	public void setPtnr1_label_seq_id(String ptnr1_label_seq_id) {
		this.ptnr1_label_seq_id = ptnr1_label_seq_id;
	}
	/**
	 * @return the ptnr1_label_atom_id
	 */
	public String getPtnr1_label_atom_id() {
		return ptnr1_label_atom_id;
	}
	/**
	 * @param ptnr1_label_atom_id the ptnr1_label_atom_id to set
	 */
	public void setPtnr1_label_atom_id(String ptnr1_label_atom_id) {
		this.ptnr1_label_atom_id = ptnr1_label_atom_id;
	}
	/**
	 * @return the pdbx_ptnr1_label_alt_id
	 */
	public String getPdbx_ptnr1_label_alt_id() {
		return pdbx_ptnr1_label_alt_id;
	}
	/**
	 * @param pdbx_ptnr1_label_alt_id the pdbx_ptnr1_label_alt_id to set
	 */
	public void setPdbx_ptnr1_label_alt_id(String pdbx_ptnr1_label_alt_id) {
		this.pdbx_ptnr1_label_alt_id = pdbx_ptnr1_label_alt_id;
	}
	/**
	 * @return the pdbx_ptnr1_PDB_ins_code
	 */
	public String getPdbx_ptnr1_PDB_ins_code() {
		return pdbx_ptnr1_PDB_ins_code;
	}
	/**
	 * @param pdbx_ptnr1_PDB_ins_code the pdbx_ptnr1_PDB_ins_code to set
	 */
	public void setPdbx_ptnr1_PDB_ins_code(String pdbx_ptnr1_PDB_ins_code) {
		this.pdbx_ptnr1_PDB_ins_code = pdbx_ptnr1_PDB_ins_code;
	}
	/**
	 * @return the pdbx_ptnr1_standard_comp_id
	 */
	public String getPdbx_ptnr1_standard_comp_id() {
		return pdbx_ptnr1_standard_comp_id;
	}
	/**
	 * @param pdbx_ptnr1_standard_comp_id the pdbx_ptnr1_standard_comp_id to set
	 */
	public void setPdbx_ptnr1_standard_comp_id(String pdbx_ptnr1_standard_comp_id) {
		this.pdbx_ptnr1_standard_comp_id = pdbx_ptnr1_standard_comp_id;
	}
	/**
	 * @return the ptnr1_symmetry
	 */
	public String getPtnr1_symmetry() {
		return ptnr1_symmetry;
	}
	/**
	 * @param ptnr1_symmetry the ptnr1_symmetry to set
	 */
	public void setPtnr1_symmetry(String ptnr1_symmetry) {
		this.ptnr1_symmetry = ptnr1_symmetry;
	}
	/**
	 * @return the ptnr2_label_asym_id
	 */
	public String getPtnr2_label_asym_id() {
		return ptnr2_label_asym_id;
	}
	/**
	 * @param ptnr2_label_asym_id the ptnr2_label_asym_id to set
	 */
	public void setPtnr2_label_asym_id(String ptnr2_label_asym_id) {
		this.ptnr2_label_asym_id = ptnr2_label_asym_id;
	}
	/**
	 * @return the ptnr2_label_comp_id
	 */
	public String getPtnr2_label_comp_id() {
		return ptnr2_label_comp_id;
	}
	/**
	 * @param ptnr2_label_comp_id the ptnr2_label_comp_id to set
	 */
	public void setPtnr2_label_comp_id(String ptnr2_label_comp_id) {
		this.ptnr2_label_comp_id = ptnr2_label_comp_id;
	}
	/**
	 * @return the ptnr2_label_seq_id
	 */
	public String getPtnr2_label_seq_id() {
		return ptnr2_label_seq_id;
	}
	/**
	 * @param ptnr2_label_seq_id the ptnr2_label_seq_id to set
	 */
	public void setPtnr2_label_seq_id(String ptnr2_label_seq_id) {
		this.ptnr2_label_seq_id = ptnr2_label_seq_id;
	}
	/**
	 * @return the ptnr2_label_atom_id
	 */
	public String getPtnr2_label_atom_id() {
		return ptnr2_label_atom_id;
	}
	/**
	 * @param ptnr2_label_atom_id the ptnr2_label_atom_id to set
	 */
	public void setPtnr2_label_atom_id(String ptnr2_label_atom_id) {
		this.ptnr2_label_atom_id = ptnr2_label_atom_id;
	}
	/**
	 * @return the pdbx_ptnr2_label_alt_id
	 */
	public String getPdbx_ptnr2_label_alt_id() {
		return pdbx_ptnr2_label_alt_id;
	}
	/**
	 * @param pdbx_ptnr2_label_alt_id the pdbx_ptnr2_label_alt_id to set
	 */
	public void setPdbx_ptnr2_label_alt_id(String pdbx_ptnr2_label_alt_id) {
		this.pdbx_ptnr2_label_alt_id = pdbx_ptnr2_label_alt_id;
	}
	/**
	 * @return the pdbx_ptnr2_PDB_ins_code
	 */
	public String getPdbx_ptnr2_PDB_ins_code() {
		return pdbx_ptnr2_PDB_ins_code;
	}
	/**
	 * @param pdbx_ptnr2_PDB_ins_code the pdbx_ptnr2_PDB_ins_code to set
	 */
	public void setPdbx_ptnr2_PDB_ins_code(String pdbx_ptnr2_PDB_ins_code) {
		this.pdbx_ptnr2_PDB_ins_code = pdbx_ptnr2_PDB_ins_code;
	}
	/**
	 * @return the ptnr1_auth_asym_id
	 */
	public String getPtnr1_auth_asym_id() {
		return ptnr1_auth_asym_id;
	}
	/**
	 * @param ptnr1_auth_asym_id the ptnr1_auth_asym_id to set
	 */
	public void setPtnr1_auth_asym_id(String ptnr1_auth_asym_id) {
		this.ptnr1_auth_asym_id = ptnr1_auth_asym_id;
	}
	/**
	 * @return the ptnr1_auth_comp_id
	 */
	public String getPtnr1_auth_comp_id() {
		return ptnr1_auth_comp_id;
	}
	/**
	 * @param ptnr1_auth_comp_id the ptnr1_auth_comp_id to set
	 */
	public void setPtnr1_auth_comp_id(String ptnr1_auth_comp_id) {
		this.ptnr1_auth_comp_id = ptnr1_auth_comp_id;
	}
	/**
	 * @return the ptnr1_auth_seq_id
	 */
	public String getPtnr1_auth_seq_id() {
		return ptnr1_auth_seq_id;
	}
	/**
	 * @param ptnr1_auth_seq_id the ptnr1_auth_seq_id to set
	 */
	public void setPtnr1_auth_seq_id(String ptnr1_auth_seq_id) {
		this.ptnr1_auth_seq_id = ptnr1_auth_seq_id;
	}
	/**
	 * @return the ptnr2_auth_asym_id
	 */
	public String getPtnr2_auth_asym_id() {
		return ptnr2_auth_asym_id;
	}
	/**
	 * @param ptnr2_auth_asym_id the ptnr2_auth_asym_id to set
	 */
	public void setPtnr2_auth_asym_id(String ptnr2_auth_asym_id) {
		this.ptnr2_auth_asym_id = ptnr2_auth_asym_id;
	}
	/**
	 * @return the ptnr2_auth_comp_id
	 */
	public String getPtnr2_auth_comp_id() {
		return ptnr2_auth_comp_id;
	}
	/**
	 * @param ptnr2_auth_comp_id the ptnr2_auth_comp_id to set
	 */
	public void setPtnr2_auth_comp_id(String ptnr2_auth_comp_id) {
		this.ptnr2_auth_comp_id = ptnr2_auth_comp_id;
	}
	/**
	 * @return the ptnr2_auth_seq_id
	 */
	public String getPtnr2_auth_seq_id() {
		return ptnr2_auth_seq_id;
	}
	/**
	 * @param ptnr2_auth_seq_id the ptnr2_auth_seq_id to set
	 */
	public void setPtnr2_auth_seq_id(String ptnr2_auth_seq_id) {
		this.ptnr2_auth_seq_id = ptnr2_auth_seq_id;
	}
	/**
	 * @return the ptnr2_symmetry
	 */
	public String getPtnr2_symmetry() {
		return ptnr2_symmetry;
	}
	/**
	 * @param ptnr2_symmetry the ptnr2_symmetry to set
	 */
	public void setPtnr2_symmetry(String ptnr2_symmetry) {
		this.ptnr2_symmetry = ptnr2_symmetry;
	}
	/**
	 * @return the pdbx_ptnr3_label_atom_id
	 */
	public String getPdbx_ptnr3_label_atom_id() {
		return pdbx_ptnr3_label_atom_id;
	}
	/**
	 * @param pdbx_ptnr3_label_atom_id the pdbx_ptnr3_label_atom_id to set
	 */
	public void setPdbx_ptnr3_label_atom_id(String pdbx_ptnr3_label_atom_id) {
		this.pdbx_ptnr3_label_atom_id = pdbx_ptnr3_label_atom_id;
	}
	/**
	 * @return the pdbx_ptnr3_label_seq_id
	 */
	public String getPdbx_ptnr3_label_seq_id() {
		return pdbx_ptnr3_label_seq_id;
	}
	/**
	 * @param pdbx_ptnr3_label_seq_id the pdbx_ptnr3_label_seq_id to set
	 */
	public void setPdbx_ptnr3_label_seq_id(String pdbx_ptnr3_label_seq_id) {
		this.pdbx_ptnr3_label_seq_id = pdbx_ptnr3_label_seq_id;
	}
	/**
	 * @return the pdbx_ptnr3_label_comp_id
	 */
	public String getPdbx_ptnr3_label_comp_id() {
		return pdbx_ptnr3_label_comp_id;
	}
	/**
	 * @param pdbx_ptnr3_label_comp_id the pdbx_ptnr3_label_comp_id to set
	 */
	public void setPdbx_ptnr3_label_comp_id(String pdbx_ptnr3_label_comp_id) {
		this.pdbx_ptnr3_label_comp_id = pdbx_ptnr3_label_comp_id;
	}
	/**
	 * @return the pdbx_ptnr3_label_asym_id
	 */
	public String getPdbx_ptnr3_label_asym_id() {
		return pdbx_ptnr3_label_asym_id;
	}
	/**
	 * @param pdbx_ptnr3_label_asym_id the pdbx_ptnr3_label_asym_id to set
	 */
	public void setPdbx_ptnr3_label_asym_id(String pdbx_ptnr3_label_asym_id) {
		this.pdbx_ptnr3_label_asym_id = pdbx_ptnr3_label_asym_id;
	}
	/**
	 * @return the pdbx_ptnr3_label_alt_id
	 */
	public String getPdbx_ptnr3_label_alt_id() {
		return pdbx_ptnr3_label_alt_id;
	}
	/**
	 * @param pdbx_ptnr3_label_alt_id the pdbx_ptnr3_label_alt_id to set
	 */
	public void setPdbx_ptnr3_label_alt_id(String pdbx_ptnr3_label_alt_id) {
		this.pdbx_ptnr3_label_alt_id = pdbx_ptnr3_label_alt_id;
	}
	/**
	 * @return the pdbx_ptnr3_PDB_ins_code
	 */
	public String getPdbx_ptnr3_PDB_ins_code() {
		return pdbx_ptnr3_PDB_ins_code;
	}
	/**
	 * @param pdbx_ptnr3_PDB_ins_code the pdbx_ptnr3_PDB_ins_code to set
	 */
	public void setPdbx_ptnr3_PDB_ins_code(String pdbx_ptnr3_PDB_ins_code) {
		this.pdbx_ptnr3_PDB_ins_code = pdbx_ptnr3_PDB_ins_code;
	}
	/**
	 * @return the details
	 */
	public String getDetails() {
		return details;
	}
	/**
	 * @param details the details to set
	 */
	public void setDetails(String details) {
		this.details = details;
	}
	/**
	 * @return the pdbx_dist_value
	 */
	public String getPdbx_dist_value() {
		return pdbx_dist_value;
	}
	/**
	 * @param pdbx_dist_value the pdbx_dist_value to set
	 */
	public void setPdbx_dist_value(String pdbx_dist_value) {
		this.pdbx_dist_value = pdbx_dist_value;
	}
	/**
	 * @return the pdbx_value_order
	 */
	public String getPdbx_value_order() {
		return pdbx_value_order;
	}
	/**
	 * @param pdbx_value_order the pdbx_value_order to set
	 */
	public void setPdbx_value_order(String pdbx_value_order) {
		this.pdbx_value_order = pdbx_value_order;
	}

	public String getPdbx_leaving_atom_flag() {
		return pdbx_leaving_atom_flag;
	}

	public void setPdbx_leaving_atom_flag(String pdbx_leaving_atom_flag) {
		this.pdbx_leaving_atom_flag = pdbx_leaving_atom_flag;
	}
}
