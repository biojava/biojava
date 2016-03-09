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
 */
package org.biojava.nbio.structure.io.mmcif.model;

/**
 * Created by Matt on 11/1/2015.
 */
public class StructSite {
	String id;
	String details;
	String pdbx_evidence_code;
	String pdbx_auth_asym_id;
	String pdbx_auth_comp_id;
	String pdbx_auth_seq_id;
	String pdbx_num_residues;

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getPdbx_evidence_code() {
		return pdbx_evidence_code;
	}

	public void setPdbx_evidence_code(String pdbx_evidence_code) {
		this.pdbx_evidence_code = pdbx_evidence_code;
	}

	/**
	 * @return the pdbx_auth_asym_id
	 */
	public String getPdbx_auth_asym_id() {
		return pdbx_auth_asym_id;
	}

	/**
	 * @param pdbx_auth_asym_id the pdbx_auth_asym_id to set
	 */
	public void setPdbx_auth_asym_id(String pdbx_auth_asym_id) {
		this.pdbx_auth_asym_id = pdbx_auth_asym_id;
	}

	/**
	 * @return the pdbx_auth_comp_id
	 */
	public String getPdbx_auth_comp_id() {
		return pdbx_auth_comp_id;
	}

	/**
	 * @param pdbx_auth_comp_id the pdbx_auth_comp_id to set
	 */
	public void setPdbx_auth_comp_id(String pdbx_auth_comp_id) {
		this.pdbx_auth_comp_id = pdbx_auth_comp_id;
	}

	/**
	 * @return the pdbx_auth_seq_id
	 */
	public String getPdbx_auth_seq_id() {
		return pdbx_auth_seq_id;
	}

	/**
	 * @param pdbx_auth_seq_id the pdbx_auth_seq_id to set
	 */
	public void setPdbx_auth_seq_id(String pdbx_auth_seq_id) {
		this.pdbx_auth_seq_id = pdbx_auth_seq_id;
	}

	/**
	 * @return the pdbx_num_residues
	 */
	public String getPdbx_num_residues() {
		return pdbx_num_residues;
	}

	/**
	 * @param pdbx_num_residues the pdbx_num_residues to set
	 */
	public void setPdbx_num_residues(String pdbx_num_residues) {
		this.pdbx_num_residues = pdbx_num_residues;
	}
}
