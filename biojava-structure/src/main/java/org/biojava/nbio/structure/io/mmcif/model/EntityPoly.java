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
 * created at Jun 1, 2008
 */
package org.biojava.nbio.structure.io.mmcif.model;


/** 
 * Container for _entity_poly records
 *
 *
 * @since 5.0
 * @author Jose Duarte
 */
public class EntityPoly extends AbstractBean{
	String entity_id;
	String type;
	String nstd_chirality;
	String nstd_linkage;
	String nstd_monomer;
	String type_details;
	String pdbx_seq_one_letter_code;
	String pdbx_seq_one_letter_code_can;
	String pdbx_strand_id;
	String pdbx_target_identifier;
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	/**
	 * @return the type
	 */
	public String getType() {
		return type;
	}
	/**
	 * @param type the type to set
	 */
	public void setType(String type) {
		this.type = type;
	}
	/**
	 * @return the nstd_chirality
	 */
	public String getNstd_chirality() {
		return nstd_chirality;
	}
	/**
	 * @param nstd_chirality the nstd_chirality to set
	 */
	public void setNstd_chirality(String nstd_chirality) {
		this.nstd_chirality = nstd_chirality;
	}
	/**
	 * @return the nstd_linkage
	 */
	public String getNstd_linkage() {
		return nstd_linkage;
	}
	/**
	 * @param nstd_linkage the nstd_linkage to set
	 */
	public void setNstd_linkage(String nstd_linkage) {
		this.nstd_linkage = nstd_linkage;
	}
	/**
	 * @return the nstd_monomer
	 */
	public String getNstd_monomer() {
		return nstd_monomer;
	}
	/**
	 * @param nstd_monomer the nstd_monomer to set
	 */
	public void setNstd_monomer(String nstd_monomer) {
		this.nstd_monomer = nstd_monomer;
	}
	/**
	 * @return the type_details
	 */
	public String getType_details() {
		return type_details;
	}
	/**
	 * @param type_details the type_details to set
	 */
	public void setType_details(String type_details) {
		this.type_details = type_details;
	}
	/**
	 * @return the pdbx_seq_one_letter_code
	 */
	public String getPdbx_seq_one_letter_code() {
		return pdbx_seq_one_letter_code;
	}
	/**
	 * @param pdbx_seq_one_letter_code the pdbx_seq_one_letter_code to set
	 */
	public void setPdbx_seq_one_letter_code(String pdbx_seq_one_letter_code) {
		this.pdbx_seq_one_letter_code = pdbx_seq_one_letter_code;
	}
	/**
	 * @return the pdbx_seq_one_letter_code_can
	 */
	public String getPdbx_seq_one_letter_code_can() {
		return pdbx_seq_one_letter_code_can;
	}
	/**
	 * @param pdbx_seq_one_letter_code_can the pdbx_seq_one_letter_code_can to set
	 */
	public void setPdbx_seq_one_letter_code_can(String pdbx_seq_one_letter_code_can) {
		this.pdbx_seq_one_letter_code_can = pdbx_seq_one_letter_code_can;
	}
	/**
	 * @return the pdbx_strand_id
	 */
	public String getPdbx_strand_id() {
		return pdbx_strand_id;
	}
	/**
	 * @param pdbx_strand_id the pdbx_strand_id to set
	 */
	public void setPdbx_strand_id(String pdbx_strand_id) {
		this.pdbx_strand_id = pdbx_strand_id;
	}
	public String getPdbx_target_identifier() {
		return pdbx_target_identifier;
	}
	public void setPdbx_target_identifier(String pdbx_target_identifier) {
		this.pdbx_target_identifier = pdbx_target_identifier;
	}

}
