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
 * created at May 31, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

/** A class to containt the _struct_ref field data
 * 
 * @author Andreas Prlic
 *
 */
public class StructRef extends AbstractBean {
	String id;
	String db_name;
	String db_code;
	String entity_id;
	String pdbx_db_accession;
	String pdbx_align_begin;
	String pdbx_seq_one_letter_code;
	String biol_id;
	public String getBiol_id() {
		return biol_id;
	}
	public void setBiol_id(String biol_id) {
		this.biol_id = biol_id;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getDb_name() {
		return db_name;
	}
	public void setDb_name(String db_name) {
		this.db_name = db_name;
	}
	public String getDb_code() {
		return db_code;
	}
	public void setDb_code(String db_code) {
		this.db_code = db_code;
	}
	public String getEntity_id() {
		return entity_id;
	}
	public void setEntity_id(String entity_id) {
		this.entity_id = entity_id;
	}
	public String getPdbx_db_accession() {
		return pdbx_db_accession;
	}
	public void setPdbx_db_accession(String pdbx_db_accession) {
		this.pdbx_db_accession = pdbx_db_accession;
	}
	public String getPdbx_align_begin() {
		return pdbx_align_begin;
	}
	public void setPdbx_align_begin(String pdbx_align_begin) {
		this.pdbx_align_begin = pdbx_align_begin;
	}
	public String getPdbx_seq_one_letter_code() {
		return pdbx_seq_one_letter_code;
	}
	public void setPdbx_seq_one_letter_code(String pdbx_seq_one_letter_code) {
		this.pdbx_seq_one_letter_code = pdbx_seq_one_letter_code;
	}
	
	

}
