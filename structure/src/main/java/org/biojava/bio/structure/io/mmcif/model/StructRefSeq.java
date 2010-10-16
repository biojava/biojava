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

public class StructRefSeq extends AbstractBean{
	String align_id;
	String ref_id;
	String pdbx_PDB_id_code;
	String pdbx_strand_id;
	String seq_align_beg;
	String pdbx_seq_align_beg_ins_code;
	String seq_align_end;
	String pdbx_seq_align_end_ins_code;
	String pdbx_db_accession;
	String db_align_beg;
	String pdbx_db_align_beg_ins_code;
	String db_align_end;
	String pdbx_db_align_end_ins_code;
	String pdbx_auth_seq_align_beg;
	String pdbx_auth_seq_align_end;
	String details;

	public StructRefSeq(){
		super();
		pdbx_db_align_beg_ins_code = "?";
		pdbx_db_align_end_ins_code = "?";

	}

	public String getAlign_id() {
		return align_id;
	}
	public void setAlign_id(String align_id) {
		this.align_id = align_id;
	}
	public String getRef_id() {
		return ref_id;
	}
	public void setRef_id(String ref_id) {
		this.ref_id = ref_id;
	}
	public String getPdbx_PDB_id_code() {
		return pdbx_PDB_id_code;
	}
	public void setPdbx_PDB_id_code(String pdbx_PDB_id_code) {
		this.pdbx_PDB_id_code = pdbx_PDB_id_code;
	}
	public String getPdbx_strand_id() {
		return pdbx_strand_id;
	}
	public void setPdbx_strand_id(String pdbx_strand_id) {
		this.pdbx_strand_id = pdbx_strand_id;
	}
	public String getSeq_align_beg() {
		return seq_align_beg;
	}
	public void setSeq_align_beg(String seq_align_beg) {
		this.seq_align_beg = seq_align_beg;
	}
	public String getPdbx_seq_align_beg_ins_code() {
		return pdbx_seq_align_beg_ins_code;
	}
	public void setPdbx_seq_align_beg_ins_code(String pdbx_seq_align_beg_ins_code) {
		this.pdbx_seq_align_beg_ins_code = pdbx_seq_align_beg_ins_code;
	}
	public String getSeq_align_end() {
		return seq_align_end;
	}
	public void setSeq_align_end(String seq_align_end) {
		this.seq_align_end = seq_align_end;
	}
	public String getPdbx_seq_align_end_ins_code() {
		return pdbx_seq_align_end_ins_code;
	}
	public void setPdbx_seq_align_end_ins_code(String pdbx_seq_align_end_ins_code) {
		this.pdbx_seq_align_end_ins_code = pdbx_seq_align_end_ins_code;
	}
	public String getPdbx_db_accession() {
		return pdbx_db_accession;
	}
	public void setPdbx_db_accession(String pdbx_db_accession) {
		this.pdbx_db_accession = pdbx_db_accession;
	}
	public String getDb_align_beg() {
		return db_align_beg;
	}
	public void setDb_align_beg(String db_align_beg) {
		this.db_align_beg = db_align_beg;
	}
	public String getPdbx_db_align_beg_ins_code() {
		return pdbx_db_align_beg_ins_code;
	}
	public void setPdbx_db_align_beg_ins_code(String pdbx_db_align_beg_ins_code) {
		this.pdbx_db_align_beg_ins_code = pdbx_db_align_beg_ins_code;
	}
	public String getDb_align_end() {
		return db_align_end;
	}
	public void setDb_align_end(String db_align_end) {
		this.db_align_end = db_align_end;
	}
	public String getPdbx_db_align_end_ins_code() {
		return pdbx_db_align_end_ins_code;
	}
	public void setPdbx_db_align_end_ins_code(String pdbx_db_align_end_ins_code) {
		this.pdbx_db_align_end_ins_code = pdbx_db_align_end_ins_code;
	}
	public String getPdbx_auth_seq_align_beg() {
		return pdbx_auth_seq_align_beg;
	}
	public void setPdbx_auth_seq_align_beg(String pdbx_auth_seq_align_beg) {
		this.pdbx_auth_seq_align_beg = pdbx_auth_seq_align_beg;
	}
	public String getPdbx_auth_seq_align_end() {
		return pdbx_auth_seq_align_end;
	}
	public void setPdbx_auth_seq_align_end(String pdbx_auth_seq_align_end) {
		this.pdbx_auth_seq_align_end = pdbx_auth_seq_align_end;
	}
	public String getDetails() {
		return details;
	}
	public void setDetails(String details) {
		this.details = details;
	}

}
