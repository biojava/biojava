package org.biojava.nbio.structure.io.mmcif.model;

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
 * Created by andreas on 9/11/15.
 */

/** A class to store sequence mismatch annotations
 *
 */
public class StructRefSeqDif {

	String align_id;
	String pdbx_pdb_id_code;
	String mon_id;
	String pdbx_pdb_strand_id;
	Integer seq_num;
	String pdbx_pdb_ins_code;
	String pdbx_seq_db_name;
	String pdbx_seq_db_accession_code;
	String db_mon_id;
	String pdbx_seq_db_seq_num;
	String details;
	String pdbx_auth_seq_num;
	Integer pdbx_ordinal;

	public String getAlign_id() {
		return align_id;
	}

	public void setAlign_id(String align_id) {
		this.align_id = align_id;
	}

	public String getPdbx_pdb_id_code() {
		return pdbx_pdb_id_code;
	}

	public void setPdbx_pdb_id_code(String pdbx_pdb_id_code) {
		this.pdbx_pdb_id_code = pdbx_pdb_id_code;
	}

	public String getMon_id() {
		return mon_id;
	}

	public void setMon_id(String mon_id) {
		this.mon_id = mon_id;
	}

	public String getPdbx_pdb_strand_id() {
		return pdbx_pdb_strand_id;
	}

	public void setPdbx_pdb_strand_id(String pdbx_pdb_strand_id) {
		this.pdbx_pdb_strand_id = pdbx_pdb_strand_id;
	}

	public Integer getSeq_num() {
		return seq_num;
	}

	public void setSeq_num(Integer seq_num) {
		this.seq_num = seq_num;
	}

	public String getPdbx_pdb_ins_code() {
		return pdbx_pdb_ins_code;
	}

	public void setPdbx_pdb_ins_code(String pdbx_pdb_ins_code) {
		this.pdbx_pdb_ins_code = pdbx_pdb_ins_code;
	}

	public String getPdbx_seq_db_name() {
		return pdbx_seq_db_name;
	}

	public void setPdbx_seq_db_name(String pdbx_seq_db_name) {
		this.pdbx_seq_db_name = pdbx_seq_db_name;
	}

	public String getPdbx_seq_db_accession_code() {
		return pdbx_seq_db_accession_code;
	}

	public void setPdbx_seq_db_accession_code(String pdbx_seq_db_accession_code) {
		this.pdbx_seq_db_accession_code = pdbx_seq_db_accession_code;
	}

	public String getDb_mon_id() {
		return db_mon_id;
	}

	public void setDb_mon_id(String db_mon_id) {
		this.db_mon_id = db_mon_id;
	}

	public String getPdbx_seq_db_seq_num() {
		return pdbx_seq_db_seq_num;
	}

	public void setPdbx_seq_db_seq_num(String pdbx_seq_db_seq_num) {
		this.pdbx_seq_db_seq_num = pdbx_seq_db_seq_num;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getPdbx_auth_seq_num() {
		return pdbx_auth_seq_num;
	}

	public void setPdbx_auth_seq_num(String pdbx_auth_seq_num) {
		this.pdbx_auth_seq_num = pdbx_auth_seq_num;
	}

	public Integer getPdbx_ordinal() {
		return pdbx_ordinal;
	}

	public void setPdbx_ordinal(Integer pdbx_ordinal) {
		this.pdbx_ordinal = pdbx_ordinal;
	}
}
