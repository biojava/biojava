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
 * Bean to hold data for _pdbx_database_status mmCIF category.
 * 
 * @author Peter Rose
 * @since 5.0
 */
public class PdbxDatabaseStatus extends AbstractBean {
	private String status_code;
	private String entry_id;
	private String recvd_initial_deposition_date;
	private String deposit_site;
	private String process_site;
	private String SG_entry;
	private String pdb_format_compatible;
	private String status_code_mr;
	private String status_code_sf;
	private String status_code_cs;
	
	public String getStatus_code() {
		return status_code;
	}
	public void setStatus_code(String status_code) {
		this.status_code = status_code;
	}
	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	public String getRecvd_initial_deposition_date() {
		return recvd_initial_deposition_date;
	}
	public void setRecvd_initial_deposition_date(String recvd_initial_deposition_date) {
		this.recvd_initial_deposition_date = recvd_initial_deposition_date;
	}
	public String getDeposit_site() {
		return deposit_site;
	}
	public void setDeposit_site(String deposit_site) {
		this.deposit_site = deposit_site;
	}
	public String getProcess_site() {
		return process_site;
	}
	public void setProcess_site(String process_site) {
		this.process_site = process_site;
	}
	public String getSG_entry() {
		return SG_entry;
	}
	public void setSG_entry(String sG_entry) {
		SG_entry = sG_entry;
	}
	public String getPdb_format_compatible() {
		return pdb_format_compatible;
	}
	public void setPdb_format_compatible(String pdb_format_compatible) {
		this.pdb_format_compatible = pdb_format_compatible;
	}
	public String getStatus_code_mr() {
		return status_code_mr;
	}
	public void setStatus_code_mr(String status_code_mr) {
		this.status_code_mr = status_code_mr;
	}
	public String getStatus_code_sf() {
		return status_code_sf;
	}
	public void setStatus_code_sf(String status_code_sf) {
		this.status_code_sf = status_code_sf;
	}
	public String getStatus_code_cs() {
		return status_code_cs;
	}
	public void setStatus_code_cs(String status_code_cs) {
		this.status_code_cs = status_code_cs;
	}  
}
