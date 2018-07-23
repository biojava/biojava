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
 * Bean to hold data for _pdbx_audit_revision_history mmCIF category.
 * 
 * @author Peter Rose
 * @since 5.0
 */
public class PdbxAuditRevisionHistory extends AbstractBean {
	private String ordinal;
	private String data_content_type;
	private String major_revision;
	private String minor_revision;
	private String revision_date;
	
	public String getOrdinal() {
		return ordinal;
	}
	public void setOrdinal(String ordinal) {
		this.ordinal = ordinal;
	}
	public String getData_content_type() {
		return data_content_type;
	}
	public void setData_content_type(String data_content_type) {
		this.data_content_type = data_content_type;
	}
	public String getMajor_revision() {
		return major_revision;
	}
	public void setMajor_revision(String major_revision) {
		this.major_revision = major_revision;
	}
	public String getMinor_revision() {
		return minor_revision;
	}
	public void setMinor_revision(String minor_revision) {
		this.minor_revision = minor_revision;
	}
	public String getRevision_date() {
		return revision_date;
	}
	public void setRevision_date(String revision_date) {
		this.revision_date = revision_date;
	}
}
