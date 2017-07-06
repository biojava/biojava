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
