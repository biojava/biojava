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

public class Exptl extends AbstractBean{
	String entry_id;
	String method;
	String crystals_number;
	String absorpt_coefficient_mu;
	String absorpt_correction_T_max;
	String absorpt_correction_T_min ;
	String absorpt_correction_type ;
	String absorpt_process_details ;
	String details;
	String method_details;

	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	public String getMethod() {
		return method;
	}
	public void setMethod(String method) {
		this.method = method;
	}
	public String getCrystals_number() {
		return crystals_number;
	}
	public void setCrystals_number(String crystals_number) {
		this.crystals_number = crystals_number;
	}
	public String getAbsorpt_coefficient_mu() {
		return absorpt_coefficient_mu;
	}
	public void setAbsorpt_coefficient_mu(String absorpt_coefficient_mu) {
		this.absorpt_coefficient_mu = absorpt_coefficient_mu;
	}
	public String getAbsorpt_correction_T_max() {
		return absorpt_correction_T_max;
	}
	public void setAbsorpt_correction_T_max(String absorpt_correction_T_max) {
		this.absorpt_correction_T_max = absorpt_correction_T_max;
	}
	public String getAbsorpt_correction_T_min() {
		return absorpt_correction_T_min;
	}
	public void setAbsorpt_correction_T_min(String absorpt_correction_T_min) {
		this.absorpt_correction_T_min = absorpt_correction_T_min;
	}
	public String getAbsorpt_correction_type() {
		return absorpt_correction_type;
	}
	public void setAbsorpt_correction_type(String absorpt_correction_type) {
		this.absorpt_correction_type = absorpt_correction_type;
	}
	public String getAbsorpt_process_details() {
		return absorpt_process_details;
	}
	public void setAbsorpt_process_details(String absorpt_process_details) {
		this.absorpt_process_details = absorpt_process_details;
	}
	public String getDetails() {
		return details;
	}
	public void setDetails(String details) {
		this.details = details;
	}
	public String getMethod_details() {
		return method_details;
	}
	public void setMethod_details(String method_details) {
		this.method_details = method_details;
	}


}
