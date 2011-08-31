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
 * created at Apr 26, 2008
 */
package org.biojava.bio.structure.io.mmcif.model;

/** a bean to contain the data of the _struct lines
 *
 * @author Andreas Prlic
 *
 */
public class Struct {
	String entry_id;
	String title;
	String pdbx_descriptor;
	String pdbx_model_details;
	String pdbx_model_type_details;
	String pdbx_CASP_flag;

	public String toString(){
		return "entry_id:" +entry_id + " title:" + title + " pdbx_descriptor:" +pdbx_descriptor + " pdbx_model_details:"+pdbx_model_details;
	}

	public String getEntry_id() {
		return entry_id;
	}
	public void setEntry_id(String entry_id) {
		this.entry_id = entry_id;
	}
	public String getTitle() {
		return title;
	}
	public void setTitle(String title) {
		this.title = title;
	}
	public String getPdbx_descriptor() {
		return pdbx_descriptor;
	}
	public void setPdbx_descriptor(String pdbx_descriptor) {
		this.pdbx_descriptor = pdbx_descriptor;
	}
	public String getPdbx_model_details() {
		return pdbx_model_details;
	}
	public void setPdbx_model_details(String pdbx_model_details) {
		this.pdbx_model_details = pdbx_model_details;
	}

	public String getPdbx_model_type_details() {
		return pdbx_model_type_details;
	}

	public void setPdbx_model_type_details(String pdbx_model_type_details) {
		this.pdbx_model_type_details = pdbx_model_type_details;
	}

	public String getPdbx_CASP_flag() {
		return pdbx_CASP_flag;
	}

	public void setPdbx_CASP_flag(String pdbx_CASP_flag) {
		this.pdbx_CASP_flag = pdbx_CASP_flag;
	}




}
