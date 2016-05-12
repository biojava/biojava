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


import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import java.io.Serializable;



@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class PdbxStructAssembly implements Serializable{

	/**
	 *
	 */
	private static final long serialVersionUID = 3104504686693887219L;

	String id;
	String details;
	String method_details;
	String oligomeric_details;
	String oligomeric_count ;
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
	public String getMethod_details() {
		return method_details;
	}
	public void setMethod_details(String method_details) {
		this.method_details = method_details;
	}
	public String getOligomeric_details() {
		return oligomeric_details;
	}
	public void setOligomeric_details(String oligomeric_details) {
		this.oligomeric_details = oligomeric_details;
	}
	public String getOligomeric_count() {
		return oligomeric_count;
	}
	public void setOligomeric_count(String oligomeric_count) {
		this.oligomeric_count = oligomeric_count;
	}
	@Override
	public String toString() {
		return "PdbxStructAssembly [id=" + id + ", details=" + details
				+ ", method_details=" + method_details
				+ ", oligomeric_details=" + oligomeric_details
				+ ", oligomeric_count=" + oligomeric_count + "]";
	}





}
