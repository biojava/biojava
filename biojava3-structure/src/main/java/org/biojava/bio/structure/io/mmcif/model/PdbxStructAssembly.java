package org.biojava.bio.structure.io.mmcif.model;


import java.io.Serializable;


import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;



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
