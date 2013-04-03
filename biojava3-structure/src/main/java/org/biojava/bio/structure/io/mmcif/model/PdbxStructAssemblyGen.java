package org.biojava.bio.structure.io.mmcif.model;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class PdbxStructAssemblyGen implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 6739568389242514332L;
	String assembly_id; 
	String oper_expression; 
	String asym_id_list;
	
	
	public String getAssembly_id() {
		return assembly_id;
	}
	public void setAssembly_id(String assembly_id) {
		this.assembly_id = assembly_id;
	}
	public String getOper_expression() {
		return oper_expression;
	}
	public void setOper_expression(String oper_expression) {
		this.oper_expression = oper_expression;
	}
	public String getAsym_id_list() {
		return asym_id_list;
	}
	public void setAsym_id_list(String asym_id_list) {
		this.asym_id_list = asym_id_list;
	}
	@Override
	public String toString() {
		return "PdbxStructAssemblyGen [assembly_id=" + assembly_id
				+ ", oper_expression=" + oper_expression + ", asym_id_list="
				+ asym_id_list + "]";
	} 
	
	

	
	
}
