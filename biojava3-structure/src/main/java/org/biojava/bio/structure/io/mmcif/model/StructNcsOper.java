package org.biojava.bio.structure.io.mmcif.model;

import javax.vecmath.Matrix4d;

/**
 * A class containing the _struct_ncs_oper data
 * @author duarte_j
 *
 */
public class StructNcsOper extends AbstractBean {

	private int id; 
	private String code; 
	private String details; 
	private Matrix4d operator;
	
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public String getCode() {
		return code;
	}
	
	public void setCode(String code) {
		this.code = code;
	}
	
	public String getDetails() {
		return details;
	}
	
	public void setDetails(String details) {
		this.details = details;
	}
	
	public Matrix4d getOperator() {
		return operator;
	}
	
	public void setOperator(Matrix4d operator) {
		this.operator = operator;
	}
}
