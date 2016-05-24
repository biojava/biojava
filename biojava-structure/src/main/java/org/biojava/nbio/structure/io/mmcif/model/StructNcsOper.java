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
 * A class containing the _struct_ncs_oper data
 * 
 * <pre>
 *  _struct_ncs_oper.id 
 *  _struct_ncs_oper.code 
 *  _struct_ncs_oper.details 
 * 	_struct_ncs_oper.matrix[1][1] 
 *	_struct_ncs_oper.matrix[1][2] 
 *	_struct_ncs_oper.matrix[1][3] 
 *	_struct_ncs_oper.matrix[2][1] 
 *	_struct_ncs_oper.matrix[2][2] 
 *	_struct_ncs_oper.matrix[2][3] 
 *	_struct_ncs_oper.matrix[3][1] 
 *	_struct_ncs_oper.matrix[3][2] 
 *	_struct_ncs_oper.matrix[3][3] 
 *	_struct_ncs_oper.vector[1] 
 *	_struct_ncs_oper.vector[2] 
 *	_struct_ncs_oper.vector[3] 
 * </pre>
 * 
 * @author Jose Duarte
 */
public class StructNcsOper extends AbstractBean {

	private String id;
	private String code;
	private String details;
	
	@CIFLabel(label="matrix[1][1]")
	private String matrix11;
	
	@CIFLabel(label="matrix[1][2]")
	private String matrix12;
	
	@CIFLabel(label="matrix[1][3]")
	private String matrix13;
	
	@CIFLabel(label="matrix[2][1]")
	private String matrix21;
	
	@CIFLabel(label="matrix[2][2]")
	private String matrix22;
	
	@CIFLabel(label="matrix[2][3]")
	private String matrix23;
	
	@CIFLabel(label="matrix[3][1]")
	private String matrix31;
	
	@CIFLabel(label="matrix[3][2]")
	private String matrix32;
	
	@CIFLabel(label="matrix[3][3]")
	private String matrix33;
	
	@CIFLabel(label="vector[1]")
	private String vector1;
	
	@CIFLabel(label="vector[2]")
	private String vector2;
	
	@CIFLabel(label="vector[3]")
	private String vector3;
	
	
	public String getId() {
		return id;
	}

	public void setId(String id) {
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

	/**
	 * @return the matrix11
	 */
	public String getMatrix11() {
		return matrix11;
	}

	/**
	 * @param matrix11 the matrix11 to set
	 */
	public void setMatrix11(String matrix11) {
		this.matrix11 = matrix11;
	}

	/**
	 * @return the matrix12
	 */
	public String getMatrix12() {
		return matrix12;
	}

	/**
	 * @param matrix12 the matrix12 to set
	 */
	public void setMatrix12(String matrix12) {
		this.matrix12 = matrix12;
	}

	/**
	 * @return the matrix13
	 */
	public String getMatrix13() {
		return matrix13;
	}

	/**
	 * @param matrix13 the matrix13 to set
	 */
	public void setMatrix13(String matrix13) {
		this.matrix13 = matrix13;
	}

	/**
	 * @return the matrix21
	 */
	public String getMatrix21() {
		return matrix21;
	}

	/**
	 * @param matrix21 the matrix21 to set
	 */
	public void setMatrix21(String matrix21) {
		this.matrix21 = matrix21;
	}

	/**
	 * @return the matrix22
	 */
	public String getMatrix22() {
		return matrix22;
	}

	/**
	 * @param matrix22 the matrix22 to set
	 */
	public void setMatrix22(String matrix22) {
		this.matrix22 = matrix22;
	}

	/**
	 * @return the matrix23
	 */
	public String getMatrix23() {
		return matrix23;
	}

	/**
	 * @param matrix23 the matrix23 to set
	 */
	public void setMatrix23(String matrix23) {
		this.matrix23 = matrix23;
	}

	/**
	 * @return the matrix31
	 */
	public String getMatrix31() {
		return matrix31;
	}

	/**
	 * @param matrix31 the matrix31 to set
	 */
	public void setMatrix31(String matrix31) {
		this.matrix31 = matrix31;
	}

	/**
	 * @return the matrix32
	 */
	public String getMatrix32() {
		return matrix32;
	}

	/**
	 * @param matrix32 the matrix32 to set
	 */
	public void setMatrix32(String matrix32) {
		this.matrix32 = matrix32;
	}

	/**
	 * @return the matrix33
	 */
	public String getMatrix33() {
		return matrix33;
	}

	/**
	 * @param matrix33 the matrix33 to set
	 */
	public void setMatrix33(String matrix33) {
		this.matrix33 = matrix33;
	}

	/**
	 * @return the vector1
	 */
	public String getVector1() {
		return vector1;
	}

	/**
	 * @param vector1 the vector1 to set
	 */
	public void setVector1(String vector1) {
		this.vector1 = vector1;
	}

	/**
	 * @return the vector2
	 */
	public String getVector2() {
		return vector2;
	}

	/**
	 * @param vector2 the vector2 to set
	 */
	public void setVector2(String vector2) {
		this.vector2 = vector2;
	}

	/**
	 * @return the vector3
	 */
	public String getVector3() {
		return vector3;
	}

	/**
	 * @param vector3 the vector3 to set
	 */
	public void setVector3(String vector3) {
		this.vector3 = vector3;
	}

}
