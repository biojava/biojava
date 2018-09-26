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
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import java.io.Serializable;

/**
 * The bean for pdbx_struct_oper_list category
 * <pre>
 * _pdbx_struct_oper_list.id 
 * _pdbx_struct_oper_list.type 
 * _pdbx_struct_oper_list.symmetry_operation
 * _pdbx_struct_oper_list.matrix[1][1] 
 * _pdbx_struct_oper_list.matrix[1][2] 
 * _pdbx_struct_oper_list.matrix[1][3] 
 * _pdbx_struct_oper_list.vector[1] 
 * _pdbx_struct_oper_list.matrix[2][1] 
 * _pdbx_struct_oper_list.matrix[2][2] 
 * _pdbx_struct_oper_list.matrix[2][3] 
 * _pdbx_struct_oper_list.vector[2] 
 * _pdbx_struct_oper_list.matrix[3][1] 
 * _pdbx_struct_oper_list.matrix[3][2] 
 * _pdbx_struct_oper_list.matrix[3][3] 
 * _pdbx_struct_oper_list.vector[3] 
 * _pdbx_struct_oper_list.name 
 * </pre>
 */
@XmlAccessorType(XmlAccessType.PROPERTY)
public class PdbxStructOperList implements Serializable{

	
	private static final long serialVersionUID = 8933552854747969787L;

	@Override
	public String toString() {
		return "PdbxStructOperList [id=" + id + ", type=" + type + "]";
	}


	private String id;

	private String type;
	
	private String symmetry_operation;
	
	@CIFLabel(label="matrix[1][1]")
	String matrix11;
	@CIFLabel(label="matrix[1][2]")
	String matrix12;
	@CIFLabel(label="matrix[1][3]")
	String matrix13;
	
	@CIFLabel(label="vector[1]")
	String vector1;

	@CIFLabel(label="matrix[2][1]")
	String matrix21;
	@CIFLabel(label="matrix[2][2]")
	String matrix22;
	@CIFLabel(label="matrix[2][3]")
	String matrix23;
	
	@CIFLabel(label="vector[2]")
	String vector2;
	
	@CIFLabel(label="matrix[3][1]")
	String matrix31;
	@CIFLabel(label="matrix[3][2]")
	String matrix32;
	@CIFLabel(label="matrix[3][3]")
	String matrix33;
	
	@CIFLabel(label="vector[3]")
	String vector3;

	String name;
	
	
	// from here fields that are not in the cif category


	public PdbxStructOperList(){

	}

	@XmlAttribute
	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	@XmlAttribute
	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setMatrix11(String val){
		this.matrix11 = val;
	}
	public void setMatrix21(String val){
		this.matrix21 = val;
	}
	public void setMatrix31(String val){
		this.matrix31 = val;
	}

	public void setMatrix12(String val){
		this.matrix12 = val;
	}
	public void setMatrix22(String val){
		this.matrix22 = val;
	}
	public void setMatrix32(String val){
		this.matrix32 = val;
	}
	public void setMatrix13(String val){
		this.matrix13 = val;
	}
	public void setMatrix23(String val){
		this.matrix23 = val;
	}
	public void setMatrix33(String val){
		this.matrix33 =val;
	}
	
	public void setName(String name) {
		this.name = name;
	}

	public String getVector1() {
		return vector1;
	}
	public void setVector1(String vector1) {
		this.vector1 = vector1;
	}
	public String getVector2() {
		return vector2;
	}
	public void setVector2(String vector2) {
		this.vector2 = vector2;
	}
	public String getVector3() {
		return vector3;
	}
	public void setVector3(String vector3) {
		this.vector3 = vector3;
	}
	public String getName() {
		return name;
	}
	public String getSymmetry_operation() {
		return symmetry_operation;
	}
	public void setSymmetry_operation(String symmetry_operation) {
		this.symmetry_operation = symmetry_operation;
	}
	@XmlElement
	public String getMatrix11(){
		return matrix11;
	}
	@XmlElement
	public String getMatrix21(){
		return matrix21;
	}
	@XmlElement
	public String getMatrix31(){
		return matrix31;
	}
	@XmlElement
	public String getMatrix12(){
		return matrix12;
	}
	@XmlElement
	public String getMatrix22(){
		return matrix22;
	}
	@XmlElement
	public String getMatrix32(){
		return matrix32;
	}
	@XmlElement
	public String getMatrix13(){
		return matrix13;
	}
	@XmlElement
	public String getMatrix23(){
		return matrix23;
	}
	@XmlElement
	public String getMatrix33(){
		return matrix33;
	}
}
