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


import org.biojava.nbio.structure.jama.Matrix;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import java.io.Serializable;
import java.util.Arrays;

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
		return "PdbxStructOperList [id=" + id + ", type=" + type + ", matrix="
				+ matrix + ", vector=" + Arrays.toString(vector) + "]";
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
	
	@IgnoreField
	private Matrix matrix;

	@IgnoreField
	private double[] vector;

	public PdbxStructOperList(){
		matrix =  Matrix.identity(3,3);
		vector = new double[3];

	}
	@XmlAttribute
	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public Matrix getMatrix() {
		return matrix;
	}

	public void setMatrix(Matrix matrix) {
		this.matrix = matrix;
	}
	@XmlAttribute
	public double[] getVector() {
		return vector;
	}

	public void setVector(double[] vector) {
		this.vector = vector;
	}
	@XmlAttribute
	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public void setMatrix11(String val){
		matrix.set(0,0,Double.parseDouble(val));
	}
	public void setMatrix21(String val){
		matrix.set(1,0,Double.parseDouble(val));
	}
	public void setMatrix31(String val){
		matrix.set(2,0,Double.parseDouble(val));
	}

	public void setMatrix12(String val){
		matrix.set(0,1,Double.parseDouble(val));
	}
	public void setMatrix22(String val){
		matrix.set(1,1,Double.parseDouble(val));
	}
	public void setMatrix32(String val){
		matrix.set(2,1,Double.parseDouble(val));
	}
	public void setMatrix13(String val){
		matrix.set(0,2,Double.parseDouble(val));
	}
	public void setMatrix23(String val){
		matrix.set(1,2,Double.parseDouble(val));
	}
	public void setMatrix33(String val){
		matrix.set(2,2,Double.parseDouble(val));
	}
	
	public void setName(String name) {
		this.name = name;
	}

	public String getVector1() {
		return vector1;
	}
	public void setVector1(String vector1) {
		vector[0] = Double.parseDouble(vector1);
	}
	public String getVector2() {
		return vector2;
	}
	public void setVector2(String vector2) {
		vector[1] = Double.parseDouble(vector2);
	}
	public String getVector3() {
		return vector3;
	}
	public void setVector3(String vector3) {
		vector[2] = Double.parseDouble(vector3);
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
	public double getMatrix11(){
		return matrix.get(0,0);
	}
	@XmlElement
	public double getMatrix21(){
		return matrix.get(1,0);
	}
	@XmlElement
	public double getMatrix31(){
		return matrix.get(2,0);
	}
	@XmlElement
	public double getMatrix12(){
		return matrix.get(0,1);
	}
	@XmlElement
	public double getMatrix22(){
		return matrix.get(1,1);
	}
	@XmlElement
	public double getMatrix32(){
		return matrix.get(2,1);
	}
	@XmlElement
	public double getMatrix13(){
		return matrix.get(0,2);
	}
	@XmlElement
	public double getMatrix23(){
		return matrix.get(1,2);
	}
	@XmlElement
	public double getMatrix33(){
		return matrix.get(2,2);
	}
}
