package org.biojava.bio.structure.io.mmcif.model;


import java.io.Serializable;
import java.util.Arrays;


import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;

import org.biojava.bio.structure.jama.Matrix;

@XmlAccessorType(XmlAccessType.PROPERTY)
public class PdbxStructOperList implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 8933552854747969787L;

	@Override
	public String toString() {
		return "PdbxStructOperList [id=" + id + ", type=" + type + ", matrix="
				+ matrix + ", vector=" + Arrays.toString(vector) + "]";
	}

	
	private String id;
	
	private String type;
	
	private Matrix matrix;
	
	
	private double[] vector;

	public PdbxStructOperList(){
		matrix =  Matrix.identity(3,3);
		
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
	
	public void setMatrix11(double val){
		matrix.set(0,0,val);
	}
	public void setMatrix21(double val){
		matrix.set(1,0,val);
	}
	public void setMatrix31(double val){
		matrix.set(2,0,val);
	}
	
	public void setMatrix12(double val){
		matrix.set(0,1,val);
	}
	public void setMatrix22(double val){
		matrix.set(1,1,val);
	}
	public void setMatrix32(double val){
		matrix.set(2,1,val);
	}
	public void setMatrix13(double val){
		matrix.set(0,2,val);
	}
	public void setMatrix23(double val){
		matrix.set(1,2,val);
	}
	public void setMatrix33(double val){
		matrix.set(2,2,val);
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
