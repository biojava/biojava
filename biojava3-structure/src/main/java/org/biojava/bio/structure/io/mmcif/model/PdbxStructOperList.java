package org.biojava.bio.structure.io.mmcif.model;



import java.util.Arrays;

import org.biojava.bio.structure.jama.Matrix;

public class PdbxStructOperList {

	@Override
	public String toString() {
		return "PdbxStructOperList [id=" + id + ", type=" + type + ", matrix="
				+ matrix + ", vector=" + Arrays.toString(vector) + "]";
	}

	private String id;
	private String type;
	private Matrix matrix;

	private double[] vector;




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

	public double[] getVector() {
		return vector;
	}

	public void setVector(double[] vector) {
		this.vector = vector;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}






}
