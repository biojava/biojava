package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.jama.Matrix;

/**
 * A Pose is a Data Structure that stores the coordinates of a 3D superimposition of the structures in a multiple alignment.
 * It is associated with a {@link BlockSet} instance, named as parent. It is used as a cash variable in the BlockSet.
 *
 * @author Aleix Lafita
 * 
 */
public class Pose implements Serializable, Cloneable{

	private static final long serialVersionUID = -8309408466325388360L;
	
	BlockSet parent;

	List<Matrix> transforms;				//transform Matrix for every structure to calculate the aligned 3D superimposition
													//SIZE: n Matrices of s*s size (s=structure length)
	List<List<Matrix>> distanceMatrix; 					//A (n*l)*(n*l) matrix that stores the distance between every pair of aligned residues
													//n=nr. structures; l=alignment length
	
	
	/**
	 * Constructor
	 */
	public Pose(BlockSet blockSet) {
		
		parent = blockSet;
		transforms = new ArrayList<Matrix>();
		distanceMatrix = new ArrayList<List<Matrix>>();
		
		//Initialize the distanceMatrix with matrices of zeros
		for (int i=0; i<parent.getMultipleAlignment().getSize(); i++){
			distanceMatrix.add(new ArrayList<Matrix>());
			for (int j=0; j<parent.getMultipleAlignment().getSize(); j++){
				int len = parent.getLength();
				distanceMatrix.get(i).add(new Matrix(len,len));
			}
		}
		
	}
	
	/**
	 * Copy constructor
	 */
	public Pose(Pose p) {
		
		this.parent = p.parent;
		transforms = new ArrayList<Matrix>(p.transforms);
		distanceMatrix = new ArrayList<List<Matrix>>(p.distanceMatrix);
		
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	public Pose clone(){
		
		return new Pose(this);
	}

	@Override
	public String toString() {
		return "Pose [parent=" + parent + ", transforms=" + transforms
				+ ", distanceMatrix=" + distanceMatrix + "]";
	}

	//Getters and Setters **************************************************************************************
	
	public BlockSet getBlockSet() {
		return parent;
	}

	public void setBlockSet(BlockSet parent) {
		this.parent = parent;
	}

	public List<Matrix> getTransforms() {
		return transforms;
	}

	public void setTransforms(List<Matrix> transforms) {
		this.transforms = transforms;
	}

	public List<List<Matrix>> getDistanceMatrix() {
		return distanceMatrix;
	}

	public void setDistanceMatrix(List<List<Matrix>> distanceMatrix) {
		this.distanceMatrix = distanceMatrix;
	}
}
