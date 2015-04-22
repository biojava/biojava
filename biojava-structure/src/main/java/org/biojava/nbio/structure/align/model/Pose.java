package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
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

	List<Matrix> rotationMatrix;			//rotation Matrix for every structure to calculate the aligned 3D superimposition
													//SIZE: n Matrices of s*s size (s=structure length) - one of the structures is the reference (transform = 0)
	List<Atom> translation;					//shift matrix for the atoms of every protein. There are n entries in the list.
	List<Matrix> distanceTables; 			//A list of n (l*l) matrices that store the distance between every pair of residues for every protein
													//n=nr. structures; l=alignment length 
													//(This variable does not change, because proteins stay the same. Should it be moved to the MultipleAlignment class?) TODO
	
	
	/**
	 * Constructor
	 */
	public Pose(BlockSet blockSet) {
		
		parent = blockSet;
		rotationMatrix = new ArrayList<Matrix>();
		translation = new ArrayList<Atom>();
		distanceTables = new ArrayList<Matrix>();
		
		//Initialize the distanceMatrix with matrices of zeros
		/*for (int j=0; j<parent.getMultipleAlignment().getSize(); j++){
			int len = parent.getMultipleAlignment().getAtomArrays().get(j).length;
			distanceTables.add(new Matrix(len,len));
		}*/
		
	}
	
	/**
	 * Copy constructor
	 */
	public Pose(Pose p) {
		
		this.parent = p.parent;
		rotationMatrix = new ArrayList<Matrix>(p.rotationMatrix);
		translation = new ArrayList<Atom>(p.translation);
		distanceTables = new ArrayList<Matrix>(p.distanceTables);
		
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	public Pose clone(){
		
		return new Pose(this);
	}

	@Override
	public String toString() {
		return "Pose [parent=" + parent + ", transforms=" + rotationMatrix
				+ ", distanceMatrix=" + distanceTables + "]";
	}

	//Getters and Setters **************************************************************************************
	
	public BlockSet getBlockSet() {
		return parent;
	}

	public void setBlockSet(BlockSet parent) {
		this.parent = parent;
	}

	public List<Matrix> getRotationMatrix() {
		return rotationMatrix;
	}

	public void setRotationMatrix(List<Matrix> transforms) {
		this.rotationMatrix = transforms;
	}

	public List<Matrix> getDistanceMatrix() {
		return distanceTables;
	}

	public void setDistanceMatrix(List<Matrix> distanceMatrix) {
		this.distanceTables = distanceMatrix;
	}

	public List<Atom> getTranslation() {
		return translation;
	}

	public void setTranslation(List<Atom> translation) {
		this.translation = translation;
	}
}
