package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A Pose is a Data Structure that stores the coordinates of a 3D superimposition of the structures in a multiple alignment.
 * It is associated with a {@link BlockSet} instance, named as parent. It is used as a cache variable in the BlockSet.
 *
 * @author Aleix Lafita
 * 
 */
public class Pose implements Serializable, Cloneable{

	private static final long serialVersionUID = -8309408466325388360L;
	
	BlockSet parent;

	List<Matrix> rotation;			//rotation Matrix for every structure to calculate the 3D superimposition.
	List<Atom> translation;			//translation Vectors for the atoms of every structure.

	/**
	 * Constructor
	 */
	public Pose(BlockSet blockSet) {
		
		parent = blockSet;
		rotation = new ArrayList<Matrix>();
		translation = new ArrayList<Atom>();
		
	}
	
	/**
	 * Copy constructor
	 */
	public Pose(Pose p) {
		
		this.parent = p.parent;
		rotation = new ArrayList<Matrix>(p.rotation);
		translation = new ArrayList<Atom>(p.translation);
		
	}
	
	/**
	 * Creates and returns a copy of this object. Uses the copy constructor.
	 */
	@Override
	public Pose clone(){
		
		return new Pose(this);
	}

	@Override
	public String toString() {
		return "Pose [parent=" + parent + ", rotation=" + rotation
				+ ", translation=" + translation + "]";
	}

	//Getters and Setters **************************************************************************************
	
	/**
	 * Get the parent BlockSet that generated this Pose instance.
	 */
	public BlockSet getBlockSet() {
		return parent;
	}

	/**
	 * Set the link to the parent BlockSet of this Pose instance.
	 * @param parent
	 */
	public void setBlockSet(BlockSet parent) {
		this.parent = parent;
	}

	/**
	 * Get the list of rotation matrices. One Matrix for the atom rotation of each structure.
	 */
	public List<Matrix> getRotationMatrix() {
		return rotation;
	}

	/**
	 * Set the list of rotation matrices. One Matrix for the atom rotation of each structure.
	 * @param rotation
	 */
	public void setRotationMatrix(List<Matrix> rotationMatrix) {
		this.rotation = rotationMatrix;
	}

	/**
	 * Get the list of translation vectors as Atom coordinates. One vector for every structure.
	 */
	public List<Atom> getTranslation() {
		return translation;
	}

	/**
	 * Set the list of translation vectors as Atom coordinates. One vector for every structure.
	 * @param translation
	 */
	public void setTranslation(List<Atom> translation) {
		this.translation = translation;
	}
}
