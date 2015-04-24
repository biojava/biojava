package org.biojava.nbio.structure.align.model;

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
public interface Pose extends Cloneable{

	/**
	 * Creates and returns an identical copy of this object.
	 * @return Pose identical copy of this object.
	 */
	public Object clone();
	
	/** 
     * Set the back-reference to its parent BlockSet.
     * @param parent the parent BlockSet.
     * @see #getBlockSet()
     */
	public void setBlockSet(BlockSet parent);

	/** 
     * Returns the parent BlockSet of the Pose.
     * Returns null if there is no referenced object. 
     * @return BlockSet the parent BlockSet of the Pose, or null.
     * @see #setBlockSet(BlockSet)
     */
	public BlockSet getBlockSet();

	/**
	 * Returns the List of rotation Matrices. One Matrix for each structure.
	 * @return List of rotation Matrix
	 * @see #setRotation(List)
	 */
	public List<Matrix> getRotation();

	/**
	 * Set the List of rotation Matrices. One Matrix for each structure.
	 * @param rotation the List of rotation Matrix
	 * @see #getRotation()
	 */
	public void setRotation(List<Matrix> rotation);

	/**
	 * Return the List of translation Vectors as Atom coordinates. One vector for every structure.
	 * @return List of translation Atoms
	 * @see #setTranslation(List)
	 */
	public List<Atom> getTranslation();

	/**
	 * Set the List of translation Vectors as Atom coordinates. One vector for every structure.
	 * @param translation the List of translation Atoms
	 * @see #getTranslation()
	 */
	public void setTranslation(List<Atom> translation);

	/**
	 * The background distances are the distances between every Atom of two aligned structures (Dot-Plot).
	 * Return the background distance Matrices as a double List. One Matrix for every structure pairwise combination.
	 * @return List double with the background distance Matrices
	 * @see #setBackDistMatrix(List)
	 */
	public List<List<Matrix>> getBackDistMatrix();

	/**
	 * The background distances are the distances between every Atom of two aligned structures (Dot-Plot).
	 * Set the background distance Matrices as a double List. One Matrix for every structure pairwise combination.
	 * @param backDistMatrix List double with the background distance Matrices
	 * @see #getBackDistMatrix(List)
	 */
	public void setBackDistMatrix(List<List<Matrix>> backDistMatrix);
	
	/**
	 * Returns the number of aligned structures in the Pose.
	 * @return int number of aligned structures
	 */
	public int size();
	
}
