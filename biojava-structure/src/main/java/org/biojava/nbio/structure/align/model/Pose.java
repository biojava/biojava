package org.biojava.nbio.structure.align.model;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A Pose is a Data Structure that stores information of the 3D superimposition of the structures in a multiple alignment.
 * It can be associated to a {@link BlockSet} instance, its parent, to describe the superimposition of one of the
 * flexible parts of a flexible alignment. It can also be associated to a {@link MultipleAlignment} instance, its parent,
 * to describe the global (final) superimposition of the structures.
 * The idea of a Pose is to be used as a cache variable for the 3D information and contain all the methods to calculate
 * (update) the 3D superimposition given the aligned residues. It only contains Getter methods, the variables are
 * automatically calculated when the Pose is updated (so the do not need to be set).
 *
 * @author Aleix Lafita
 * 
 */
public interface Pose extends Cloneable{

	/**
	 * Method used to calculate the 3D superposition.
	 */
	public static enum PoseMethod {
		REFERENCE,	//align everything to the first structure, the master. N superpositions.
		MEDIAN,		//take the closest structure to all others in average as the master and align everything to it. N^2 superpositions.
		CONSENSUS;	//build a consensus structure and align everything to it. N^2 superpositions, slower than MEDIAN.
		public static PoseMethod DEFAULT = REFERENCE;
	}
	
	/**
	 * Creates and returns an identical copy of this object.
	 * @return Pose identical copy of this object.
	 */
	public Object clone();

	/** 
     * Returns the parent Object of the Pose.
     * @return parent Object of the Pose.
     */
	public Object getParent();

	/**
	 * Returns the List of rotation Matrices. One Matrix for each structure.
	 * TODO is it better to update the Pose or to throw an exception if it is null?)
	 * @return List of rotation Matrix.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 */
	public List<Matrix> getRotation() throws StructureAlignmentException;
	
	/**
	 * Return the List of translation Vectors as Atom coordinates. One vector for every structure.
	 * @return List of translation Atoms.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 */
	public List<Atom> getTranslation() throws StructureAlignmentException;

	/**
	 * The background distances are the distances between every Atom of two aligned structures (Dot-Plot).
	 * Return the background distance Matrices as a double List. One Matrix for every structure pairwise combination.
	 * @return double List with the background distance Matrices.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 */
	public List<List<Matrix>> getBackDistMatrix() throws StructureAlignmentException;
	
	/**
	 * Returns the global RMSD of the 3D superimposition.
	 * @return double RMSD of the 3D superimposition.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 * @see #getRMSD(int)
	 */
	public double getRMSD() throws StructureAlignmentException;
	
	/**
	 * Returns the RMSD of the specified structure against all others.
	 * @return double RMSD of the specified structure against the others.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 * @see #getRMSD()
	 */
	public double getRMSD(int structureNr) throws StructureAlignmentException;
	
	/**
	 * Returns the global TM-score of the 3D superimposition.
	 * @return double TM-score of the 3D superimposition.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 * @see #getTMscore(int)
	 */
	public double getTMscore() throws StructureAlignmentException;
	
	/**
	 * Returns the TM-score of the specified structure against all others.
	 * @return double TM-score of the specified structure against the others.
	 * @throws StructureAlignmentException if Pose is empty. Call {@link #updatePose(PoseMethod)} first.
	 * @see #getTMscore()
	 */
	public double getTMscore(int structureNr) throws StructureAlignmentException;
	
	/**
	 * Calculates and sets all the Pose variables from the parent information.
	 * Methods: REFERENCE (align everything to the first structure, the master - T(n)=n), 
	 * 			MEDIAN (take the closest structure to all others in average as the master and align everything to it - T(n)=n^2),
	 * 			CONSENSUS (build a consensus structure and align everything to it - T(n)=n^2+C)
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superimposition.
	 * @throws StructureException
	 * @throws StructureAlignmentException
	 */
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException;
	
	/**
	 * Returns the number of aligned structures in the Pose.
	 * @return int number of aligned structures.
	 * @throws StructureAlignmentException
	 */
	public int size() throws StructureAlignmentException;

}
