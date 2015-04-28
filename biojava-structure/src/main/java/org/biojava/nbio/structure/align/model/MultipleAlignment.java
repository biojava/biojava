package org.biojava.nbio.structure.align.model;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A MultipleAlignment is a Data Structure to store the core information of a multiple structure alignment, as a return type.
 * Each alignment is described as a collection of {@link BlockSet}s that define the aligned positions, 
 * a collection of structure identifiers (i,e. Atom arrays), information about the 3D superimposition in {@link Pose},
 * and creation properties (algorithm, version, etc).
 * A collection of MultipleAlignments that share the Atom arrays and creation properties form a {@link MultipleAlignmentEnsemble}.
 * Every MultipleAlignment has a {@link MultipleAlignmentEnsemble} as its parent.
 *
 * @author Aleix Lafita
 * 
 */
public interface MultipleAlignment extends Cloneable{
	
	/**
	 * Creates and returns an identical copy of this object.
	 * @return MultipleAlignment identical copy of this object.
	 */
	public Object clone();
	
	/** 
     * Returns the parent Ensemble of the MultipleAlignment.
     * Returns null if there is no referenced object.
     * @return MultipleAlignmentEnsemble the parent MultipleAlignment of the BlockSet, or null.
     * @see #setParent(MultipleAlignmentEnsemble)
     */
	public MultipleAlignmentEnsemble getParent();
	
	/** 
     * Set the back-reference to its parent Ensemble.
     * @param parent the parent MultipleAlignmentEnsemble.
     * @see #getParent()
     */
	public void setParent(MultipleAlignmentEnsemble parent);
	
	/**
	 * Returns the name of the multiple structure alignment algorithm that created this MSTA object.
	 * @return String name of the algorithm.
	 */
	public String getAlgorithmName();

	/**
	 * Returns the version of the algorithm used to generate this MSTA object.
	 * @return String version of the algorithm.
	 */
	public String getVersion();

	/**
	 * Returns the creation time of this MSTA object, in milliseconds.
	 * @return long creation time.
	 */
	public long getIoTime();

	/**
	 * Returns the running time of the structure alignment calculation, in milliseconds.
	 * @return long running time of the calculation.
	 * @see #getIoTime()
	 */
	public long getCalculationTime();

	/**
	 * Returns the structure alignment object ID.
	 * @return long structure alignment ID.
	 * @see #setId(long)
	 */
	public long getId();
	
	/**
	 * Returns a List containing the names of the structures aligned (i.e.: PDB code, SCOP domain, etc.).
	 * They are in the same order as in the Atom List and alignment List (same index number for same structure).
	 * @return List of String names of the structures
	 * @see #getAtomArrays()
	 */
	public List<String> getStructureNames();

	/**
	 * Returns the List of Atom arrays. Every structure has an Atom array associated.
	 * @return List of Atom[].
	 * @see #getStructureNames()
	 */
	public List<Atom[]> getAtomArrays();
	
	/**
	 * Returns the BlockSet List of the multiple structure alignment.
	 * Initializes the variable if it is null.
	 * @return List of BlockSets that describe the aligned residues of all the structures.
	 * @see #getBlockSetNum()
	 * @see #setBlockSets(List)
	 */
	public List<BlockSet> getBlockSets();

	/**
	 * Sets the List of BlockSet List of the specified alignment. The optimal alignment is always stored at position 0.
	 * @param blockSets the List of BlockSets that describe the aligned residues.
	 * @see #getBlockSets()
	 */
	public void setBlockSets(List<BlockSet> blockSets);

	/**
	 * Returns the List of Strings that represent the multiple sequence alignment of all the structures.
	 * @return List of Strings multiple sequence alignment
	 * @see #updateAlnSequences()
	 */
	public List<String> getAlnSequences();

	/**
	 * Calculates and sets the List of Strings that represent the multiple sequence alignment of all the structures.
	 * @see #getAlnSequences()
	 */
	public void updateAlnSequences();
	
	/**
	 * Returns the Pose object that stores the 3D superimposition information of the multiple structure alignment.
	 * @return Pose the 3D superimposition information of the alignment
	 * @throws StructureAlignmentException 
	 * @see #updatePose(PoseMethod)
	 */
	public Pose getPose() throws StructureAlignmentException;
	
	/**
	 * Calculates and sets the new 3D superimposition information in the Pose of the multiple alignment.
	 * Methods: REFERENCE (align everything to the first structure, the master), 
	 * 			MEDIAN (take the closest structure to all others in average as the master and align everything to it),
	 * 			CONSENSUS (build a consensus structure and align everything to it)
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superimposition.
	 * @throws StructureException
	 * @throws StructureAlignmentException
	 * @see #getPose()
	 */
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException;
	
	/**
	 * Returns the score of the multiple alignment used by the algorithm.
	 * @return double score of the multiple alignment used by the algorithm.
	 * @see #setAlgScore(double)
	 */
	public double getAlgScore();

	/**
	 * Sets the score of the multiple alignment used by the algorithm.
	 * @param algScore score of the multiple alignment used by the algorithm.
	 * @see #getAlgScore()
	 */
	public void setAlgScore(double algScore);

	/**
	 * Returns the probability used by the algorithm.
	 * @return double probability
	 * @see #setProbability(double)
	 */
	public double getProbability();

	/**
	 * Sets the probability used by the algorithm.
	 * @param probability
	 * @see #getProbability()
	 */
	public void setProbability(double probability);
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignment.
	 * @return int number of aligned structures
	 * @throws StructureAlignmentException 
	 * @see #length()
	 * @see #getCoreLength()
	 * @see #getBlockSetNum()
	 */
	public int size() throws StructureAlignmentException;
	
	/**
	 * Returns the total number of aligned residues (columns) in the multiple alignment: the sum of all BlockSet lengths.
	 * @return int the total number of aligned residues in the alignment.
	 * @throws StructureAlignmentException if there are no BlockSets.
	 * @see #updateLength()
	 * @see #getCoreLength()
	 * @see #size()
	 * @see #getBlockSetNum()
	 */
	public int length() throws StructureAlignmentException;
	
	/**
	 * Calculates and sets the total number of aligned residues (columns) in the alignment: the sum of all BlockSet lengths.
	 * @throws StructureAlignmentException
	 * @see #length()
	 * @see #updateCoreLength()
	 * @see #getCoreLength()
	 */
	public void updateLength() throws StructureAlignmentException;
	
	/**
	 * Returns the number of aligned residues (columns) without gaps in the alignment: the sum of all BlockSet core lengths.
	 * @return int the total number of aligned residues.
	 * @throws StructureAlignmentException if there are no BlockSets.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int getCoreLength() throws StructureAlignmentException;
	
	/**
	 * Calculates and sets the number of aligned residues (columns) without gaps in the alignment: the sum of all BlockSet core lengths.
	 * @throws StructureAlignmentException
	 * @see #getCoreLength()
	 * @see #length()
	 * @see #updateLength()
	 */
	public void updateCoreLength() throws StructureAlignmentException;
	
	/**
	 * Returns the number of BlockSets in the multiple structure alignment.
	 * @return int number of BlockSets.
	 * @throws StructureAlignmentException if the MultipleAlignment is empty.
	 * @see #length()
	 * @see #size()
	 */
	public int getBlockSetNum() throws StructureAlignmentException;

	/**
	 * Returns the List containing the interatomic distance Matrix of each structure.
	 * @return List of Matrix interatomic distance matrices.
	 * @throws StructureAlignmentException 
	 */
	public List<Matrix> getDistanceMatrix() throws StructureAlignmentException;
	
	/**
	 * Calls all the update methods for the Cache variables.
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superimposition.
	 * @throws StructureAlignmentException
	 * @throws StructureException 
	 * @see #updateLength()
	 * @see #updateCoreLength()
	 * @see #updatePose(PoseMethod)
	 * @see #updateAlnSequences()
	 */
	public void updateCache(PoseMethod method) throws StructureAlignmentException, StructureException;

}
