package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A MultipleAlignment is a Data Structure to store the core information of a 
 * multiple structure alignment, as a return type.
 * <p>
 * Each alignment is described as a collection of:
 * <ul><li>{@link BlockSet}s that define the aligned positions,
 * <li>Structure identifiers (i,e. Atom arrays, structure names),
 * <li>Information about the 3D superimposition in a set of 4D transformation 
 * matrices,
 * <li>Creation properties (algorithm, version, etc).
 * </ul>
 * A collection of MultipleAlignments that share the same structures and 
 * creation properties are part of the same {@link MultipleAlignmentEnsemble}.
 * Every MultipleAlignment has a {@link MultipleAlignmentEnsemble} as its 
 * parent, which contains the shared meta-information.
 *
 * @author Aleix Lafita
 * @author Spencer Bliven
 * @since 4.1.0
 * 
 */
public interface MultipleAlignment extends ScoresCache {
	
	/**
	 * Creates and returns an identical copy of this alignment, including a 
	 * deep copy of all constituent BlockSets.
	 * 
	 * @return MultipleAlignment identical copy of this object.
	 */
	public MultipleAlignment clone();
	
	/** 
     * Returns the parent Ensemble of the MultipleAlignment.
     * Returns null if there is no referenced object.
     * 
     * @return MultipleAlignmentEnsemble the parent MultipleAlignment of the 
     * BlockSet, or null.
     * @see #setEnsemble(MultipleAlignmentEnsemble)
     */
	public MultipleAlignmentEnsemble getEnsemble();
	
	/** 
     * Set the back-reference to its parent Ensemble.
     * <p>
     * Neither removes this alignment from its previous ensemble, if any, nor
     * adds it to the new parent. Calling code should assure that links to
     * and from the ensemble are consistent and free of memory leaks.
     * 
     * @param parent the parent MultipleAlignmentEnsemble.
     * @see #getEnsemble()
     */
	public void setEnsemble(MultipleAlignmentEnsemble parent);

	/**
	 * Returns the BlockSet List of the multiple structure alignment.
	 * Initializes the variable if it is null.
	 * 
	 * @return List of BlockSets that describe the aligned residues of all the
	 * structures.
	 * @see #getBlocks()
	 * @see #setBlockSets(List)
	 */
	public List<BlockSet> getBlockSets();

	/**
	 * Sets the List of BlockSet List of the specified alignment.
	 * 
	 * @param blockSets the List of BlockSets that describe the aligned 
	 * residues.
	 * @see #getBlockSets()
	 */
	public void setBlockSets(List<BlockSet> blockSets);

	/**
	 * Convenience method to get a List of all Blocks from all BlockSets.
	 * 
	 * @return List of Blocks
	 * @see #getBlockSets()
	 */
	public List<Block> getBlocks();

	/**
	 * Returns a transformation 4D matrix for each structure giving the
	 * 3D superposition information of the multiple structure alignment.
	 * <p>
	 * Individual BlockSets may override the transformation matrix for 
	 * particular parts of the alignment. Flexible alignments will generally
	 * return null from this method, while rigid-body methods would typically
	 * store the transformation matrices here as well as in the only BlockSet: 
	 * {@link BlockSet#getTransformations()}.
	 * 
	 * @return the 3D superposition information of the alignment or null
	 * 			if flexible
	 */
	public List<Matrix4d> getTransformations();
	
	/**
	 * Set a new superposition for the structures.
	 * <p>
	 * This may trigger other properties to update which depend on the 
	 * superposition. In particular, the list of scores should be reset by
	 * implementations after changing the transformation matrices.
	 * 
	 * @param matrices 4D
	 * @throws IllegalArgumentException when the size of the alignment and 
	 * the size of transformations do not match.
	 */
	public void setTransformations(List<Matrix4d> transformations);
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignment.
	 * 
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getCoreLength()
	 */
	public int size();

	/**
	 * Returns the total number of aligned residues (columns) in the multiple 
	 * alignment: the sum of all BlockSet lengths.
	 * 
	 * @return int the total number of aligned residues in the alignment.
	 * @see #getCoreLength()
	 * @see #size()
	 */
	public int length();

	/**
	 * Returns the number of aligned residues (columns) without gaps in the 
	 * alignment: the sum of all BlockSet core lengths.
	 * 
	 * @return int the total number of aligned residues.
	 * @see #length()
	 * @see #size()
	 */
	public int getCoreLength();
	
	/**
	 * Clear scores and other properties which depend on the specific 
	 * alignment. This frees memory and ensures consistency of the cached 
	 * variables.<p>
	 * Recursively clears member BlockSets.
	 */
	public void clear();
	
	/**
	 * Return a summary of the MultipleAlignment, containing the structures, 
	 * the lengths and the cached scores. Can be used as a header for the 
	 * differnt display options.
	 * 
	 * @return String header summary of the MultipleAlignment
	 */
	@Override
	public String toString();
}
