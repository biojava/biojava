package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A MultipleAlignment is a Data Structure to store the core information of a multiple structure alignment, 
 * as a return type.
 * Each alignment is described as a collection of: <p>
 * - {@link BlockSet}s that define the aligned positions, <p>
 * - Structure identifiers (i,e. Atom arrays, structure names), <p>
 * - Information about the 3D superimposition in a 4D transformation matrix, <p>
 * - Creation properties (algorithm, version, etc). <p>
 * A collection of MultipleAlignments that share the same structures and creation properties 
 * are part of the same {@link MultipleAlignmentEnsemble}.
 * Every MultipleAlignment has a {@link MultipleAlignmentEnsemble} as its parent, which contains this information
 *
 * @author Aleix Lafita
 * 
 */
public interface MultipleAlignment extends ScoresCache {
	
	/**
	 * Creates and returns an identical copy of this alignment, including a deep
	 * clone of all constituent blocks sets.
	 * @return MultipleAlignment identical copy of this object.
	 */
	public MultipleAlignment clone();
	
	/** 
     * Returns the parent Ensemble of the MultipleAlignment.
     * Returns null if there is no referenced object.
     * @return MultipleAlignmentEnsemble the parent MultipleAlignment of the BlockSet, or null.
     * @see #setParent(MultipleAlignmentEnsemble)
     */
	public MultipleAlignmentEnsemble getEnsemble();
	
	/** 
     * Set the back-reference to its parent Ensemble.
     * <p>
     * Neither removes this alignment from its previous ensemble, if any, nor
     * adds it to the new parent. Calling code should assure that links to
     * and from the ensemble are consistent and free of memory leaks.
     * @param parent the parent MultipleAlignmentEnsemble.
     * @see #getEnsemble()
     */
	public void setEnsemble(MultipleAlignmentEnsemble parent);

	/**
	 * Returns the BlockSet List of the multiple structure alignment.
	 * Initializes the variable if it is null.
	 * @return List of BlockSets that describe the aligned residues of all the structures.
	 * @see #getBlockSetNum()
	 * @see #setBlockSets(List)
	 */
	public List<BlockSet> getBlockSets();

	/**
	 * Sets the List of BlockSet List of the specified alignment.
	 * @param blockSets the List of BlockSets that describe the aligned residues.
	 * @see #getBlockSets()
	 */
	public void setBlockSets(List<BlockSet> blockSets);

	/**
	 * Convenience method to get a list of all blocks from all blocksets
	 * @return List of alignment Blocks
	 */
	public List<Block> getBlocks();

	/**
	 * Returns a transformation 4D matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 * <p>
	 * Individual BlockSets may override the transformation matrix for particular
	 * substructures. Flexible alignments will generally return null from
	 * this method, while rigid-body methods would typically store the global
	 * matrices here and return null for {@link BlockSet#getTransformations()}.
	 * @return the 3D superimposition information of the alignment
	 */
	public List<Matrix4d> getTransformations();
	
	/**
	 * Set a new superposition for the structures.
	 * <p>
	 * This may trigger other properties to update which depend on the superposition.
	 * In particular, the list of scores should be reset by implementations after
	 * changing the transformation matrices.
	 * @param matrices 4D
	 * @throws IllegalArgumentException when the size of the alignment and the size of transformations do not match.
	 */
	public void setTransformations(List<Matrix4d> transformations);
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignment.
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getCoreLength()
	 * @see #getBlockSetNum()
	 */
	public int size();

	/**
	 * Returns the total number of aligned residues (columns) in the multiple alignment: 
	 * the sum of all BlockSet lengths.
	 * @return int the total number of aligned residues in the alignment.
	 * @see #updateLength()
	 * @see #getCoreLength()
	 * @see #size()
	 * @see #getBlockSetNum()
	 */
	public int length();

	/**
	 * Returns the number of aligned residues (columns) without gaps in the alignment: 
	 * the sum of all BlockSet core lengths.
	 * @return int the total number of aligned residues.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int getCoreLength();
	
	/**
	 * Clear scores and other properties which depend on the specific alignment.
	 * <p>
	 * This can free memory and ensures consistency for cached variables.
	 */
	public void clear();
}
