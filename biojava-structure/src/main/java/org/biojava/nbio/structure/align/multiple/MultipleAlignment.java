package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import org.biojava.nbio.structure.Atom;

/**
 * A MultipleAlignment is a Data Structure to store the core information of a 
 * multiple structure alignment, as a return type.
 * <p>
 * Each alignment is described as a collection of:
 * <ul><li>{@link BlockSet}s that define the aligned positions and 3D 
 * 		superposition,
 * <li>Structure identifiers (i,e. Atom arrays, structure names),
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
	 * Returns the BlockSet with the specified index of the MultipleAlignment.
	 * Throws an Exception if the index is out of bounds, like accessing a 
	 * normal List.
	 * 
	 * @param index of the BlockSet
	 * @return BlockSets at the specified index
	 * @see #getBlocks()
	 * @see #getBlockSets()
	 */
	public BlockSet getBlockSet(int index);

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
	 * Modifications of this List will not alter the MultipleAlignment,
	 * but modifications to the Blocks will.
	 * 
	 * @return List of Blocks
	 * @see #getBlockSets()
	 */
	public List<Block> getBlocks();

	/**
	 * Returns the Block with the specified index of the MultipleAlignment.
	 * Throws an Exception if the index is out of bounds, like accessing a 
	 * normal List.
	 * 
	 * @param index of the BlockSet
	 * @return Block at the specified index
	 * @see #getBlocks()
	 * @see #getBlockSets()
	 */
	public Block getBlock(int index);

	/**
	 * Returns the array of Atoms for each structure from its parent
	 * Ensemble.
	 * Throws an Exception if the parent ensemble is null or the Atom 
	 * variables are not previously set.
	 * 
	 * @return List of Atom arrays
	 * @see #getEnsemble()
	 */
	public List<Atom[]> getAtomArrays();

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
	 * variables.
	 * <p>
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
