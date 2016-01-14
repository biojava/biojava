/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.align.multiple;

import java.util.List;

/**
 * A Block is a Data Structure that stores aligned positions of a 
 * {@link MultipleAlignment} with the condition that residues are in a 
 * sequential order.
 * <p>
 * A collection of Blocks, named {@link BlockSet}, allows the description of 
 * circular permutations (CP) and non-topological alignments.
 * <p>
 * Every Block object is part of a {@link BlockSet} instance, its parent, 
 * which has in turn a {@link MultipleAlignment} instance as parent.
 * 
 * @author Aleix Lafita
 * @author Spencer Bliven
 * @since 4.1.0
 * 
 */
public interface Block extends ScoresCache {

	/**
	 * Creates and returns an identical copy of this block.
	 * 
	 * @return Block identical copy of this object.
	 */
	public Block clone();

	/** 
	 * Set the back-reference to its parent BlockSet.
	 * 
	 * @param parent the parent BlockSet.
	 * @see #getBlockSet()
	 */
	public void setBlockSet(BlockSet parent);

	/** 
	 * Returns the parent BlockSet of the Block.
	 * Returns null if there is no referenced object. 
	 * 
	 * @return BlockSet the parent BlockSet of the Block, or null.
	 * @see #setBlockSet(BlockSet)
	 */
	public BlockSet getBlockSet();

	/**
	 * Returns the double List containing the aligned residues for each 
	 * structure.<p>
	 * alignRes.get(structure).get(residue) = alignRes.get(size).get(length).
	 * 
	 * @return List a double List of aligned residues for each structure.
	 * @see #setAlignRes()
	 */
	public List<List<Integer>> getAlignRes();

	/**
	 * Set the double List containing the aligned residues for each structure.
	 * 
	 * @param alignRes a double List of Integers with the aligned residues.
	 * @see #getAlignRes()
	 */
	public void setAlignRes(List<List<Integer>> alignRes);

	/**
	 * Returns the total number of aligned positions (columns) in the Block.
	 * 
	 * @return int number of aligned residues.
	 * @see #getCoreLength();
	 * @see #size()
	 */
	public int length();

	/**
	 * Returns the number of aligned structures (rows) in the Block.
	 * 
	 * @return int number of aligned structures.
	 * @see #length()
	 * @see #getCoreLength()
	 */
	public int size();

	/**
	 * Returns the number of aligned positions (columns) without gaps in the 
	 * Block.
	 * 
	 * @return int number of aligned residues.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 */
	public int getCoreLength();

	/**
	 * Clear scores and other properties which depend on the specific 
	 * alignment. This frees memory and ensures consistency of the cached 
	 * variables.
	 */
	public void clear();

}
