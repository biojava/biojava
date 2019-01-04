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

import javax.vecmath.Matrix4d;

/**
 * A BlockSet is a Data Structure to store a flexible alignment part of a
 * multiple alignment. It is described by a collection of {@link Block}s and
 * a transformation matrix for every structure.
 * <p>
 * It allows the description of non-topological and circularly permutated
 * flexible parts, thus being as general as possible, thanks to the multiple
 * {@link Block} format.
 * <p>
 * Every BlockSet has a unique transformation 4D matrix for every structure,
 * which describes the 3D superimposition of the structures in this particular
 * part of the alignment.
 * <p>
 * A collection of BlockSets, in a {@link MultipleAlignment}, allows the
 * description of alignments with several flexible parts.
 * Every BlockSet object is part of a {@link MultipleAlignment} instance,
 * its parent.
 *
 * @author Aleix Lafita
 * @author Spencer Bliven
 * @since 4.1.0
 *
 */
public interface BlockSet extends ScoresCache {

	/**
	 * Creates and returns an identical copy of this blockset,
	 * including a deep copy of all constituent {@link Block}s.
	 *
	 * @return BlockSet identical copy of this object.
	 */
	public BlockSet clone();

	/**
	 * Returns the parent MultipleAlignment of the BlockSet.
	 * Returns null if there is no referenced object.
	 *
	 * @return MultipleAlignment the parent MultipleAlignment of the BlockSet,
	 * or null.
	 * @see #setMultipleAlignment(MultipleAlignment)
	 */
	public MultipleAlignment getMultipleAlignment();

	/**
	 * Set the back-reference to its parent MultipleAlignment.
	 * <p>
	 * Neither removes this BlockSet from its previous alignment, if any, nor
	 * adds it to the new parent. Calling code should assure that links to
	 * and from the ensemble are consistent and free of memory leaks.
	 *
	 * @param parent the parent MultipleAlignment.
	 * @see #getMultipleAlignment()
	 */
	public void setMultipleAlignment(MultipleAlignment parent);

	/**
	 * Returns the List of alignment Blocks of the BlockSet.
	 * It initializes a new List of Blocks if it is null.
	 * @return List of alignment Blocks.
	 *
	 * @see #setBlocks(List)
	 */
	public List<Block> getBlocks();

	/**
	 * Set the List of alignment Blocks of the BlockSet.
	 * <p>
	 * Also calls {@link Block#setBlockSet(BlockSet)} for each argument
	 *
	 * @param blocks List of alignment Blocks.
	 * @see #getBlocks()
	 */
	public void setBlocks(List<Block> blocks);

	/**
	 * Returns a transformation matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 *
	 * @return the 3D superimposition information of the alignment
	 */
	public List<Matrix4d> getTransformations();

	/**
	 * Set a new superposition for the structures.
	 * This may trigger other properties to update which depend on the
	 * superposition.
	 *
	 * @param matrices
	 */
	public void setTransformations(List<Matrix4d> transformations);

	/**
	 * Returns the total number of aligned residues (columns) in the alignment:
	 * the sum of all Block lengths.
	 *
	 * @return int the total number of aligned residues.
	 * @see #getCoreLength()
	 * @see #size()
	 */
	public int length();

	/**
	 * Returns the number of aligned residues (columns) without gaps in the
	 * alignment: the sum of all Block core lengths.
	 *
	 * @return int the total number of aligned residues.
	 * @see #length()
	 * @see #size()
	 */
	public int getCoreLength();

	/**
	 * Returns the number of non null positions (residues) of each structure in
	 * the alignment Block Set. The values can be used to compute the coverages.
	 *
	 * @return List of residue counts for each structure
	 */
	public List<Integer> getAlignResCounts();

	/**
	 * Returns the number of aligned structures in the BlockSet.
	 *
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getCoreLength()
	 */
	public int size();

	/**
	 * Clear scores and other properties which depend on the specific
	 * alignment. This frees memory and ensures consistency of the cached
	 * variables.<p>
	 * Recursively clears the memeber Blocks.
	 */
	public void clear();
}
