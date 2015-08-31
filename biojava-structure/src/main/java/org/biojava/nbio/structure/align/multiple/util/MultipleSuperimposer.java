package org.biojava.nbio.structure.align.multiple.util;

import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;

/**
 * Interface for Multiple Alignment superposition algorithms.
 * <p>
 * There can be several implementations of a superimposer, 
 * because there is not a unique superposition of multiple 
 * structures given their residue equivalencies, and their
 * running times are asymptotically different.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public interface MultipleSuperimposer {
	
	/**
	 * Superimpose all structures from a {@link MultipleAlignment}. The 
	 * superposition is done for all individual BlockSets. If there is 
	 * only one BlockSet.
	 * <p>
	 * At a minimum, this should set the transformation matrices for 
	 * the individual {@link BlockSet#setTransformations(List) BlockSet}s.
	 * <p>
	 * This method only calculates and sets the transformation 4D Matrices.
	 * If any score is needed it should be calculated and set separately 
	 * afterwards with {@link MultipleAlignmentScorer}.
	 * 
	 * @param alignment MultipleAlignment specifying the aligned residues (via
	 *  		the {@link MultipleAlignment#getBlockSets() blocksets}) and the 
	 *  		atoms to align (via the {@link MultipleAlignment#getEnsemble() 
	 *  		ensemble}).
	 * 
	 * @throws StructureException
	 */
	public void superimpose(MultipleAlignment alignment) 
			throws StructureException;
}
