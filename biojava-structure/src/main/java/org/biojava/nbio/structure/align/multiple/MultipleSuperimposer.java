package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import org.biojava.nbio.structure.StructureException;

public interface MultipleSuperimposer {
	/**
	 * Superimpose all structures from a {@link MultipleAlignment}.
	 * <p>
	 * At a minimum, this should set the transformation matrices for either
	 * the {@link MultipleAlignment#setTransformations(List) MultipleAlignment}
	 * (for rigid-body superpositions) or the individual {@link
	 * BlockSet#setTransformations(List) BlockSet}s.
	 * Implementations may also set various scores
	 * @param alignment MultipleAlignment specifying the aligned residues (via
	 *  the {@link MultipleAlignment#getBlockSets() blocksets}) and the atoms to
	 *  align (via the {@link MultipleAlignment#getEnsemble() ensemble}).
	 * @throws StructureException
	 * @throws StructureAlignmentException
	 */
	public void superimpose(MultipleAlignment alignment) throws StructureException, StructureAlignmentException;
}
