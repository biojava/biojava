package org.biojava.nbio.structure.align.superimpose;

import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.BlockSet;
import org.biojava.nbio.structure.align.model.MultipleAlignment;
import org.biojava.nbio.structure.align.model.StructureAlignmentException;

public interface MultipleSuperimposer {
	/**
	 * Superimpose all structures from a multiple alignment.
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
