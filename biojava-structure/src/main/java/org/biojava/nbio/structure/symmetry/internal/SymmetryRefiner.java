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
package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;

/**
 * Interface for all symmetry refinement implementations.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public interface SymmetryRefiner {

	/**
	 * Returns a refined symmetry alignment, where the repeat residues are
	 * aligned consistently in a MultipleAlignment.
	 *
	 * @param selfAlignment
	 *            optimal self-alignment calculated by CeSymm
	 * @param atoms
	 *            coordinates of the structure
	 * @param order
	 *            order of symmetry to use
	 * @return MultipleAlignment refined symmetry alignment
	 * @throws RefinerFailedException
	 * @throws StructureException
	 */
	public MultipleAlignment refine(AFPChain selfAlignment, Atom[] atoms, int order)
			throws RefinerFailedException, StructureException;

}
