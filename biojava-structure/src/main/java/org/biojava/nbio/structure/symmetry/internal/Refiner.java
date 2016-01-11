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

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * A method to refine the AFP alignment from one or more alternative 
 * self-alignments in order to define consistent subunits.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 * 
 */
public interface Refiner {

	/**
	 * Returns a refined symmetry alignment, where the subunit 
	 * residues are aligned consistently and separated into the 
	 * blocks of the AFPChain.
	 * 
	 * @param afpAlignments List of returned self-alignments in CeSymm
	 * @param atoms coordinates of the structure
	 * @param order order of symmetry to use.
	 * @return AFPChain refined symmetry alignment
	 * @throws RefinerFailedException
	 * @throws StructureException
	 */
	public AFPChain refine(List<AFPChain> afpAlignments, Atom[] atoms,
			int order) throws RefinerFailedException, StructureException;
	
}
