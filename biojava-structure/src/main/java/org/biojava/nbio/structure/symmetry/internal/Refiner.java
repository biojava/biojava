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
