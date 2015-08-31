package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * A method to decide the order of symmetry (number of subunits) 
 * given a structure self-alignment, calculated by CE-Symm.
 * 
 * @author dmyersturnbull
 * @since 4.2.0
 * 
 */
public interface OrderDetector {

	public int calculateOrder(AFPChain afpChain, Atom[] ca) 
			throws RefinerFailedException;
	
}
