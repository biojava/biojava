package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * A method to decide the order of symmetry given a self-alignment.
 * 
 * @author dmyersturnbull
 *
 */
public interface OrderDetector {

	public int calculateOrder(AFPChain afpChain, Atom[] ca) 
			throws OrderDetectionFailedException;
	
}
