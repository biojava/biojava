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
