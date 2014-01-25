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
 * @since 3.0.8
 */
package org.biojava.bio.structure.align.symm.order;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * A method to decide the order of symmetry given a self-alignment.
 * @author dmyersturnbull
 *
 */
public interface OrderDetector {

	int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException;
	
}
