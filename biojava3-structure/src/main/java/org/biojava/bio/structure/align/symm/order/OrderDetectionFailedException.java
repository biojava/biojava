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

/**
 * Order-detection failed.
 * @author dmyersturnbull
 */
public class OrderDetectionFailedException extends Exception {
	private static final long serialVersionUID = -7040421412578699838L;

	public OrderDetectionFailedException() {
		super();
	}

	public OrderDetectionFailedException(String message, Throwable cause) {
		super(message, cause);
	}

	public OrderDetectionFailedException(String message) {
		super(message);
	}

	public OrderDetectionFailedException(Throwable cause) {
		super(cause);
	}

}
