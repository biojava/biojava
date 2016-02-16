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

/**
 * Refinement of the self-alignment failed.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 * 
 */
public class RefinerFailedException extends Exception {
	
	private static final long serialVersionUID = -3592155787060329421L;

	public RefinerFailedException() {
		super();
	}

	public RefinerFailedException(String message, Throwable cause) {
		super(message, cause);
	}

	public RefinerFailedException(String message) {
		super(message);
	}

	public RefinerFailedException(Throwable cause) {
		super(cause);
	}

}
