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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

/**
 * A failure to load a testing resource in {@link ResourceList}.
 * @author dmyerstu
 */
public class ResourceException extends RuntimeException {

	public ResourceException() {
		super();
	}
	
	public ResourceException(String message) {
		super(message);
	}
	
	public ResourceException(String message, Exception cause) {
		super(message, cause);
	}

	public ResourceException(Exception cause) {
		super(cause);
	}
	
}
