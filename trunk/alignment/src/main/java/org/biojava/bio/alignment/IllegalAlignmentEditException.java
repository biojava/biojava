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

package org.biojava.bio.alignment;

import org.biojava.bio.BioRuntimeException;

/**
 * 
 * <p>
 * 
 * The usual reason for throwing an IllegalAlignmentEditException is that you
 * are
 * 
 * trying to shift a group of bases in such a way that it would require deleting
 * bases.
 * 
 * Only gaps can be deleted.
 * 
 * </p>
 * 
 * 
 * 
 * @author David Waring
 */

public class IllegalAlignmentEditException extends BioRuntimeException {

	/**
	 * Generated Serial Version ID.
	 */
	private static final long serialVersionUID = 4662428071284597336L;

	/**
	 * 
	 * Just make the exception.
	 */

	public IllegalAlignmentEditException() {
		super();
	}

	/**
	 * 
	 * Make the exception with a message.
	 */

	public IllegalAlignmentEditException(String message) {
		super(message);
	}

	public IllegalAlignmentEditException(Throwable t) {
		super(t);
	}

	public IllegalAlignmentEditException(Throwable t, String message) {
		super(message, t);
	}

}
