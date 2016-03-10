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

package org.biojava.nbio.ontology;


/**
 * Thrown to indicate that an ontology term is not acceptable or
 * appropriate in a given context
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.4
 */

public class InvalidTermException
		extends OntologyException
{

	private static final long serialVersionUID = 1L;

public InvalidTermException() {
		super();
	}

	public InvalidTermException(String message) {
		super(message);
	}

	public InvalidTermException(Throwable cause) {
		super(cause);
	}

	public InvalidTermException(String message, Throwable cause) {
		super(message, cause);
	}
}
