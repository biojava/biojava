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
 * Thrown to indicate an error in an Ontology object
 *
 * @author Matthew Pocock
 * @since 1.4
 */

public class OntologyException
extends Exception {

	private static final long serialVersionUID = 1L;

public OntologyException() {
		super();
	}

	public OntologyException(String message) {
		super(message);
	}

	public OntologyException(Throwable cause) {
		super(cause);
	}

	public OntologyException(String message, Throwable cause) {
		super(message, cause);
	}
}
