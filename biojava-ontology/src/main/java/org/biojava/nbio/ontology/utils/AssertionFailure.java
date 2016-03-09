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
package org.biojava.nbio.ontology.utils;

/**
 * An unchecked exception representing an Assertion failure.
 *
 * <p>Assertion failures should be raised when code finds itself in a state that
 * should be impossible. It should not be raised in response to any predictable
 * error condition. Assertion failures indicate that something has gone
 * badly wrong, and that the assumptions under which library code has been
 * developed are not holding.</p>
 *
 * <p>This extends {@link java.lang.AssertionError}, adding convenient
 * constructors with messages and causes.</p>
 *
 *
 * Your application may exit due to one of these being thrown. This usualy
 * indicates that something is badly wrong with library code. It should never
 * be raised in response to invalid arguments to methods, or incorrectly
 * formatted data. It is not your fault. Report the error to the mailing list,
 * or who ever else is responsible for the library code you are using.
 *
 *
 * Under some rare circumstances, you may wish to catch assertion failures. For
 * example, when debugging library code, or when the success or failure of an
 * opperation is utterly inconsequential. Ignoring assertion failures
 * out-of-hand is a sure-fire way to make your code buggy.
 *
 *
 * Raise AssertionFailure in your code when something that should be impossible
 * has happened. For example, if you have checked the alphabet of a symbol list
 * you are working with, and somewhere further down an IllegalSymbolException
 * is raised, then this is an assertion failure.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class AssertionFailure
extends AssertionError {
	public AssertionFailure(String message) {
		super(message);
	}

	public AssertionFailure(Throwable cause) {
		initCause(cause);
	}

	public AssertionFailure(String message, Throwable cause) {
		this(message);
		initCause(cause);
	}
}
