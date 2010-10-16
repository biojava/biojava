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

package org.biojava.utils.lsid;

/**
 * Exception thrown in the event of an error in
 * parsing a LSID-formatted string.
 *
 * @author Michael Heuer
 */
public class LifeScienceIdentifierParseException
    extends Exception {

    /**
     * Construct a new parse exception with no message.
     */
    public LifeScienceIdentifierParseException() {
        super();
    }

    /**
     * Construct a new parse exception with the
     * specified error message.
     *
     * @param message error message
     */
    public LifeScienceIdentifierParseException(String message) {
        super(message);
    }
}
