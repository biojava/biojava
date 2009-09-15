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

package org.biojava.directory;

/**
 * <code>ProviderNotFoundException</code> is thrown when a sequence
 * database provider can not be located.
 *
 * @author Keith James
 * @author Brian Gilman
 */
public class ProviderNotFoundException extends RegistryException {

    /**
     * Creates a <code>ProviderNotFoundException</code> with no detail
     * message.
     */
    public ProviderNotFoundException() {
	super();
    }

    /**
     * Creates a <code>ProviderNotFoundException</code> with the
     * specified detail message.
     *
     * @param message a <code>String</code>.
     */
    public ProviderNotFoundException(String message) {
	super(message);
    }
}
    



