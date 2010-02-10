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
 * A <code>RegistryException</code> thrown when the registry cannot
 * find an implementation of a requested <code>SequenceDB</code>.
 *
 * @author Brian Gilman
 * @author Thomas Down
 * @author Keith James
 * @author Matthew Pocock
 * @version $Revision$
 */
public class RegistryException extends Exception {

    /**
     * Creates a new <code>RegistryException</code> with no detail
     * message.
     */
    public RegistryException(){
	super();
    }

    /**
     * Creates a new <code>RegistryException</code> with the specified
     * detail message.
     *
     * @param message a <code>String</code>.
     */
    public RegistryException(String message){
	super(message);
    }

    /**
     * Creates a new <code>RegistryException</code> with no detail
     * message, wrapping another <code>Throwable</code>.
     *
     * @param t a <code>Throwable</code>.
     */
    public RegistryException(Throwable t) {
	super(t);
    }

    /**
     * Creates a new <code>RegistryException</code> with the specified
     * detail message, wrapping another <code>Throwable</code>.
     *
     * @param t a <code>Throwable</code>.
     * @param message a <code>String</code>.
     * @deprecated use new RegistryException(message, cause)
     */
    public RegistryException(Throwable t, String message) {
	this(message, t);
    }
    
    public RegistryException(String message, Throwable t) {
      super(message, t);
    }
}
