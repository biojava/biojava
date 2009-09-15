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

package org.biojava.bio.program.das.dasalignment;

import org.biojava.bio.BioException;

/**
 * An exception of one of the DAS classes.
 *
 * @author Andreas Prlic, Thomas Down, Benjamin Schuster-Böckler
 */

public class DASException extends BioException {
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
     * Constructs a DASException object.
     *
     * @param s  a String ...
     */

    public DASException(String s) {
	super(s);
    }

    /**
     * Constructs a DASException object.
     *
     * @param t  a Throwable object
     * @param s  a String ...
     */
    public DASException (Throwable t, String s) {
	super(s, t);
    }

    /**
     * Constructs a DASException object.
     *
     * @param t  a Throwable object
     */
    public DASException (Throwable t) {
	super(t);
    }
}
