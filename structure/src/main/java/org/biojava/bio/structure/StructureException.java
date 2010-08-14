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

package org.biojava.bio.structure;


/**
 * An exception during the parsing of a PDB file.
 *
 * @author Andreas Prlic, Thomas Down, Benjamin Schuster-Boeckler
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */

public class StructureException extends Exception {
    
    public static final long serialVersionUID = 62946250941673093l;
    /**
     * Constructs a StructureException object.
     *
     * @param s  a String ...
     */
    public StructureException(String s) {
	super(s);
    }
    
    /**
     * Constructs a StructureException object.
     *
     * @param t  a Throwable object
     * @param message  a String ...
     */ 
    public StructureException(String message, Throwable t){
        super(message, t);
    }
    
    /**
     * Constructs a StructureException object.
     *
     * @param t  a Throwable object
     */
    public StructureException (Throwable t) {
	super(t);
    }
}
