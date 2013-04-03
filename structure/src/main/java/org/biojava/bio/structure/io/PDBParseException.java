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

package org.biojava.bio.structure.io;

/**
 * An exception during the parsing of a PDB file.
 *
 * @author Andreas Prlic, Thomas Down, Benjamin Schuster-B&ouml;ckler
 */

public class PDBParseException extends Exception{
    public static final long serialVersionUID = 219047230178423923l;
    /**
     * Constructs a PDBParseException object.
     *
     * @param s  a String ...
     */
    
    public PDBParseException(String s) {
        super(s);
    }
    /**
     * Constructs a PDBParseException object.
     *
     * @param t  a Throwable object
     * @param s  a String ...
     */ 
    public PDBParseException ( String s,Throwable t) {
        super(s, t);
    }
    /**
     * Constructs a PDBParseException object.
     *
     * @param t  a Throwable object
     */
    public PDBParseException (Throwable t) {
        super(t);
    }
}
