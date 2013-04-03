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
 * An exception for use during translating amino acids in a PDB file.
 * @author Jules Jacobsen
 */
public class UnknownPdbAminoAcidException extends Exception {

    /**
	 * 
	 */
	private static final long serialVersionUID = -5571696240026118421L;
	/**
     * Constructs a PDBParseException object.
     *
     * @param s  a String ...
     */

    public UnknownPdbAminoAcidException(String s) {
        super(s);
    }
    /**
     * Constructs a UnknownPdbAminoAcidException object.
     *
     * @param t  a Throwable object
     * @param s  a String ...
     */
    public UnknownPdbAminoAcidException ( String s,Throwable t) {
        super(s, t);
    }
    /**
     * Constructs a UnknownPdbAminoAcidException object.
     *
     * @param t  a Throwable object
     */
    public UnknownPdbAminoAcidException (Throwable t) {
        super(t);
    }

}
