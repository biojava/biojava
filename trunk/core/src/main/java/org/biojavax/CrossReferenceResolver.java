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

package org.biojavax;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.BioEntry;

/**
 * This interface returns symbols or sequence for a given cross-reference.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public interface CrossReferenceResolver {
    
    /**
     * Given a cross reference, return the corresponding symbol list.
     * @param cr the cross reference to look up.
     * @param a the alphabet to construct the infinitely ambiguous symbol list
     *        over if it cannot be found.
     * @return the symbol list matching it. If none, return an
     * infintely-ambiguous symbol list rather than null.
     */
    public SymbolList getRemoteSymbolList(CrossRef cr, Alphabet a);
    
    
    /**
     * Given the <code>CrossRef</code> return the corresponding
     * <code>BioEntry</code>
     * @param cr the cross reference
     * @throws org.biojavax.CrossReferenceResolutionException if a problem occurs
     *         during resolution.
     * @return The cross referenced entry (most likely a <code>RichSequence</code>
     *         or null if it cannot be found.
     */
    public BioEntry getRemoteBioEntry(CrossRef cr) 
    throws CrossReferenceResolutionException;
}
