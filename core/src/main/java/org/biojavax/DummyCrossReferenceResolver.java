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
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.BioEntry;
import org.biojavax.bio.seq.InfinitelyAmbiguousSymbolList;


/**
 * A simple implementation of CrossReferenceResolver.
 * @author Richard Holland
 * @author Mark Schreiber
 * @since 1.5
 */
public class DummyCrossReferenceResolver implements CrossReferenceResolver {
        
    /**
     * {@inheritDoc}
     * All responses are instances of InfinitelyAmbiguousSymbolList.
     */
    public SymbolList getRemoteSymbolList(CrossRef cr, Alphabet a) {
        if (!(a instanceof FiniteAlphabet)) throw new IllegalArgumentException("Cannot construct dummy symbol list for a non-finite alphabet");
        return new InfinitelyAmbiguousSymbolList((FiniteAlphabet)a);
    }
    
    /**
     * {@inheritDoc}
     * All responses are null.
     */
    public BioEntry getRemoteBioEntry(CrossRef cr){
        return null;
    }
}

