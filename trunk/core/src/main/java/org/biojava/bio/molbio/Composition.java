/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.bio.molbio;

import java.util.Iterator;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Computes composition statistics about a <code>SymbolList</code>.
 * Essentially a conveniece wrapper around a Distribution.
 * @author Mark Schreiber
 * @since 1.6
 */
public class Composition {
    private SymbolList symbolList;
    private Distribution distribution;
    private DistributionTrainerContext dtc;

    
    /**
     * Set the <code>SymbolList</code> to calculation the composition of.
     * @param symbolList a <code>SymbolList</code> from the DNA Alphabet.
     * @throws org.biojava.bio.symbol.IllegalSymbolException if <code>symbolList</code>
     * is not DNA.
     */
    public synchronized void setSymbolList(final SymbolList symbolList) throws IllegalSymbolException {
        this.symbolList = symbolList;
        train(symbolList);
    }

    private synchronized void train(final SymbolList symbolList) throws IllegalSymbolException {
        dtc.registerDistribution(getDistribution());
        for(Iterator i = symbolList.iterator(); i.hasNext();){
            dtc.addCount(getDistribution(), (Symbol)i.next(), 1.0);
        }
        dtc.train();
        dtc.clearCounts();
    }

    /**
     * Returns the distribution backing this class.
     * @return a <code>Distribution</code>
     */
    public Distribution getDistribution() {
        return distribution;
    }
}
