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
 */


package org.biojava.bio.dist;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Mark Schreiber
 * @since 1.1
 */
class GeneralNthOrderDistribution extends AbstractOrderNDistribution implements Serializable{
    private Map dists;
    private static final long serialVersionUID = 42388921; //Change this value if internal implementation changes significantly

    GeneralNthOrderDistribution(Alphabet alpha, DistributionFactory df)
        throws IllegalAlphabetException
    {
        super(alpha);
        dists = new HashMap();

        for(Iterator i = ((FiniteAlphabet) getConditioningAlphabet()).iterator(); i.hasNext(); ) {
            Symbol si = (Symbol) i.next();
            dists.put(si.getName(), df.createDistribution(getConditionedAlphabet()));
        }
    }

    public Distribution getDistribution(Symbol sym)
        throws IllegalSymbolException
    {
        Distribution d = (Distribution) dists.get(sym.getName());
        if(d == null) {
            getConditioningAlphabet().validate(sym);
        }
        return d;
    }

    public void setDistribution(Symbol sym, Distribution dist)
        throws IllegalSymbolException, IllegalAlphabetException
    {
        getConditioningAlphabet().validate(sym);
        if(dist.getAlphabet() != getConditionedAlphabet()) {
            throw new IllegalAlphabetException(
                    "The distribution must be over " + getConditionedAlphabet() +
                    ", not " + dist.getAlphabet());
        }

        Distribution old = (Distribution) dists.get(sym);
        if( (old != null) && (weightForwarder != null) ) {
            old.removeChangeListener(weightForwarder);
        }

        if(weightForwarder != null) {
            dist.addChangeListener(weightForwarder);
        }

        dists.put(sym.getName(), dist);
    }

    public Collection conditionedDistributions() {
        return dists.values();
    }
}

