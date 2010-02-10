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


package org.biojava.bio.dist;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * Default factory for Order-N distributions.
 *
 * @author Thomas Down
 * @author Mark Schreiber
 * @since 1.1
 */

public class OrderNDistributionFactory implements DistributionFactory {
    /**
     * Factory which used DistributionFactory.DEFAULT to create conditioned
     * distributions.
     */

    public static final DistributionFactory DEFAULT;

    static {
        DEFAULT = new OrderNDistributionFactory(DistributionFactory.DEFAULT);
    }

    private final DistributionFactory df;

    /**
     * Construct a new OrderNDistributionFactory with a specified factory
     * for conditioned distributions.
     *
     * @param df The DistributionFactory used for construction new conditioned
     *           distributions.
     */

    public OrderNDistributionFactory(DistributionFactory df) {
        this.df = df;
    }


    /**
     * Creates an OrderNDistribution of the appropriate type.
     *
     * @param alpha the Alphabet should be in a form that clearly indicates the
     * conditioning and the conditioned alphabet unless it is very obvious. For
     * example (DNA x DNA) is obvious, ((DNA x DNA x DNA) x DNA) indicates that
     * (DNA x DNA x DNA) is the conditioning <code>Alphabet</code> and DNA is the
     * conditioned <code>Alphabet</code>. (DNA x DNA x DNA x DNA) doesn't but
     * for compatibility with biojava 1.2 this is allowed in the constructor.
     * As from biojava 1.2.3 or greater this will be internally converted to
     * ((DNA x DNA x DNA) x DNA) which was the convention implied by biojava 1.2
     * Calls to the returned <code>Distribution</code>s <code>getAlphabet()</code>
     * method will return the converted <code>Alphabet</code>.
     *
     * @return An OrderNDistribution
     * @throws IllegalAlphabetException if a Distribution cannot be made with
     * that <code>Alphabet</code>.
     */
    public Distribution createDistribution(Alphabet alpha)
        throws IllegalAlphabetException
    {
        List aList = alpha.getAlphabets();
        if (
          aList.size() == 2 &&
          aList.get(0) == org.biojava.bio.seq.DNATools.getDNA()
        ) {
            return new IndexedNthOrderDistribution(alpha, df);
        } else {
            //convert things like (DNA x DNA x DNA) to ((DNA x DNA) x DNA)
            Alphabet conditioned = (Alphabet)aList.get(aList.size()-1);
            Alphabet conditioning =
                AlphabetManager.getCrossProductAlphabet(aList.subList(0,aList.size()-1));
            List l = new ArrayList();
            l.add(conditioning);
            l.add(conditioned);
            alpha = AlphabetManager.getCrossProductAlphabet(l);
            //System.out.println(alpha.getName());
            return new GeneralNthOrderDistribution(alpha, df);
        }
    }
}
