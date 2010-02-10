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


package org.biojava.bio.dp.twohead;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MagicalState;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ListTools;

/**
 * Cache for columns of emission probabilities in pair-wise alignment
 * algorithms. 
 *
 * @author Matthew Pocock
 * @author David Huen (fixes for magical state)
 */
public class EmissionCache implements Serializable{
  private final Map eMap;
  private final Alphabet alpha;
  private final State[] states;
  private final int dsi;
  private final ScoreType scoreType;
  private final Symbol[] gap;
  
  public EmissionCache(
    Alphabet alpha,
    State[] states,
    int dsi,
    ScoreType scoreType
  ) {
    this.eMap = new HashMap();
    this.alpha = alpha;
    this.states = states;
    this.dsi = dsi;
    this.scoreType = scoreType;
    
    List alphas = alpha.getAlphabets();
    this.gap = new Symbol[alphas.size()];
    for(int i = 0; i < this.gap.length; i++) {
      this.gap[i] = ((Alphabet) alphas.get(i)).getGapSymbol();
    }
  }

    public final double [] getEmissions(List symList)
        throws IllegalSymbolException
    {
        return getEmissions(symList, true);
    }
 
    /**
     * Retrieve the emission scores from the cache for every EmissionState
     * for the specified symbols.  The scores
     *  will be computed and cached if not already available.
     * <p>
     * The correct behaviour should be to exclude emission from the
     * MagicalState except at the origin and bottom-right of the DP
     * matrix.  As such, the emission vector is only cached when it
     * it is computed for exorcise set to true. The vector is computed
     * afresh every time when exorcise is false.  it should be noted that
     * if exorcise is <b>NOT<b> set to false with the emission vector
     * at the ends of the matrix, the model will fail as the start and
     * and end transitions become impossible.
     *
     * @param symList a list of the symbols in each head that require a lookup for emission probabilities.
     * @param exorcise Prevents emission from the MagicalState.
     */
    public final double [] getEmissions(List symList, boolean exorcise)
        throws IllegalSymbolException 
    {
        double [] emission;
        if (exorcise) {
            emission = (double []) eMap.get(symList);
            if(emission == null) {
                emission = computeEmissions(symList, true);

                eMap.put(ListTools.createList(symList), emission);
            } else {
                //System.out.print("-");
            }
        }
        else {
            emission = computeEmissions(symList, false);
        }
        return emission;    
    }

    private double [] computeEmissions(List symList, boolean exorcise)
        throws IllegalSymbolException
    {
        // create a Symbol array to permit lookup of
        // the Symbol to query the emission Distributions with
        // for different values of advance for an EmissionState.
        Symbol sym[][] = new Symbol[2][2];
//        List ll = ListTools.createList(symList);
        sym[0][0] = AlphabetManager.getGapSymbol();
        sym[1][1] = alpha.getSymbol(Arrays.asList(new Symbol [] {
            (Symbol) symList.get(0),
            (Symbol) symList.get(1)
            }));
        sym[1][0] = alpha.getSymbol(Arrays.asList(new Symbol [] {
            (Symbol) symList.get(0),
            gap[1]
            }));
        sym[0][1] = alpha.getSymbol(Arrays.asList(new Symbol [] {
            gap[0],
            (Symbol) symList.get(1)
            }));

        double [] emission = new double[dsi];

        // compute the log(emission probability)
        for(int i = 0; i < dsi; i++) {
            if (exorcise && (states[i] instanceof MagicalState)) {
                // exclude emission from MagicalState
                emission[i] = Double.NEGATIVE_INFINITY;
            }
            else {
                EmissionState es = (EmissionState) states[i];
                int [] advance = es.getAdvance();
                Distribution dis = es.getDistribution();
                Symbol s = sym[advance[0]][advance[1]];
                emission[i] = Math.log(scoreType.calculateScore(dis, s));
            }
        }

        return emission;
    }

  public void clear() {
    eMap.clear();
  }
}
