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

import java.util.Collection;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * Provides an N'th order distribution.  This is a distribution over one
 * alphabet which is conditioned on having previously observed one or
 * more other symbols (potentially from different alphabets).
 *
 * <p>
 * Order-N distributions are always over a CrossProductAlphabet.
 * </p>
 *
 * <p>
 * <strong>Note:</strong> Unlike normal distributions, the total weights for
 * all symbols in the overall alphabet do <em>not</em> sum to 1.0.  Instead,
 * the weights of each sub-distribution should sum to 1.0.
 * </p>
 *
 * <p>
 * This would typically be used in conjunction with an OrderNSymbolList.
 * </p>
 *
 * @author Thomas Down
 * @author Samiul Hasan
 * @author Matthew Pocock
 * @since 1.0
 */

public interface OrderNDistribution extends Distribution {
    /**
     * Get the conditioning alphabet of this distribution.  If the `overall'
     * alphabet is a cross-product of two alphabets, this will be the first 
     * of those alphabets.  If it is a cross-product of more than two alphabets,
     * the conditioning alphabet is the cross-product of all but the last
     * alphabet.
     *
     * @return the conditioning Alphabet
     */

    public Alphabet getConditioningAlphabet();

    /**
     * Get the conditioned alphabet.  This is the last alphabet in the
     * distribution's overall cross-product.  It will be the alphabet of
     * all the sub-distributions contained within this OrderNDistribution.
     *
     * @return the conditioned Alphabet
     */

    public Alphabet getConditionedAlphabet();

  /**
   * Get the conditioned distributions.
   *
   * @return  the conditioned distributions
   */
  //fixme: I'm sure this should return a list
    public Collection conditionedDistributions();

  /**
   * Set the distribution assocated with a symbol.
   *
   * @param sym   the symbol in the conditioning Alphabet
   * @param dist  a distribution over the conditioned Alphabet
   * @throws IllegalSymbolException   if sym is not in the conditioning Alphabet
   * @throws IllegalAlphabetException if dist is not over the conditioned
   *    Alphabet
   */
      public abstract void setDistribution(Symbol sym, Distribution dist)
      throws IllegalSymbolException, IllegalAlphabetException;
  
  public abstract Distribution getDistribution(Symbol sym)
      throws IllegalSymbolException;
}
