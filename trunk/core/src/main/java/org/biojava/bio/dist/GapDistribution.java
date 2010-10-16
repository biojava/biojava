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

import java.io.Serializable;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

/**
 * This distribution emits gap symbols.
 * <p>
 * It is a useful thing to have around for pair-wise alignment, as you can
 * build a PairDistribution that emits gaps in one sequence and Symbols in the
 * other. The GapDistriution will always emit with a probability of 1, as every
 * symbol has a matches alphabet that contains the empty set. Or is this so?
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public class GapDistribution
  extends
    Unchangeable
  implements
    Distribution,
    Serializable
{
  private final Alphabet alpha;
    private static final long serialVersionUID = 88622317;
    
  public double getWeight(Symbol sym) throws IllegalSymbolException {
    return 1.0;
  }

  public void setWeight(Symbol s, double w) throws IllegalSymbolException,
  UnsupportedOperationException {
    getAlphabet().validate(s);
    throw new UnsupportedOperationException(
      "The weights are immutable: " + s.getName() + " -> " + w
    );
  }

  public Alphabet getAlphabet() {
    return alpha;
  }
    
  public Symbol sampleSymbol() {
    return getAlphabet().getGapSymbol();
  }
  
  public Distribution getNullModel() {
    return this;
  }

  public void setNullModel(Distribution nullModel)
  throws IllegalAlphabetException, ChangeVetoException {
    throw new ChangeVetoException("Can't change null model for GapDistribution");
  }
  
  public void registerWithTrainer(DistributionTrainerContext dtc) {
    dtc.registerTrainer(this, IgnoreCountsTrainer.getInstance());
  }

  /**
   * Get a GapDistribution for an alphabet.
   *
   * @param alpha  the Alphabet that this distribution is over
   */
  public GapDistribution(Alphabet alpha) {
    this.alpha = alpha;
  }
}
  

