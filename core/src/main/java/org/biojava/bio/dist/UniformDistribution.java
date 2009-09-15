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
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

/**
 * An implementation of an uniform distribution
 * @serial WARNING serialized versions of this class may be incompatible with future versions of BioJava
 * @author Matthew Pocock
 * @author Mark Schreiber
 * @author Thomas Down
 * @since 1.0
 */

public class UniformDistribution
  extends
    AbstractDistribution
  implements
    Serializable
{
  private final FiniteAlphabet alphabet;
  private Distribution nullModel;
  private static final long serialVersionUID = 88454038;

  public Alphabet getAlphabet() {
    return alphabet;
  }

  public Distribution getNullModel() {
    return nullModel;
  }

  /**
   * Assign a background distribution.
   *
   * @param nullModel the background distribution to assign
   * @throws IllegalAlphabetException if nullModel is over an incompattible
   *    alphabet
   */
  protected void setNullModelImpl(Distribution nullModel)
  throws IllegalAlphabetException {
    this.nullModel = nullModel;
  }

  protected double getWeightImpl(AtomicSymbol s)
  throws IllegalSymbolException {
    return (1.0 / alphabet.size());
  }

  protected void setWeightImpl(AtomicSymbol sym, double weight)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "Can't change the weights in a UniformDistribution"
    );
  }

  public void registerWithTrainer(DistributionTrainerContext dtc) {
    dtc.registerTrainer(this, IgnoreCountsTrainer.getInstance());
  }

  /**
   * Create a new UniformDistribution.
   *
   * @param alphabet  the finite alphabet to be over
   */
  public UniformDistribution(FiniteAlphabet alphabet) {
    this.alphabet = alphabet;
    this.nullModel = this;
  }
}
