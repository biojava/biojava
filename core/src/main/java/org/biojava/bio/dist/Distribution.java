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

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <p>
 * An encapsulation of a probability distribution over the Symbols within an
 * alphabet.
 * </p>
 *
 * <p>
 * A distribution can be implemented as a map from symbol to probability. It is
 * more correct to think of them as being integrals or sums over probability
 * dencity funcitons. In this world view, getWeight should look at the
 * getMatches of the symbol it is given and then perform the apropreate sum or
 * integral to return the probability of something within that set of symbols
 * being emitted.
 * </p>
 *
 * <p>
 * This interface should handle the case of emitting an ambiguity symbol.
 * This should be just the sum of the probabiltiy of emitting each matching
 * symbol. It is up to the code using the Distribution instance to divide out
 * the null model appropreately.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public interface Distribution extends Changeable {
  /**
   * <p>
   * Whenever a distribution changes the values that would be returned by
   * getWeight, they should fire a ChangeEvent with this object as the type.
   * </p>
   *
   * <p>
   * If the whole distribution changes, then the change and previous fields of
   * the ChangeEvent should be left null. If only a single weight is modified,
   * then change should be of the form Object[] { symbol, new Double(newVal) }
   * and previous should be of the form Object[] { symbol, new Double(oldVal) }
   * </p>
   */
  public static final ChangeType WEIGHTS = new ChangeType(
    "distribution weights changed",
    "org.biojava.bio.dist.Distribution",
    "WEIGHTS"
  );

  /**
   * <p>
   * Whenever the null model distribution changes the values that would be
   * returned by getWeight, either by being edited or by being replaced, a
   * ChangeEvent with this object as the type should be thrown.
   * </p>
   *
   * <p>
   * If the null model has changed its weights, then the ChangeEvent should
   * refer back to the ChangeEvent from the null model.
   * </p>
   */
  public static final ChangeType NULL_MODEL = new ChangeType(
    "distribution null model changed",
    "org.biojava.bio.dist.Distribution",
    "NULL_MODEL"
  );
  
  /**
   * The alphabet from which this spectrum emits symbols.
   *
   * @return  the Alphabet associated with this spectrum
   */
  Alphabet getAlphabet();
    
  /**
   * <p>
   * Return the probability that Symbol s is emitted by this spectrum.
   * </p>
   *
   * <p>
   * If the symbol is  ambiguou, then it is the sum of the probability that
   * each one of the matching symbols was emitted.
   * </p>
   *
   * @param s the Symbol emitted
   * @return  the probability of emitting that symbol
   * @throws IllegalSymbolException if s is not from this state's alphabet
   */
  double getWeight(Symbol s) throws IllegalSymbolException;
  
  /**
   * Set the probability or odds that Symbol s is emitted by this state.
   *
   * @param s the Symbol emitted
   * @param w  the probability of emitting that symbol
   * @throws IllegalSymbolException if s is not from this state's alphabet, or
   *         if it is an ambiguity symbol and the implementation can't handle
   *         this case
   * @throws ChangeVetoException if this state does not allow weights
   *         to be tampered with, or if one of the listeners vetoed this change
   */
  void setWeight(Symbol s, double w)
  throws IllegalSymbolException, ChangeVetoException;

  /**
   * Sample a symbol from this state's probability distribution.
   *
   * @return the symbol sampled
   */
  Symbol sampleSymbol();
  
  /**
   * Retrieve the null model Distribution that this Distribution recognizes.
   *
   * @return  the apropriate null model
   */
  Distribution getNullModel();
  
  /**
   * Set the null model Distribution that this Distribution recognizes.
   *
   * @param nullDist  the new null model Distribution
   * @throws IllegalAlphabetException if the null model has the wrong alphabet
   * @throws ChangeVetoException  if this Distirbution doesn't support setting
   *         the null model, or if one of its listeners objects
   */
  void setNullModel(Distribution nullDist)
  throws IllegalAlphabetException, ChangeVetoException;
  
  /**
   * <p>
   * Register this distribution with a training context.
   * </p>
   *
   * <p>
   * This should be invoked from within dtc.addDistribution(). This method
   * is responsible for constructing a suitable DistributionTrainer instance
   * and registering it by calling
   * dtc.registerDistributionTrainer(this, trainer). If the distribution is a
   * view onto another distribution, it can force the other to be registered by
   * calling dtc.addDistribution(other), and can then get on with registering
   * it's own trainer.
   * </p>
   *
   * @param dtc the DistributionTrainerContext with witch to register a trainer
   */
  void registerWithTrainer(DistributionTrainerContext dtc);
  
  /**
   * This listens to the null model distribution events and converts them into
   * NULL_MODEL events.
   *
   * @author Matthew Pocock
   * @since 1.1
   * @deprecated use
   *    <code>new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY)</code>
   *    instead
   */
  public class NullModelForwarder extends ChangeForwarder {
    /**
     * Create a new forwarder.
     *
     * @param source  the Object who events are forwarded on behalf of
     * @param cs      the change support that manages the listeners
     */
    public NullModelForwarder(Object source, ChangeSupport cs) {
      super(source, cs);
    }
    
    protected ChangeEvent generateEvent(ChangeEvent ce) {
      if(ce.getType() == WEIGHTS) {
        return new ChangeEvent(
          getSource(),
          NULL_MODEL,
          null, null,
          ce
        );
      }
      return null;
    }
  }
}
