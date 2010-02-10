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

import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * An abstract implementation of Distribution.
 * <p>
 * You will need to override <code>getWeight()</code> for a simple
 * implementation. You may also wish to override the other methods if the
 * default implementation is not suitable.
 * </p>
 *
 * <p>
 * The <code>registerWithTrainer</code> method registers
 * an <code>IgnoreCountsTrainer</code>.  To make an <code>AbstractDistribution</code>
 * subclass trainable, this method must be overridden.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber (serialization support)
 * @author Greg Cox
 * @since 1.0
 */

public abstract class AbstractDistribution
  extends
    AbstractChangeable
  implements
    Distribution
{
  /**
   * Forwarder for modifications to the null model.
   */
  protected transient ChangeForwarder nullModelForwarder = null;

  protected ChangeSupport getChangeSupport(ChangeType ct) {
      ChangeSupport changeSupport = super.getChangeSupport(ct);

    if(
      ((Distribution.NULL_MODEL.isMatchingType(ct)) || (ct.isMatchingType(Distribution.NULL_MODEL))) &&
      nullModelForwarder == null
    ) {
      nullModelForwarder =
              new ChangeForwarder.Retyper(this, changeSupport, Distribution.NULL_MODEL);
      getNullModel().addChangeListener(nullModelForwarder, Distribution.WEIGHTS);
    }

    return changeSupport;
  }

  /**
   * Implement this to actually set the weight.
   *
   * <p>
   * Do not inform any listeners. This has already been done for you. Just
   * update state.
   * </p>
   *
   * @param sym     the AtomicSymbol to update for
   * @param weight  the new weight for that symbol
   * @throws IllegalSymbolException if the symbol is not known
   * @throws ChangeVetoException    if the change is to be prevented
   */
  abstract protected void setWeightImpl(AtomicSymbol sym, double weight)
  throws IllegalSymbolException, ChangeVetoException;

  /**
   * Set the weight of a given symbol in this distribution.
   * <P>
   * This implementation informs all listeners of the change, and then calls
   * setWeightImpl to make the actual change. Sub-classes should over-ride
   * setWeightImpl to implement the actual storage of the weights.
   *
   * @param sym  the Symbol to set the weight for
   * @param weight  it's new weight
   * @throws IllegalSymbolException if sym is not known
   * @throws ChangeVetoException    if the update was prevented
   */
  final public void setWeight(Symbol sym, double weight)
  throws IllegalSymbolException, ChangeVetoException {
    if(!hasListeners()) {
      doSetWeight(sym, weight);
    } else {
      ChangeEvent ce = new ChangeEvent(
        this,
        Distribution.WEIGHTS,
        new Object[] {sym, new Double(weight)},
        new Object[] {sym, new Double(getWeight(sym))}
      );
      ChangeSupport changeSupport = super.getChangeSupport(Distribution.WEIGHTS);
      synchronized(changeSupport) {
        changeSupport.firePreChangeEvent(ce);
        doSetWeight(sym, weight);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  private void doSetWeight(Symbol sym, double weight)
  throws IllegalSymbolException, ChangeVetoException {
    if(sym instanceof AtomicSymbol) {
      setWeightImpl((AtomicSymbol) sym, weight);
    } else {
      //need to divide the weight up amongst the atomic symbols according
      //to the null model
      FiniteAlphabet fa = (FiniteAlphabet) sym.getMatches();
      double totalNullWeight = this.getNullModel().getWeight(sym);
      for(Iterator si = fa.iterator(); si.hasNext(); ) {
        AtomicSymbol as = (AtomicSymbol) si.next();
        double symNullWeight = this.getNullModel().getWeight(as);
        setWeightImpl(as, weight * symNullWeight / totalNullWeight);
      }
    }
  }

  /**
   * Implement this to set the null model.
   *
   * <p>
   * You should not inform any change listeners in this method. All of that
   * work has been done for you.
   * </p>
   *
   * @param nullModel  the new null model Distribution
   * @throws IllegalAlphabetException if the null model is for the wrong alphabet
   * @throws ChangeVetoException  if your implementation wishes to block this
   *    opperation
   */
  abstract protected void setNullModelImpl(Distribution nullModel)
  throws IllegalAlphabetException, ChangeVetoException;

  final public void setNullModel(Distribution nullModel)
  throws IllegalAlphabetException, ChangeVetoException {
    if(nullModel.getAlphabet() != getAlphabet()) {
      throw new IllegalAlphabetException(
        "Could not use distribution " + nullModel +
        " as its alphabet is " + nullModel.getAlphabet().getName() +
        " and this distribution's alphabet is " + getAlphabet().getName()
      );
    }
    Distribution oldModel = getNullModel();
    if(nullModelForwarder != null) {
      if(oldModel != null) {
        oldModel.removeChangeListener(nullModelForwarder);
      }
      nullModel.addChangeListener(nullModelForwarder);
    }
    if(!hasListeners()) {
      // if there are no listeners yet, don't go through the overhead of
      // synchronized regions or of trying to inform them.
      setNullModelImpl(nullModel);
    } else {
      // OK - so somebody is intereted in me. Do it properly this time.
      ChangeEvent ce = new ChangeEvent(
        this,
        Distribution.NULL_MODEL,
        nullModel,
        oldModel
      );
      ChangeSupport changeSupport = super.getChangeSupport(Distribution.NULL_MODEL);
      synchronized(changeSupport) {
        changeSupport.firePreChangeEvent(ce);
        setNullModelImpl(nullModel);
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  /**
   * Retrieve the weight for this distribution.
   * <P>
   * Performs the standard munge to handle ambiguity symbols. The actual weights
   * for each atomic symbol should be calculated by the getWeightImpl
   * functions.
   *
   * @param sym the Symbol to find the probability of
   * @return the probability that one of the symbols matching amb was emitted
   * @throws IllegalSymbolException if for any reason the symbols within amb
   *         are not recognized by this state
   */
  public final double getWeight(Symbol sym)
  throws IllegalSymbolException {
    if(sym instanceof AtomicSymbol) {
      return getWeightImpl((AtomicSymbol) sym);
    } else {
      Alphabet ambA = sym.getMatches();
      if(((FiniteAlphabet) ambA).size() == 0) { // a gap
        getAlphabet().validate(sym);

        double totalWeight = 0.0;
        for (Iterator i = ((FiniteAlphabet)getAlphabet()).iterator();
                          i.hasNext(); ) {

          Symbol s = (Symbol)i.next();
          totalWeight += getWeight(s);
        }
        return 1.0 - totalWeight;
      }
      if(ambA instanceof FiniteAlphabet) {
        FiniteAlphabet fa = (FiniteAlphabet) ambA;
        double sum = 0.0;
        for(Iterator i = fa.iterator(); i.hasNext(); ) {
          Object obj = i.next();
          if(!(obj instanceof AtomicSymbol)) {
            throw new BioError(
              "Assertion Failure: Not an instance of AtomicSymbol: " +
              obj
            );
          }
          AtomicSymbol as = (AtomicSymbol) obj;
          sum += getWeightImpl(as);
        }
        return sum;
      } else {
        throw new IllegalSymbolException(
           "Can't find weight for infinite set of symbols matched by " +
           sym.getName()
        );
      }
    }
  }

  /**
   * Override this method to implement getting the weight for an atomic
   * symbol. You should just do what is necessary to fetch state. All the work
   * with exceptions and listeners will have been handled for you.
   *
   * @param sym   the AtomicSymbol to get the weight for
   * @return      the weight
   * @throws IllegalSymbolException if sym is not known
   */
  protected abstract double getWeightImpl(AtomicSymbol sym)
  throws IllegalSymbolException;

  public Symbol sampleSymbol() {
    double p = Math.random();
    try {
      for(Iterator i = ((FiniteAlphabet) getAlphabet()).iterator(); i.hasNext(); ) {
        AtomicSymbol s = (AtomicSymbol) i.next();
        p -= getWeight(s);
        if( p <= 0) {
          return s;
        }
      }
      return getAlphabet().getGapSymbol();

//      StringBuffer sb = new StringBuffer();
//      for(Iterator i = ((FiniteAlphabet) this.getAlphabet()).iterator(); i.hasNext(); ) {
//        AtomicSymbol s = (AtomicSymbol) i.next();
//        double w = getWeight(s);
//        sb.append("\t" + s.getName() + " -> " + w + "\n");
//      }
//      throw new BioError(
//        "Could not find a symbol to emit from alphabet " + getAlphabet() +
//        ". Do the probabilities sum to 1?" + "\np=" + p + "\n" + sb.substring(0)
//      );
    } catch (IllegalSymbolException ire) {
      throw new BioError(
        "Unable to iterate over all symbols in alphabet - " +
        "things changed beneath me!", ire
      );
    }
  }

  /**
   * Register an IgnoreCountsTrainer instance as the trainer for this
   * distribution.  Override this if you wish to implement a trainable
   * distribution.
   *
   * @param dtc  the context to register with
   */
  public void registerWithTrainer(DistributionTrainerContext dtc) {
    dtc.registerTrainer(this, IgnoreCountsTrainer.getInstance());
  }
  
  public int hashCode() {
      int hc = 0;
      try {
          for (Iterator i = ((FiniteAlphabet) getAlphabet()).iterator(); i.hasNext(); ) {
              Symbol s = (Symbol) i.next();
              hc = hc ^ (int) (getWeight(s) * (1 << 30));
          }
      } catch (IllegalSymbolException ex) {
          throw new BioError("Assertion failed", ex);
      }
      return hc;
  }
  
  public boolean equals(Object o) {
      if (o instanceof Distribution) {
          Distribution d = (Distribution) o;
          if (d.getAlphabet() != this.getAlphabet()) {
              return false;
          }
          try {
              for (Iterator i = ((FiniteAlphabet) getAlphabet()).iterator(); i.hasNext(); ) {
                  Symbol s = (Symbol) i.next();
                  if (this.getWeight(s) != d.getWeight(s)) {
                      return false;
                  }
              }
              return true;
          } catch (IllegalSymbolException ex) {
              throw new BioError(ex);
          }
      }
      return false;
  }
}
