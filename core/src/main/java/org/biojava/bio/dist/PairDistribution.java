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
import java.lang.ref.SoftReference;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ListTools;

/**
 * Class for pairing up two independant distributions.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Samiul Hasan
 * @since 1.1
 */


public class PairDistribution
extends AbstractChangeable
implements Serializable, Distribution {
  private static Map cache;

  static {
    cache = new HashMap();
  }

  /**
   * Get a uniform null model over a PairDistribution over [first,second].
   *
   * @param first   the first Alphabet
   * @param second  the second Alphabet
   * @return    a Distribution that is a uniform distribution over the product
   *    of first and second
   */
  protected static Distribution getNullModel(Distribution first, Distribution second) {
    synchronized(cache) {
      first = first.getNullModel();
      second = second.getNullModel();
      List distL = new ListTools.Doublet(first, second);
      SoftReference ref = (SoftReference) cache.get(distL);
      Distribution dist;
      if(ref == null) {
        dist = new PairDistribution(first, second);
        cache.put(distL, new SoftReference(dist));
      } else {
        dist = (Distribution) ref.get();
        if(dist == null) {
          dist = new PairDistribution(first, second);
          cache.put(distL, new SoftReference(dist));
        }
      }
      return dist;
    }
  }

  private Distribution first;
  private Distribution second;
  private Alphabet alphabet;

  public Alphabet getAlphabet() {
    return alphabet;
  }

  public Distribution getNullModel() {
    return getNullModel(first, second);
  }

  public void setNullModel(Distribution nullModel)
  throws IllegalAlphabetException, ChangeVetoException {
    throw new ChangeVetoException(
      "PairDistribution objects can't have their null models changed."
    );
  }

  /**
   * Register this paired distribution with a model trainer.
   * @param trainer the trainer to register this distribution with.
   */
  public void registerWithTrainer(org.biojava.bio.dp.ModelTrainer trainer) {
    trainer.registerDistribution(first);
    trainer.registerDistribution(second);

    trainer.registerTrainer(this, new PairTrainer());
  }

  public double getWeight(Symbol sym)
  throws IllegalSymbolException {
    if(sym instanceof BasisSymbol) {
      List symL = ((BasisSymbol) sym).getSymbols();
      Symbol f = (Symbol) symL.get(0);
      Symbol s = (Symbol) symL.get(1);

      return first.getWeight(f) * second.getWeight(s);
    } else {
      double score = 0.0;
      for(Iterator i = ((FiniteAlphabet) sym.getMatches()).iterator();
      i.hasNext(); ) {
        AtomicSymbol s = (AtomicSymbol) i.next();
        score += getWeight(s);
      }
      return score;
    }
  }

  public void setWeight(Symbol sym, double weight)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "Can't set the weight directly in a PairDistribution. " +
      "You must set the weights in the underlying distributions."
    );
  }

  /**
   * Create a new PairDistribution that represents the product of two other
   * distributions. The alphabet will be the product of the first and seccond
   * distribution's alphabets, and the weights will be the products of the
   * weights for the first and seccond distributions given the first and second
   * component of the symbol respectively.
   *
   * @param first   the first distribution
   * @param second  the second distribution
   */
  public PairDistribution(Distribution first, Distribution second) {
    this.first = first;
    this.second = second;
    this.alphabet = AlphabetManager.getCrossProductAlphabet(
      Arrays.asList(new Alphabet[] {
        first.getAlphabet(), second.getAlphabet()
      })
    );
  }

  public void registerWithTrainer(DistributionTrainerContext dtc) {
    dtc.registerTrainer(this, new PairTrainer());
  }

  private class PairTrainer
  extends IgnoreCountsTrainer
  implements Serializable {
    public double getCount(DistributionTrainerContext dtc, AtomicSymbol as)
    throws IllegalSymbolException {
      getAlphabet().validate(as);

      List symL = as.getSymbols();
      Symbol f = (Symbol) symL.get(0);
      Symbol s = (Symbol) symL.get(1);

      // I don't think this is correct. Pants!
      return
        (dtc.getCount(first, f) + dtc.getCount(second, s)) * 0.5;

    }

    public void addCount(
      DistributionTrainerContext dtc, Symbol sym, double times
    ) throws IllegalSymbolException {
      getAlphabet().validate(sym);
      if(!(sym instanceof AtomicSymbol)) {
        throw new IllegalSymbolException(
          "Can't add counts for ambiguity symbols. Got: " +
          sym.getName()
        );
      }
      // FIXME: should get matches for symbol &
      // divide count by null model ratioes.
      List symL = ((BasisSymbol) sym).getSymbols();
      Symbol f = (Symbol) symL.get(0);
      Symbol s = (Symbol) symL.get(1);

      dtc.addCount(first, f, times);
      dtc.addCount(second, s, times);
    }
  }

  public Symbol sampleSymbol() {
    try {
      return getAlphabet().getSymbol(Arrays.asList( new Symbol[] {
        first.sampleSymbol(),
        second.sampleSymbol()
      }));
    } catch (IllegalSymbolException ise) {
      throw new BioError("Couldn't sample symbol", ise);
    }
  }
}
