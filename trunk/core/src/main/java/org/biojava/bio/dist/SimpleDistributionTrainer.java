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
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * A simple implemenation of a distribution trainer.
 * <p>
 * This requires the distribuiton being trained to have a working setWeight
 * method that doesn't throw an UnsupportedOperationExcepiton.
 * </p>
 *
 * @author Matthew Pocock
 * @deprecated  Distribution impls should be providing custom trainers.
 */

public final class SimpleDistributionTrainer
implements DistributionTrainer, Serializable {
  private final Distribution dis;
  private final Map c;

  {
    this.c = new HashMap();
  }

  public void addCount(
    DistributionTrainerContext dtc,
    AtomicSymbol sym,
    double count
  ) throws IllegalSymbolException {
    Double d = (Double) c.get(sym);
    if (d == null) {
      throw new IllegalSymbolException(
        "Symbol " + sym.getName() +
        " not found in " + dis.getAlphabet().getName()
      );
    }
    if(count < 0) {
      throw new Error(
        "Can't add a negative count for " + sym.getName() +
        " of " + count
      );
    }
    c.put(sym, new Double(d.doubleValue() + count));
  }

  public double getCount(
    DistributionTrainerContext dtc,
    AtomicSymbol sym
  ) throws IllegalSymbolException {
    Double d = (Double) c.get(sym);
    if (d == null) {
      throw new IllegalSymbolException(
        "Symbol " + sym.getName() +
        " not found in " + dis.getAlphabet().getName()
      );
    }
    return ((Double) c.get(sym)).doubleValue();
  }

  public void train(DistributionTrainerContext dtc, double weight)
  throws ChangeVetoException {
    Distribution nullModel = dis.getNullModel();
    double sum = 0.0;
    try {
      for(
        Iterator i = ((FiniteAlphabet) dis.getAlphabet()).iterator();
        i.hasNext();
      ) {
        Symbol s = (Symbol) i.next();
        Double d = (Double) c.get(s);
        sum += d.doubleValue() +
             nullModel.getWeight(s) * weight;
             // System.out.println(s.getName() + ": sum=" + sum);
      }
      for(
        Iterator i = ((FiniteAlphabet) dis.getAlphabet()).iterator();
        i.hasNext();
      ) {
        Symbol sym = (Symbol) i.next();
        Double d = (Double) c.get(sym);
        dis.setWeight(
          sym,
          (d.doubleValue() + nullModel.getWeight(sym) * weight) / sum
        );
      }
    } catch (IllegalSymbolException ise) {
      throw new BioError(
        "The alphabet for this distribution is not self-consistent"
      );
    }
  }

  public void clearCounts(DistributionTrainerContext dtc) {
    for(
      Iterator i = ((FiniteAlphabet) dis.getAlphabet()).iterator();
      i.hasNext();
    ) {
      c.put(i.next(), new Double(0.0));
    }
  }

  public SimpleDistributionTrainer(Distribution dis)
  throws IllegalAlphabetException {
    Alphabet a = dis.getAlphabet();
    if(! (a instanceof FiniteAlphabet)) {
      throw new IllegalAlphabetException(
        "Can't create a SimpleDistributionTrainer for non-finite alphabet " +
        a.getName() + " of type " + a.getClass()
      );
    }
    this.dis = dis;
  }
}
