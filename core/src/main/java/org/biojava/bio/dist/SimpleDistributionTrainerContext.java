/*
 *              BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *    http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *    http://www.biojava.org/
 *
 */


package org.biojava.bio.dist;

import java.io.Serializable;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * A no-frills implementation of DistributionTrainerContext.
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public class SimpleDistributionTrainerContext
implements DistributionTrainerContext, Serializable {
  private final Map distToTrainer;
  private final Set trainers;
  
  private double nullModelWeight;
  
  public double getNullModelWeight() {
    return this.nullModelWeight;
  }

  public void setNullModelWeight(double nullModelWeight) {
    this.nullModelWeight = nullModelWeight;
  }
  
  public void registerDistribution(Distribution dist) {
    if(!distToTrainer.keySet().contains(dist)) {
      dist.registerWithTrainer(this);
    }
  }
  
  public void registerTrainer(
    Distribution dist, DistributionTrainer trainer
  ) {
    distToTrainer.put(dist, trainer);
    trainers.add(trainer);
  }
  
  public DistributionTrainer getTrainer(Distribution dist) {
    return (DistributionTrainer) distToTrainer.get(dist);
  }

  public void addCount(Distribution dist, Symbol sym, double times)
  throws IllegalSymbolException {
    DistributionTrainer dt = getTrainer(dist);
    if(dt == null) {
      throw new NullPointerException(
        "No trainer associated with distribution " + dist
      );
    }
    if (sym instanceof AtomicSymbol) {
      dt.addCount(this, (AtomicSymbol) sym, times);
    } else {
//      Distribution nullModel = dist.getNullModel();
//      double totWeight = nullModel.getWeight(sym);
      for (
        Iterator asi = ((FiniteAlphabet) sym.getMatches()).iterator();
        asi.hasNext();
      ) {
        AtomicSymbol as = (AtomicSymbol) asi.next();
        //dt.addCount(this, as, times * (nullModel.getWeight(as) / totWeight));
        dt.addCount(this, as, times);
      }
    }
  }
  
  public double getCount(Distribution dist, Symbol sym)
  throws IllegalSymbolException {
    DistributionTrainer dt = getTrainer(dist);
    if(dt == null) {
      throw new NullPointerException(
        "No trainer associated with distribution " + dist
      );
    }
    if (sym instanceof AtomicSymbol) {
      return dt.getCount(this, (AtomicSymbol) sym);
    } else {
      double totWeight = 0.0;
      for (
        Iterator asi = ((FiniteAlphabet) sym.getMatches()).iterator();
        asi.hasNext();
      ) {
        AtomicSymbol as = (AtomicSymbol) asi.next();
        totWeight += dt.getCount(this, as);
      }
      return totWeight;
    }
  }
  
  public void train()
  throws ChangeVetoException {
    for(Iterator i = trainers.iterator(); i.hasNext(); ) {
      ((DistributionTrainer) i.next()).train(this, getNullModelWeight());
    }
  }
  
  public void clearCounts() {
    for(Iterator i = trainers.iterator(); i.hasNext(); ) {
      ((DistributionTrainer) i.next()).clearCounts(this);
    }
  }

  /**
   * Create a new context with no initial distributions or trainers.
   */
  public SimpleDistributionTrainerContext() {
    this.distToTrainer = new IdentityHashMap();
    this.trainers = new HashSet();
  }
}
