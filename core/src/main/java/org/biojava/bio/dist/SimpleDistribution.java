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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeAdapter;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * A simple implementation of a distribution, which works with any finite alphabet.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber
 * @since 1.0
 * @serial WARNING serialized versions of this class may not be compatible with later versions of BioJava
 */
public class SimpleDistribution
extends AbstractDistribution implements Serializable{
    static final long serialVersionUID = 7252850540926095728L;
    
    
  private transient AlphabetIndex indexer;
  private transient double[] weights = null;//because indexer is transient.
  private Distribution nullModel;
  private FiniteAlphabet alpha;
  
  private static class SymbolWeightMemento implements Serializable {
      static final long serialVersionUID = 5223128163879670657L;
      
      public final Symbol symbol;
      public final double weight;
      
      public SymbolWeightMemento(Symbol s, double weight) {
          this.symbol = s;
          this.weight = weight;
      }
  }
  
  private void writeObject(ObjectOutputStream oos)
      throws IOException
  {
      oos.defaultWriteObject();
      
      if (weights != null) {// fix for bug 2360
          SymbolWeightMemento[] swm = new SymbolWeightMemento[weights.length];
          for (int w = 0; w < swm.length; ++w) {
              swm[w] = new SymbolWeightMemento(indexer.symbolForIndex(w), weights[w]);
          }
          oos.writeObject(swm);
      }
  }

  private void readObject(ObjectInputStream stream)
    throws IOException, ClassNotFoundException
  {
    stream.defaultReadObject();
    
    //System.out.println("Alphabet for this dist is: "+alpha.getName());
    indexer = AlphabetManager.getAlphabetIndex(alpha);
    indexer.addChangeListener(
      new ChangeAdapter(){
        public void preChange(ChangeEvent ce) throws ChangeVetoException{
          if(hasWeights()){
            throw new ChangeVetoException(
              ce,
              "Can't allow the index to change as we have probabilities."
            );
          }
        }
      },AlphabetIndex.INDEX
    );
    weights = new double[alpha.size()];
    
    SymbolWeightMemento[] swm = (SymbolWeightMemento[]) stream.readObject();
    for (int m = 0; m < swm.length; ++m) {
        try {
            weights[indexer.indexForSymbol(swm[m].symbol)] = swm[m].weight;
        } catch (IllegalSymbolException ex) {
            throw new IOException("Symbol in serialized stream: "+swm[m].symbol.getName()+" can't be found in the alphabet");
        }
    }
  }

  public Alphabet getAlphabet() {
    return indexer.getAlphabet();
  }

  public Distribution getNullModel() {
    return this.nullModel;
  }



  protected void setNullModelImpl(Distribution nullModel)

  throws IllegalAlphabetException, ChangeVetoException {
    this.nullModel = nullModel;
  }


  /**
   * Indicate whether the weights array has been allocated yet.
   *
   * @return  true if the weights are allocated
   */
  protected boolean hasWeights() {
    return weights != null;
  }


  /**
   * Get the underlying array that stores the weights.
   *
   * <p>
   * Modifying this will modify the state of the distribution.
   * </p>
   *
   * @return  the weights array
   */
  protected double[] getWeights() {
    if(weights == null) {
      weights = new double[((FiniteAlphabet)getAlphabet()).size()];
      for(int i = 0; i < weights.length; i++) {
        weights[i] = Double.NaN;

      }
    }
     return weights;
  }



  public double getWeightImpl(AtomicSymbol s)

  throws IllegalSymbolException {
    if(!hasWeights()) {
      return Double.NaN;
    } else {
      int index = indexer.indexForSymbol(s);
      return weights[index];
    }
  }


  protected void setWeightImpl(AtomicSymbol s, double w)
  throws IllegalSymbolException, ChangeVetoException {
    double[] weights = getWeights();
    if(w < 0.0) {
      throw new IllegalArgumentException(
        "Can't set weight to negative score: " +
        s.getName() + " -> " + w
      );
    }
    weights[indexer.indexForSymbol(s)] = w;
  }

  private void initialise(FiniteAlphabet alphabet) {
    this.alpha = alphabet;
    this.indexer = AlphabetManager.getAlphabetIndex(alphabet);
    indexer.addChangeListener(
      new ChangeAdapter() {
        public void preChange(ChangeEvent ce) throws ChangeVetoException {
          if(hasWeights()) {
            throw new ChangeVetoException(
              ce,
              "Can't allow the index to change as we have probabilities."
            );
          }
        }
      },
      AlphabetIndex.INDEX
    );

    try {
      setNullModel(new UniformDistribution(alphabet));
    } catch (Exception e) {
      throw new BioError("This should never fail. Something is screwed!", e);
    }
  }

  /**
   * make an instance of SimpleDistribution for the specified Alphabet.
   */
  public SimpleDistribution(FiniteAlphabet alphabet)
  {
    initialise(alphabet);
  }

  /**
   * make an instance of SimpleDistribution with weights identical
   * to the specified Distribution.
   *
   * @param dist Distribution to copy the weights from.
   */
  public SimpleDistribution(Distribution dist)
  {
    try {
    initialise((FiniteAlphabet) dist.getAlphabet());

    // now copy over weights
    int alfaSize = ((FiniteAlphabet)getAlphabet()).size();

    for (int i = 0; i < alfaSize; i++) {
      weights = new double[alfaSize];
      weights[i] = dist.getWeight(indexer.symbolForIndex(i));
    }
    }
    catch (IllegalSymbolException ise) {
      System.err.println("an impossible error surely! "); ise.printStackTrace();
    }
  }

  /**
   * Register an SimpleDistribution.Trainer instance as the trainer for this distribution.
   */
  public void registerWithTrainer(DistributionTrainerContext dtc) {
   dtc.registerTrainer(this, new Trainer());
  }


  /**
   * A simple implementation of a trainer for this class.
   *
   * @author Matthew Pocock
   * @since 1.0
   */
  protected class Trainer implements DistributionTrainer {
    private final Count counts;

    /**
     * Create a new trainer.
     */
    public Trainer() {
      counts = new IndexedCount(indexer);
    }

    public void addCount(DistributionTrainerContext dtc, AtomicSymbol sym, double times)
    throws IllegalSymbolException {
      try {
          counts.increaseCount(sym, times);
      } catch (ChangeVetoException cve) {
        throw new BioError(
          "Assertion Failure: Change to Count object vetoed", cve
        );
      }
    }

    public double getCount(DistributionTrainerContext dtc, AtomicSymbol sym)
    throws IllegalSymbolException {
      return counts.getCount(sym);
    }



    public void clearCounts(DistributionTrainerContext dtc) {
      try {
        int size = ((FiniteAlphabet) counts.getAlphabet()).size();
        for(int i = 0; i < size; i++) {
          counts.zeroCounts();
        }
      } catch (ChangeVetoException cve) {
        throw new BioError(
          "Assertion Failure: Change to Count object vetoed",cve
        );
      }
    }



    public void train(DistributionTrainerContext dtc, double weight)
    throws ChangeVetoException {
      if(!hasListeners())  {
        trainImpl(dtc, weight);
      } else {
        ChangeSupport changeSupport = getChangeSupport(Distribution.WEIGHTS);
        synchronized(changeSupport) {
          ChangeEvent ce = new ChangeEvent(
            SimpleDistribution.this,
            Distribution.WEIGHTS
          );
          changeSupport.firePreChangeEvent(ce);
          trainImpl(dtc, weight);
          changeSupport.firePostChangeEvent(ce);
        }
      }
    }



    protected void trainImpl(DistributionTrainerContext dtc, double weight) {
      //System.out.println("Training");
      try {
        Distribution nullModel = getNullModel();
        double[] weights = getWeights();
        double[] total = new double[weights.length];
        double sum = 0.0;

        for(int i = 0; i < total.length; i++) {
          AtomicSymbol s = (AtomicSymbol) indexer.symbolForIndex(i);
          sum +=
            total[i] =
              getCount(dtc, s) +
              nullModel.getWeight(s) * weight;
        }
        double sum_inv = 1.0 / sum;
        for(int i = 0; i < total.length; i++) {
          //System.out.println("\t" + weights[i] + "\t" + total[i] * sum_inv);
          weights[i] = total[i] * sum_inv;
        }
      } catch (IllegalSymbolException ise) {
        throw new BioError(
          "Assertion Failure: Should be impossible to mess up the symbols.",ise
        );
      }
    }
  }
}



