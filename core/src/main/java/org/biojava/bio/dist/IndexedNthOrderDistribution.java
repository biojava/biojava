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
import java.util.Arrays;
import java.util.Collection;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeListener;

/**
 *  Description of the Class
 *
 * @author     Thomas Down
 * @author     Matthew Pocock
 * @since      1.1
 */
class IndexedNthOrderDistribution
extends AbstractOrderNDistribution implements Serializable {
    static final long serialVersionUID = 8943232063125632673L;
    
    private transient Distribution[] dists;
    private transient AlphabetIndex index;


  IndexedNthOrderDistribution(Alphabet alpha, DistributionFactory df)
       throws IllegalAlphabetException {
    super(alpha);

    FiniteAlphabet conditioning = (FiniteAlphabet) getConditioningAlphabet();
    index = AlphabetManager.getAlphabetIndex(conditioning);
    index.addChangeListener(ChangeListener.ALWAYS_VETO, AlphabetIndex.INDEX);
    // Throws if alpha isn't indexable
    dists = new Distribution[conditioning.size()];

    for(int i = 0; i < conditioning.size(); ++i) {
      dists[i] = df.createDistribution(getConditionedAlphabet());
    }
  }

  /**
   *  Sets the Distribution attribute of the IndexedNthOrderDistribution
   *  object
   *
   * @param  sym                           a symbol in the conditioning alphabet
   * @param  dist                          the new Distribution
   * @exception  IllegalSymbolException    if the symbol is not in the
   *  conditioning alpahbet
   * @exception  IllegalAlphabetException  if the Distribution is not over the
   *   conditioned alphabet
   */
  public void setDistribution(Symbol sym, Distribution dist)
       throws IllegalSymbolException, IllegalAlphabetException {
    int indx = index.indexForSymbol(sym);
    if(dist.getAlphabet() != getConditionedAlphabet()) {
      throw new IllegalAlphabetException(
          "The distribution must be over " + getConditionedAlphabet() +
          ", not " + dist.getAlphabet());
    }

    Distribution old = dists[indx];
    if((old != null) && (weightForwarder != null)) {
      old.removeChangeListener(weightForwarder);
    }

    if(weightForwarder != null) {
      dist.addChangeListener(weightForwarder);
    }

    dists[indx] = dist;
  }

  /**
   *  Gets the Distribution attribute of the IndexedNthOrderDistribution
   *  object
   *
   * @param  sym                         a symbol
   * @return                             the Distribution for that symbol
   * @exception  IllegalSymbolException  if the symbol is not acceptable
   */
  public Distribution getDistribution(Symbol sym)
       throws IllegalSymbolException {
    return dists[index.indexForSymbol(sym)];
  }

  /**
   * Retrieve all of the distributions that are conditioned upon the
   * conditioning distributions.
   *
   * @return    a Collection of conditioned Distribution instances
   */
  public Collection conditionedDistributions() {
    //fixme: I'm sure this should return a List
    return Arrays.asList(dists);
  }


  private static class SymbolDistMemento implements Serializable {
      static final long serialVersionUID = 8387939181177306749L;
      
      public final Symbol symbol;
      public final Distribution dist;
      
      public SymbolDistMemento(Symbol s, Distribution dist) {
          this.symbol = s;
          this.dist = dist;
      }
  }
  
  private void writeObject(ObjectOutputStream oos)
      throws IOException
  {
      oos.defaultWriteObject();
      
      SymbolDistMemento[] sdm = new SymbolDistMemento[dists.length];
      for (int w = 0; w < sdm.length; ++w) {
          sdm[w] = new SymbolDistMemento(index.symbolForIndex(w), dists[w]);
      }
      oos.writeObject(sdm);
  }

  private void readObject(ObjectInputStream stream)
    throws IOException, ClassNotFoundException
  {
    stream.defaultReadObject();
    
    FiniteAlphabet conditioning = (FiniteAlphabet) getConditioningAlphabet();
    index = AlphabetManager.getAlphabetIndex(conditioning);
    index.addChangeListener(ChangeListener.ALWAYS_VETO, AlphabetIndex.INDEX);
    dists = new Distribution[conditioning.size()];
    
    SymbolDistMemento[] sdm = (SymbolDistMemento[]) stream.readObject();
    for (int m = 0; m < sdm.length; ++m) {
        try {
            dists[index.indexForSymbol(sdm[m].symbol)] = sdm[m].dist;
        } catch (IllegalSymbolException ex) {
            throw new IOException("Symbol in serialized stream can't be found in the alphabet");
        }
    }
  }
}
