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
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Simple base class for OrderNDistributions.
 *
 * @author Samiul Hasan
 * @author Matthew Pocock
 * @author Thomas Down
 * @since 1.2
 */
public abstract class AbstractOrderNDistribution extends AbstractDistribution
    implements OrderNDistribution, Serializable
{
    static final long serialVersionUID = 1406135308618188893L;
    
  private Alphabet alphabet;
  private Alphabet firstA;
  private Alphabet lastA;
  private Distribution nullModel;

  /**
   * The listener that will forward events from the underlying distributions to
   * listeners for this distribution.
   */
  protected transient ChangeForwarder weightForwarder = null;
  
  protected ChangeSupport getChangeSupport(ChangeType ct) {
    ChangeSupport changeSupport = super.getChangeSupport(ct);
    
    if(
      ( (Distribution.WEIGHTS.isMatchingType(ct)) || (ct.isMatchingType(Distribution.WEIGHTS)) ) &&
      weightForwarder == null
    ) {
      weightForwarder =
              new ChangeForwarder.Retyper(this, changeSupport, Distribution.WEIGHTS);
      for(Iterator i = conditionedDistributions().iterator(); i.hasNext(); ) {
        Distribution dist = (Distribution) i.next();
        dist.addChangeListener(weightForwarder, Distribution.WEIGHTS);
      }
    }
    
    return changeSupport;
  }
  
    /**
     * Construct a new NthOrderDistribution.
     *
     * @param alpha  the Alpahbet this is over
     */

  protected AbstractOrderNDistribution(Alphabet alpha)
  throws IllegalAlphabetException  {
    this.alphabet = alpha;
    List aList = alpha.getAlphabets();
    int lb1 = aList.size() - 1;
    if(aList.size() == 2) {
      this.firstA = (Alphabet) aList.get(0);
    } else {
      this.firstA = AlphabetManager.getCrossProductAlphabet(aList.subList(0, lb1));
    }
    this.lastA = (Alphabet) aList.get(lb1); 
    this.nullModel = new UniformNullModel();
  }
  
    /**
     * Get the conditioning alphabet of this distribution.  If the `overall'
     * alphabet is a cross-product of two alphabets, this will be the first 
     * of those alphabets.  If it is a cross-product of more than two alphabets,
     * the conditioning alphabet is the cross-product of all but the last
     * alphabet.
     *
     * @return the conditioning Alphabet
     */

    public Alphabet getConditioningAlphabet() {
	return firstA;
    }

    /**
     * Get the conditioned alphabet.  This is the last alphabet in the
     * distribution's overall cross-product.  It will be the alphabet of
     * all the sub-distributions contained within this OrderNDistribution.
     */

    public Alphabet getConditionedAlphabet() {
	return lastA;
    }
   
  public Alphabet getAlphabet() {
    return alphabet;
  }
  
    /**
     * Get a weight from one of the sub-distributions, conditioned
     * on the first part of the symbol.
     *
     * @param sym  the symbol to look up
     * @return the weight
     * @throws IllegalSymbolException  if sym is not recognised
     */

  protected double getWeightImpl(AtomicSymbol sym) throws IllegalSymbolException {
    List symL = sym.getSymbols();
    int lb1 = symL.size() - 1;
    BasisSymbol firstS;
    if(symL.size() == 2) {
      firstS = (AtomicSymbol) symL.get(0);
    } else {
      firstS = (AtomicSymbol) firstA.getSymbol(symL.subList(0, lb1));
    }
    Distribution dist = getDistribution(firstS);
    return dist.getWeight((AtomicSymbol) symL.get(lb1));
  }

    /**
     * Set a weight in one of the conditioned distributions.  It is the callers
     * responsibility to ensure that all the conditioned distributions have total
     * weights which sum to 1.0.
     *
     * @param sym   the symbol to set the weight for
     * @param w     the new weight
     */

    public void setWeightImpl(AtomicSymbol sym, double w) 
    throws IllegalSymbolException, ChangeVetoException {
      List symL = sym.getSymbols();
      int lb1 = symL.size() - 1;
      Symbol firstS;
      if(symL.size() == 2) {
        firstS = (Symbol) symL.get(0);
      } else {
        firstS = firstA.getSymbol(symL.subList(0, lb1));
      }
      Distribution dist = getDistribution(firstS);
      dist.setWeight((Symbol) symL.get(lb1), w);
    }
  
  public void setNullModelImpl(Distribution nullModel) {
  	this.nullModel = nullModel;   
  }
  
  public Distribution getNullModel()  {
  	return this.nullModel;
  }
  
  public void registerWithTrainer(DistributionTrainerContext dtc) {
    for(Iterator i = conditionedDistributions().iterator(); i.hasNext(); ) {
      dtc.registerDistribution((Distribution) i.next());
    }
    dtc.registerTrainer(this, new IgnoreCountsTrainer() {
      public void addCount(
        DistributionTrainerContext dtc,
        AtomicSymbol sym,
        double count
      ) throws IllegalSymbolException {
        List symL = sym.getSymbols();
        int lb1 = symL.size() - 1;
        Symbol firstS;
        if(lb1 == 1) {
          firstS = (Symbol) symL.get(0);
        } else {
          firstS = firstA.getSymbol(symL.subList(0, lb1));
        }
        Distribution dist = getDistribution(firstS);
        dtc.addCount(dist, (Symbol) symL.get(lb1), count);
      }
      
      public double getCount(
        DistributionTrainerContext dtc,
        AtomicSymbol sym
      ) throws IllegalSymbolException {
        List symL = sym.getSymbols();
        int lb1 = symL.size() - 1;
        Symbol firstS;
        if(lb1 == 1) {
          firstS = (Symbol) symL.get(0);
        } else {
          firstS = firstA.getSymbol(symL.subList(0, lb1));
        }
        Distribution dist = getDistribution(firstS);
        return dtc.getCount(dist, (AtomicSymbol) symL.get(lb1));
      }
    });
  }
  
  private class UniformNullModel
  extends AbstractDistribution implements Serializable 
  {
    private static final long serialVersionUID = -3357793043843515032L;
    private Distribution nullModel = new UniformDistribution(
      (FiniteAlphabet) lastA
    );
    
    public Alphabet getAlphabet() {
      return AbstractOrderNDistribution.this.getAlphabet();
    }
    
    protected double getWeightImpl(AtomicSymbol sym)
    throws IllegalSymbolException {
      List symL = sym.getSymbols();
      int lb1 = symL.size() - 1;
      return nullModel.getWeight((AtomicSymbol) symL.get(lb1));
    }
    
    protected void setWeightImpl(AtomicSymbol sym, double weight)
    throws ChangeVetoException {
      throw new ChangeVetoException(
        "Can't change the weight of this null model"
      );
    }
    
    public Distribution getNullModel() {
      return this;
    }
    
    protected void setNullModelImpl(Distribution nullModel)
    throws IllegalAlphabetException, ChangeVetoException {
      throw new ChangeVetoException(
        "Can't set the null model for NthOrderDistribution.UniformNullModel"
      );
    }
  }
}
