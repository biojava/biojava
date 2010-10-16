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

import org.biojava.bio.BioError;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ReversibleTranslationTable;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A translated view of some underlying distribution.  The <code>getWeight</code>
 * method returns the result of calling <code>getWeight</code> on the underlying
 * distribution, having first translated the <code>Symbol</code> parameter using
 * the supplied <code>ReversibleTranslationTable</code>.  All changes to the
 * underlying distribution are reflected by the <code>TranslatedDistribution</code>.
 *
 * <p>
 * The <code>TranslatedDistribution</code> is not directly mutable: calling
 * <code>setWeight</code> will result in a <code>ChangeVetoException</code>.
 * However, a <code>DistributionTrainer</code> may be registered for a
 * <code>TranslatedDistribution</code>.  Any counts received by this trainer
 * are untranslated then forwarded to the underlying distribution.  It is
 * valid to add counts to both a <code>TranslatedDistribution</code> and
 * its underlying distribution in a single training session, so
 * <code>TranslatedDistribution</code> objects are useful for tying
 * parameters together when training Markov Models.
 * </p>
 *
 * <h2>Example usage</h2>
 *
 * <pre>
 * Distribution d = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
 * d.setWeight(DNATools.a(), 0.7);
 * d.setWeight(DNATools.c(), 0.1);
 * d.setWeight(DNATools.g(), 0.1);
 * d.setWeight(DNATools.t(), 0.1);
 * Distribution complemented = new TranslatedDistribution(
 *     DNATools.complementTable(),
 *     d,
 *     DistributionFactory.DEFAULT
 * );
 * System.out.println(
 *    "complemented.getWeight(DNATools.t()) = " +
 *    complemented.getWeight(DNATools.t())
 * );  // Should print 0.7
 * </pre>
 * 
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @since 1.1
 */
public class TranslatedDistribution
  extends
    AbstractChangeable
  implements
    Distribution,
    Serializable
{
  private final Distribution other;
  private final Distribution delegate;
  private final ReversibleTranslationTable table;
  private transient ChangeListener forwarder;
  private transient ChangeListener delegateUpdate;

  /**
   * Create a new TranslatedDistribution. Make these things via getDistribuiton.
   *
   * @param table    a ReversibleTranslationTable used to map the symbols
   * @param other    the underlying ditribution
   * @param distFact a DistributionFactory used to create a delegate for
   *    stooring mapped weights
   */
  public TranslatedDistribution(
    ReversibleTranslationTable table,
    Distribution other,
    DistributionFactory distFact
  ) throws IllegalAlphabetException {
    if (! (other.getAlphabet() instanceof FiniteAlphabet)) {
        throw new IllegalAlphabetException("The current implementation of TranslatedDistribution is only valid for distributions over finite alphabets");
    }
      
    if(!table.getTargetAlphabet().equals(other.getAlphabet())) {
      throw new IllegalAlphabetException(
        "Table target alphabet and distribution alphabet don't match: " +
        table.getTargetAlphabet().getName() + " and " +
        other.getAlphabet().getName() + " without symbol "
      );
    }
    this.other = other;
    this.table = table;
    this.delegate = distFact.createDistribution(table.getSourceAlphabet());
    
    syncDelegate();
    
    delegateUpdate = new ChangeListener() {
        public void preChange(ChangeEvent ce) {}
        public void postChange(ChangeEvent ce) {
            ChangeType ct = ce.getType();
            Object change = ce.getChange();
            if(ct == Distribution.WEIGHTS) {
                boolean synced = false;
                if((change != null) && (change instanceof Object[]) ) {
                    Object[] ca = (Object[]) change;
                    if( (ca.length == 2) && (ca[0] instanceof Symbol) && (ca[1] instanceof Number)) {
                        try {
                            delegate.setWeight(
                                (Symbol) ca[0],
                                ((Number) ca[1]).doubleValue()
                            );
                            synced = true;
                        } catch (Exception ise) {
                            throw new BioError("Couldn't synchronize weight", ise);
                        }
                    }
                }
                if (!synced) {
                    // Weights have changed, but we can't understand the event, so re-sync them
                    // all.
                    syncDelegate();
                }
            }
        }
    } ;
    addChangeListener(delegateUpdate);
  }
  
  private void syncDelegate() {
      for (Iterator i = ((FiniteAlphabet) delegate.getAlphabet()).iterator(); i.hasNext(); ) {
        Symbol s = (Symbol) i.next();
        try {
            delegate.setWeight(s, other.getWeight(table.untranslate(s)));
        } catch (Exception ex) {
            throw new BioError(ex, "Assertion failed: couldn't map distributions");
        }
    }
  }

  public Alphabet getAlphabet() {
    return table.getSourceAlphabet();
  }

  public double getWeight(Symbol sym)
    throws IllegalSymbolException
  {
    return delegate.getWeight(sym);
  }

  public void setWeight(Symbol sym, double weight)
    throws IllegalSymbolException, ChangeVetoException 
  {
    throw new ChangeVetoException("Can't directly edit a TranslatedDistribution");
  }

  public Symbol sampleSymbol() {
    return delegate.sampleSymbol();
  }

  public Distribution getNullModel() {
    return delegate.getNullModel();
  }

  public void setNullModel(Distribution dist)
  throws IllegalAlphabetException, ChangeVetoException {
    delegate.setNullModel(dist);
  }

  /**
   * Retrieve the translation table encapsulating the map from this emission
   * spectrum to the underlying one.
   *
   * @return a ReversibleTranslationtTable
   */
  public ReversibleTranslationTable getTable() {
    return table;
  }

  public void registerWithTrainer(DistributionTrainerContext dtc) {
    dtc.registerDistribution(other);

    dtc.registerTrainer(this, new DistributionTrainer() {
      public void addCount(
        DistributionTrainerContext dtc,
        AtomicSymbol s,
        double count
      ) throws IllegalSymbolException {
        dtc.addCount(other, table.translate(s), count);
      }

      public double getCount(
        DistributionTrainerContext dtc,
        AtomicSymbol s
      ) throws IllegalSymbolException {
        return dtc.getCount(other, table.translate(s));
      }

      public void train(DistributionTrainerContext dtc, double weight)
        throws ChangeVetoException 
      {
          // This is a no-op, since our counts have already been passed on to
          // the sister Distribution.
      }

      public void clearCounts(DistributionTrainerContext dtc) {
      }
    });
  }

  protected ChangeSupport getChangeSupport(ChangeType ct) {
    ChangeSupport cs = super.getChangeSupport(ct);

    if(forwarder == null &&
       (Distribution.WEIGHTS.isMatchingType(ct) || ct.isMatchingType(Distribution.WEIGHTS)))
    {
      forwarder = new Forwarder(this, cs);
      other.addChangeListener(forwarder, Distribution.WEIGHTS);
    }

    return cs;
  }

  private class Forwarder extends ChangeForwarder {
    public Forwarder(Object source, ChangeSupport changeSupport) {
      super(source, changeSupport);
    }

    protected ChangeEvent generateChangeEvent(ChangeEvent ce) {
      ChangeType ct = ce.getType();
      Object change = ce.getChange();
      Object previous = ce.getPrevious();
      if(ct == Distribution.WEIGHTS) {
        if( (change != null) && (change instanceof Object[]) ) {
          Object[] ca = (Object[]) change;
          if( (ca.length == 2) && (ca[0] instanceof Symbol) ) {
            try {
              change = new Object[] { table.translate((Symbol) ca[0]), ca[1] };
            } catch (IllegalSymbolException ise) {
              throw new BioError("Couldn't translate symbol", ise);
            }
          }
        }
        if( (previous != null) && (previous instanceof Object[]) ) {
          Object[] pa = (Object[]) previous;
          if( (pa.length == 2) && (pa[0] instanceof Symbol) ) {
            try {
              previous = new Object[] { table.translate((Symbol) pa[0]), pa[1] };
            } catch (IllegalSymbolException ise) {
              throw new BioError("Couldn't translate symbol", ise);
            }
          }
        }
      } else if(ct == Distribution.NULL_MODEL) {
        change = null;
        previous = null;
      }
      return new ChangeEvent(
        TranslatedDistribution.this, ct,
        change, previous, ce
      );
    }
  }
}
