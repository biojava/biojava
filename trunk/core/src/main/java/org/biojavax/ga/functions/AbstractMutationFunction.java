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

package org.biojavax.ga.functions;

import org.biojava.bio.dist.OrderNDistribution;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * Abstract implementation of <code>MutationFunction</code> all custom
 * implementations should inherit from here.
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */
public abstract class AbstractMutationFunction
    extends AbstractChangeable implements MutationFunction {

  private double[] mutationProbs;
  private OrderNDistribution mutationSpectrum;

  protected AbstractMutationFunction() {
    mutationProbs = MutationFunction.DEFAULT_MUTATION_PROBS;
  }

  public final void setMutationProbs(double[] mutationProbs) throws ChangeVetoException {
    if(!hasListeners()){
      this.mutationProbs = mutationProbs;
    }else{
      ChangeEvent ce = new ChangeEvent(this,
                                       MutationFunction.MUTATION_PROBS,
                                       mutationProbs,
                                       this.mutationProbs
                                       );
      ChangeSupport changeSupport = super.getChangeSupport(MutationFunction.MUTATION_PROBS);
      synchronized(changeSupport){
        changeSupport.firePreChangeEvent(ce);
          this.mutationProbs = mutationProbs;
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }
  public final double[] getMutationProbs() {
    return mutationProbs;
  }

  public final void setMutationSpectrum(OrderNDistribution mutationSpectrum)
      throws ChangeVetoException {

    if(!hasListeners()){
      this.mutationSpectrum = mutationSpectrum;
    }else{
      ChangeEvent ce = new ChangeEvent(this,
                                       MutationFunction.MUTATION_SPECTRUM,
                                       mutationSpectrum,
                                       this.mutationSpectrum
                                       );
      ChangeSupport changeSupport =
          super.getChangeSupport(MutationFunction.MUTATION_SPECTRUM);
      synchronized(changeSupport){
        changeSupport.firePreChangeEvent(ce);
          this.mutationSpectrum = mutationSpectrum;
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }
  public final OrderNDistribution getMutationSpectrum() {
    return mutationSpectrum;
  }
}