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

import org.biojava.bio.BioError;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * Abstract implementation of <code>CrossOverFunction</code>. All custom
 * implementations should inherit from here.
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */
public abstract class AbstractCrossOverFunction extends AbstractChangeable 
        implements CrossOverFunction {
  private int maxCross;
  private double[] crossProbs;

  protected AbstractCrossOverFunction() {
    try {
      setMaxCrossOvers(CrossOverFunction.DEFAULT_MAX_CROSS);
      setCrossOverProbs(CrossOverFunction.DEFAULT_CROSS_PROB);
    }
    catch (ChangeVetoException ex) {
      throw new BioError("Cannot set the default values of the CrossOverFunction", ex);
    }
  }

  public final void setMaxCrossOvers(int maxCrossOvers) throws ChangeVetoException {
    if(!hasListeners()){
      maxCross = maxCrossOvers;
    }else{
      ChangeEvent ce = new ChangeEvent(this,
                                       CrossOverFunction.MAX_CROSSES,
                                       new Integer(maxCrossOvers),
                                       new Integer(this.maxCross)
                                       );
      ChangeSupport changeSupport = super.getChangeSupport(CrossOverFunction.MAX_CROSSES);
      synchronized(changeSupport){
        changeSupport.firePreChangeEvent(ce);
        maxCross = maxCrossOvers;
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }

  public final int getMaxCrossOvers() {
    return maxCross;
  }

  public final void setCrossOverProbs(double[] crossOverProbs) throws ChangeVetoException {
    if(!hasListeners()){
      crossProbs = crossOverProbs;
    }else{
      ChangeEvent ce = new ChangeEvent(this,
                                       CrossOverFunction.CROSS_PROB,
                                       crossOverProbs,
                                       this.crossProbs
                                       );
      ChangeSupport changeSupport = super.getChangeSupport(CrossOverFunction.CROSS_PROB);
      synchronized(changeSupport){
        changeSupport.firePreChangeEvent(ce);
        crossProbs = crossOverProbs;
        changeSupport.firePostChangeEvent(ce);
      }
    }
  }
  public final double[] getCrossOverProbs() {
    return crossProbs;
  }
}