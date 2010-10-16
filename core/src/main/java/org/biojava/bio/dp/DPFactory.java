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


package org.biojava.bio.dp;

import org.biojava.bio.BioException;
import org.biojava.bio.dp.onehead.SingleDP;
import org.biojava.bio.dp.twohead.CellCalculatorFactoryMaker;
import org.biojava.bio.dp.twohead.DPInterpreter;
import org.biojava.bio.dp.twohead.PairwiseDP;

/**
 * The interface for objects that can generate a DP object for a MarkovModel.
 *
 * @author Matthew Pocock
 */ 
public interface DPFactory {
  public DP createDP(MarkovModel model)
  throws IllegalArgumentException, BioException;
  
  public static final DPFactory DEFAULT = new DefaultFactory(new DPInterpreter.Maker());
  
  public static class DefaultFactory implements DPFactory {
    private final CellCalculatorFactoryMaker cfFacM;
    
    public DefaultFactory(CellCalculatorFactoryMaker cfFacM) {
      this.cfFacM = cfFacM;
    }
    
    public DP createDP(MarkovModel model)
    throws IllegalArgumentException, BioException {
      int heads = model.advance().length;
      MarkovModel flat = DP.flatView(model);
      if(heads == 1) {
        return new SingleDP(flat);
      } else if(heads == 2) {
        return new PairwiseDP(flat, cfFacM);
      } else {
        throw new IllegalArgumentException(
          "Can't create DPFactory for models with " + heads + " heads"
        );
      }
    }
  }
}
