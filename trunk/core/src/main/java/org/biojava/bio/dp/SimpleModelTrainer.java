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

import java.io.Serializable;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class SimpleModelTrainer
extends SimpleDistributionTrainerContext
implements ModelTrainer, Serializable {
  private Set models = new HashSet();

  public void registerModel(MarkovModel model) {
    if(!models.contains(model)) {
      for(Iterator i = model.stateAlphabet().iterator(); i.hasNext(); ) {
        State s = (State) i.next();
        try {
          Distribution dist = model.getWeights(s);
          registerDistribution(dist);
        } catch (IllegalSymbolException ise) {
          throw new BioError("Couldn't register states from model", ise);
        }
        if(s instanceof EmissionState) {
          Distribution dist = ((EmissionState) s).getDistribution();
          registerDistribution(dist);
        }
      }
    }
  }
}
