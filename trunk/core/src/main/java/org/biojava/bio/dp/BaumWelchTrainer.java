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

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.onehead.SingleDP;
import org.biojava.bio.dp.onehead.SingleDPMatrix;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
* <p>
* Train a hidden markov model using maximum likelihood.
* </p>
*
* <p>
* Note: this class currently only works for one-head models.
* </p>
*
* @author Matthew Pocock
* @author Thomas Down
* @author Todd Riley
* @since 1.0
*/
public class BaumWelchTrainer extends AbstractTrainer implements Serializable {
 protected double singleSequenceIteration(
   ModelTrainer trainer,
   SymbolList symList
 ) throws IllegalSymbolException, IllegalTransitionException, IllegalAlphabetException {
   ScoreType scoreType = ScoreType.PROBABILITY;
   SingleDP dp = (SingleDP) getDP();
   State [] states = dp.getStates();
   int [][] backwardTransitions = dp.getBackwardTransitions();
   double [][] backwardTransitionScores = dp.getBackwardTransitionScores(scoreType);
   MarkovModel model = dp.getModel();

   SymbolList [] rll = { symList };

   // System.out.print("Forward...  ");
   SingleDPMatrix fm = (SingleDPMatrix) dp.forwardMatrix(rll, scoreType);
   double fs = fm.getScore();
   // System.out.println("Score = " + fs);

   // System.out.print("Backward... ");
   SingleDPMatrix bm = (SingleDPMatrix) dp.backwardMatrix(rll, scoreType);
   // System.out.println("Score = " + bs);

   Symbol gap = AlphabetManager.getGapSymbol();

   // state trainer
   for (int i = 1; i <= symList.length(); i++) {
     Symbol sym = symList.symbolAt(i);
     double [] fsc = fm.scores[i];
     double [] bsc = bm.scores[i];
     for (int s = 0; s < dp.getDotStatesIndex(); s++) {
       if (! (states[s] instanceof MagicalState)) {
         trainer.addCount(
           ((EmissionState) states[s]).getDistribution(),
           sym,
           mathExp(fsc[s] + bsc[s] - fs)
         );
       }
     }
   }

   // transition trainer
   for (int i = 0; i <= symList.length(); i++) {
     Symbol sym = (i < symList.length())
           ? symList.symbolAt(i + 1)
           : gap;
     double [] fsc = fm.scores[i];
     double [] bsc = bm.scores[i+1];
     double [] bsc2 = bm.scores[i];
     double[] weightVector = dp.getEmission(sym, scoreType);
     for (int s = 0; s < states.length; s++) {  // any -> emission transitions
       int [] ts = backwardTransitions[s];
       double [] tss = backwardTransitionScores[s];
       Distribution dist = model.getWeights(states[s]);
       for (int tc = 0; tc < ts.length; tc++) {
         int t = ts[tc];
         if(t < dp.getDotStatesIndex()) {
           double weight = mathExp(weightVector[t]);
           if (weight != 0.0) {
             trainer.addCount(
               dist, states[t],
               mathExp(
                 fsc[s] + tss[tc] + bsc[t]
                 -
                 fs
               ) * weight
             );
           }
         } else {
           trainer.addCount(
             dist, states[t],
             mathExp(
               fsc[s] + tss[tc] + bsc2[t]
               -
               fs
             )
           );
         }
       }
     }
   }

   return fs;
 }


   public double mathExp(double arg) {
                //Double argObj = new Double(arg);
                Double resultObj;

                if (Double.isNaN(arg)) {
                    //System.err.println("NaN encountered as arg to Math.exp in BaumWelch Training Loop");
                    arg = Double.NEGATIVE_INFINITY;
           //System.exit(-1);
                }
                resultObj = new Double(Math.exp(arg));
                return(resultObj.doubleValue());
   }



 public BaumWelchTrainer(DP dp) {
   super(dp);
 }
}


