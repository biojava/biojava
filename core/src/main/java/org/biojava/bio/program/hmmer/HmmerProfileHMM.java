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

package org.biojava.bio.program.hmmer;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DotState;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.ProfileHMM;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
 




/** This is a class for representing HMMER generated Profile HMM.
 *  It differs from the normal ProfileHMM only in the states which are connected:
 *   - there are no insert <-> delete transitions allowed
 *   - there is no iO initial insert state (between begin and initial match states)
 *   - there is not iN final insert state (between final match state and end state)
 *
 *  @author Lachlan Coin
 */
public class HmmerProfileHMM extends ProfileHMM {
  protected HmmerProfileHMM(
    Alphabet alpha, int columns,
    DistributionFactory matchFactory, DistributionFactory insertFactory, String name
  ) throws IllegalSymbolException, IllegalTransitionException,
  IllegalAlphabetException {
      super(alpha,columns, matchFactory,insertFactory,name);
  }

 

    /** This is called by constructor in setting up the allowed transitions in the model */
  protected void connectModel() throws 
	ChangeVetoException, IllegalSymbolException, IllegalTransitionException,IllegalAlphabetException{
	EmissionState mO = getMatch(0);
	removeState(getInsert(0));
	removeState(getInsert(columns()));
	DotState dO = null;
	EmissionState iO =null;
	for(int i = 1; i <= columns(); i++){
	    EmissionState mN = getMatch(i);
	    EmissionState iN = getInsert(i);
	    DotState dN = getDelete(i);

	    // from a model state
	    createTransition(mO, mN);
	    if(i < columns())
		createTransition(mN, iN);
	    createTransition(mO, dN);
      
	    // from an insert state
	    if(i<columns())
		createTransition(iN, iN);
	   
	    
	    // from a delete state
	    if(i > 1) {
		createTransition(dO, dN);
		createTransition(dO, mN);
		createTransition(iO, mN);
	    }        
	    mO = mN;
	    iO = iN;
	    dO = dN;
	    // transitions to and from magical states and all match states
	    if(i>1)
		createTransition(magicalState(),mN);
	    createTransition(mN, magicalState());
	}
	// for the transitions to end
	createTransition(dO, magicalState());
    }

    public double transScore(State from, State to, Symbol symFrom, Symbol symTo) throws IllegalSymbolException{
	    return log2(getWeights(from).getWeight(to));
    }
    
    protected static double log2(double x){
	return Math.log(x)/Math.log(2);
    }
    
    protected EmissionState makeNewInsertState(String str, Annotation ann, int[] adv, Distribution dis){
	  return new ProfileEmissionState(str, ann,adv, dis);
  }
  
   protected EmissionState makeNewMatchState(String str, Annotation ann, int[] adv, Distribution dis){
	  return new ProfileEmissionState(str, ann,adv, dis);
  }
  
}
