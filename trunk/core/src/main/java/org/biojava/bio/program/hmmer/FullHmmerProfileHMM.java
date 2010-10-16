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

//package org.biojava.bio.program.hmmer;
package org.biojava.bio.program.hmmer;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.ModelInState;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.SimpleModelInState;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeVetoException;

/** 
 * This is a class for representing the full HMMER generated Profile HMM (including loop
 * states N and C terminal looping states).
 *
 * @author Lachlan Coin
 * @since 1.3
 */
public class FullHmmerProfileHMM extends SimpleMarkovModel{

    EmissionState j;
    EmissionState c;
    EmissionState n;
    ModelInState hmmState;

  private final static int [] advance = { 1 };

 FullHmmerProfileHMM(
    HmmerProfileHMM hmm
    ) throws IllegalSymbolException, IllegalTransitionException,
  IllegalAlphabetException, ChangeVetoException {
     super(1,hmm.emissionAlphabet(),hmm.stateAlphabet().getName());

     Distribution nullDist = hmm.getInsert(1).getDistribution();
     hmmState = new SimpleModelInState(hmm, hmm.stateAlphabet().getName());

     j = new SimpleEmissionState(
        "j",
        Annotation.EMPTY_ANNOTATION,
        advance,
	nullDist
      );


     c = new SimpleEmissionState(
        "j",
        Annotation.EMPTY_ANNOTATION,
        advance,
	nullDist
      );


     n = new SimpleEmissionState(
        "j",
        Annotation.EMPTY_ANNOTATION,
        advance,
        nullDist
      );
     addState(j);
     addState(c);
     addState(n);
     addState(hmmState);




     createTransition(magicalState(), n);
     createTransition(n, hmmState);
     createTransition(n, n);
     createTransition(hmmState,c);
     createTransition(hmmState,j);
     createTransition(j,hmmState);
     createTransition(j,j);
     createTransition(c,c);
     createTransition(c,magicalState());

 }

    /** Gets the J loop state */
    public EmissionState jState(){
	return j;
    }

    /** Gets the c loop state */
    public EmissionState cState(){
	return c;
    }

    /** Gets the n loop state */
    public EmissionState nState(){
	return n;
    }

    /** Gets the inner HmmerProfileHMM state */
    public ModelInState hmm(){ 
	return hmmState;
    }

}
