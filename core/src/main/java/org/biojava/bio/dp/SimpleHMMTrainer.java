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

import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

public class SimpleHMMTrainer
    implements HMMTrainer
{
    DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
    FiniteAlphabet states;
    MarkovModel model;

    public SimpleHMMTrainer(MarkovModel model)
        throws IllegalSymbolException
    {
        this.model = model;

        // go thru model and add the Distributions
        states = model.stateAlphabet();

        for (Iterator stateI = states.iterator(); stateI.hasNext(); ) {
            State thisState = (State) stateI.next();
            // add emission Distributions
            if (thisState instanceof EmissionState) {
                EmissionState thisEmitter = (EmissionState) thisState;
                Distribution emissionDist = thisEmitter.getDistribution();
                dtc.registerDistribution(emissionDist);
                emissionDist.registerWithTrainer(dtc);
            }

            // add transition Distributions
            Distribution transDist = model.getWeights(thisState);
            dtc.registerDistribution(transDist);
            transDist.registerWithTrainer(dtc);
        }
    }

    public void startCycle()
    {
        dtc.clearCounts();
    }

    public void recordEmittedSymbol(State state, Symbol symbol, double weight)
        throws IllegalSymbolException
    {
        // look up the emission Distribution I need
        if (state instanceof EmissionState) {
            Distribution emissionDist = ((EmissionState) state).getDistribution();
            dtc.addCount(emissionDist, symbol, weight);
        }
        else throw new IllegalSymbolException("specified State is not an EmissionState.");
    }

    public void recordTransition(State source, State dest, double weight)
        throws IllegalArgumentException
    {
        // verify the transition
        try {
        if (model.containsTransition(source, dest)) {
            Distribution transDist = model.getWeights(source);
            dtc.addCount(transDist, dest, weight);
        }
        else throw new IllegalArgumentException("the specified transition is illegal for this model.");
        }
        catch (IllegalSymbolException ise) {
            throw new IllegalArgumentException("either source or destination are not valid");
        }
    }

    public void completeCycle()
        throws BioException
    {
        try {
            dtc.train();
        }
        catch (ChangeVetoException cve) {
            throw new BioException(cve);
        }
    }
}

