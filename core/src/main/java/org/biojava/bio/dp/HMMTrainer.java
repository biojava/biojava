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
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * interface implemented by objects that train HMMs.
 *
 * @author David Huen
 */
public interface HMMTrainer
{
    /**
     * called to put the trainer into an initial state for a new
     * round of training.
     */
    public void startCycle();

    /**
     * record that the specified symbol was emitted from the specified state.
     */
    public void recordEmittedSymbol(State state, Symbol symbol, double weight)
        throws IllegalSymbolException;

    /**
     * record that a transition was observed between the specified states.
     */
    public void recordTransition(State source, State dest, double weight)
        throws IllegalArgumentException;

    /**
     * indicate that a cycle of training is completed and the
     * emission/transition matrices should be updated.
     */
    public void completeCycle()
        throws BioException;
}

