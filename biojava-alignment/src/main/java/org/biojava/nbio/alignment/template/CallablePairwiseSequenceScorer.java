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
 * Created on June 17, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment.template;

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.concurrent.Callable;

/**
 * Implements a concurrency wrapper for a {@link PairwiseSequenceScorer}.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the pair is of type S
 * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
 */
public class CallablePairwiseSequenceScorer<S extends Sequence<C>, C extends Compound> implements Callable<Double> {

    private PairwiseSequenceScorer<S, C> pss;

    /**
     * Creates a pairwise sequence scoring task for simplified parallel execution.
     *
     * @param pss already initialized pairwise sequence scorer
     */
    public CallablePairwiseSequenceScorer(PairwiseSequenceScorer<S, C> pss) {
        this.pss = pss;
    }

    @Override
    public Double call() {
        return pss.getScore();
    }

}
