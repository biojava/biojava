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
 * Implements a concurrency wrapper for a {@link PairwiseSequenceAligner}.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class CallablePairwiseSequenceAligner<S extends Sequence<C>, C extends Compound>
        implements Callable<SequencePair<S, C>> {

    private PairwiseSequenceAligner<S, C> psa;

    /**
     * Creates a pairwise sequence alignment task for simplified parallel execution.
     *
     * @param psa already initialized pairwise sequence aligner
     */
    public CallablePairwiseSequenceAligner(PairwiseSequenceAligner<S, C> psa) {
        this.psa = psa;
    }

    @Override
    public SequencePair<S, C> call() {
        return psa.getPair();
    }

}
