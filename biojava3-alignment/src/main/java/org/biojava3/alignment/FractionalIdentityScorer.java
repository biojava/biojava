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
 * Created on June 21, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.PairwiseSequenceScorer;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements an algorithm which computes a score for a sequence alignment pair.  The reported score is the number of
 * alignment columns which have identical {@link Compound}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class FractionalIdentityScorer<S extends Sequence<C>, C extends Compound>
        implements PairwiseSequenceScorer<S, C> {

    private PairwiseSequenceAligner<S, C> aligner;
    private SequencePair<S, C> pair;
    private int max, score;

    /**
     * Creates a fractional identity scorer for a pair of sequences aligned by the given pairwise sequence aligner.
     *
     * @param aligner a pairwise sequence aligner
     */
    public FractionalIdentityScorer(PairwiseSequenceAligner<S, C> aligner) {
        this.aligner = aligner;
    }

    /**
     * Creates a fractional identity scorer for an aligned pair of sequences.
     *
     * @param pair an aligned pair of sequences
     */
    public FractionalIdentityScorer(SequencePair<S, C> pair) {
        set(pair);
    }

    @Override
    public SequencePair<S, C> getPair() {
        if (pair == null && aligner != null) {
            set(aligner.getPair());
        }
        return pair;
    }

    @Override
    public int getMaxScore() {
        if (pair == null && aligner != null) {
            set(aligner.getPair());
        }
        return max;
    }

    @Override
    public int getMinScore() {
        return 0;
    }

    @Override
    public int getScore() {
        if (pair == null && aligner != null) {
            set(aligner.getPair());
        }
        return score;
    }

    // helper method for initialization
    private void set(SequencePair<S, C> pair) {
        this.pair = pair;
        max = pair.getLength();
        score = pair.getNumIdenticals();
    }

}
