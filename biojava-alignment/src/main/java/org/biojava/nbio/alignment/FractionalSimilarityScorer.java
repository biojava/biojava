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

package org.biojava.nbio.alignment;

import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.*;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Implements an algorithm which computes a score for a sequence alignment pair.  The reported score is the number of
 * alignment columns which have similar {@link Compound}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class FractionalSimilarityScorer<S extends Sequence<C>, C extends Compound> extends AbstractScorer
        implements PairwiseSequenceScorer<S, C> {

    // always stored
    private S query, target;
    private int max, score;

    // optional cached input field
    private PairwiseSequenceAligner<S, C> aligner;

    /**
     * Creates a fractional similarity scorer for a pair of sequences aligned by the given pairwise sequence aligner.
     *
     * @param aligner a pairwise sequence aligner
     */
    public FractionalSimilarityScorer(PairwiseSequenceAligner<S, C> aligner) {
        query = aligner.getQuery();
        target = aligner.getTarget();
        this.aligner = aligner;
    }

    /**
     * Creates a fractional similarity scorer for an aligned pair of sequences.
     *
     * @param pair an aligned pair of sequences
     */
    public FractionalSimilarityScorer(SequencePair<S, C> pair) {
        query = pair.getQuery().getOriginalSequence();
        target = pair.getTarget().getOriginalSequence();
        max = pair.getLength();
        score = pair.getNumSimilars();
    }

    // methods for PairwiseSequenceScorer

    @Override
    public S getQuery() {
        return query;
    }

    @Override
    public S getTarget() {
        return target;
    }

    // methods for Scorer

    @Override
    public double getMaxScore() {
        if (aligner != null) {
            align();
        }
        return max;
    }

    @Override
    public double getMinScore() {
        return 0;
    }

    @Override
    public double getScore() {
        if (aligner != null) {
            align();
        }
        return score;
    }

    // helper method for initialization from an aligner
    private void align() {
        max = aligner.getPair().getLength();
        score = aligner.getPair().getNumSimilars();
        aligner = null;
    }

}
