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
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements common code for an {@link Aligner} for a pair of {@link Sequence}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractPairwiseSequenceAligner<S extends Sequence<C>, C extends Compound>
        extends AbstractMatrixAligner<S, C> implements PairwiseSequenceAligner<S, C> {

    // additional input fields
    private S query, target;

    // additional output field
    protected SequencePair<S, C> pair;

    /**
     * Before running a pairwise global sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    protected AbstractPairwiseSequenceAligner() {
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        this(query, target, gapPenalty, subMatrix, false);
    }

    /**
     * Prepares for a pairwise sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     * @param local if true, find a region of similarity rather than aligning every compound
     */
    protected AbstractPairwiseSequenceAligner(S query, S target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix, boolean local) {
        super(gapPenalty, subMatrix, local);
        this.query = query;
        this.target = target;
        reset();
    }

    /**
     * Sets the query {@link Sequence}.
     *
     * @param query the first {@link Sequence} of the pair to align
     */
    public void setQuery(S query) {
        this.query = query;
        reset();
    }

    /**
     * Sets the target {@link Sequence}.
     *
     * @param target the second {@link Sequence} of the pair to align
     */
    public void setTarget(S target) {
        this.target = target;
        reset();
    }

    // method for PairwiseSequenceAligner

    @Override
    public SequencePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
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

    // methods for AbstractMatrixAligner

    @Override
    protected CompoundSet<C> getCompoundSet() {
        return (query == null) ? null : query.getCompoundSet();
    }

    @Override
    protected List<C> getCompoundsOfQuery() {
        return (query == null) ? new ArrayList<C>() : query.getAsList();
    }

    @Override
    protected List<C> getCompoundsOfTarget() {
        return (target == null) ? new ArrayList<C>() : target.getAsList();
    }

    @Override
    protected int[] getScoreMatrixDimensions() {
        return new int[] { (query == null) ? 1 : query.getLength() + 1, (target == null) ? 1 : target.getLength() + 1,
                (getGapPenalty() == null || getGapPenalty().getType() == GapPenalty.Type.LINEAR) ? 1 : 3 };
    }

    @Override
    protected int getSubstitutionScore(int queryColumn, int targetColumn) {
        return getSubstitutionMatrix().getValue(query.getCompoundAt(queryColumn), target.getCompoundAt(targetColumn));
    }

    @Override
    protected boolean isReady() {
        return query != null && target != null && getGapPenalty() != null && getSubstitutionMatrix() != null &&
                query.getCompoundSet().equals(target.getCompoundSet());
    }

    @Override
    protected void reset() {
        super.reset();
        pair = null;
        if (query != null && target != null && getGapPenalty() != null && getSubstitutionMatrix() != null &&
                query.getCompoundSet().equals(target.getCompoundSet())) {
            int maxq = 0, maxt = 0;
            for (C c : query) {
                maxq += getSubstitutionMatrix().getValue(c, c);
            }
            for (C c : target) {
                maxt += getSubstitutionMatrix().getValue(c, c);
            }
            max = (int) Math.max(maxq, maxt);
            score = min = isLocal() ? 0 : (int) (2 * getGapPenalty().getOpenPenalty() + (query.getLength() +
                    target.getLength()) * getGapPenalty().getExtensionPenalty());
        }
    }

}
