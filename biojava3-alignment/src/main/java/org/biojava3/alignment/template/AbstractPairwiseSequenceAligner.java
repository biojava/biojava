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

package org.biojava3.alignment.template;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

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
        super(gapPenalty, subMatrix);
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

    // methods for PairwiseSequenceScorer

    @Override
    public S getQuery() {
        return query;
    }

    @Override
    public S getTarget() {
        return target;
    }

    // method for MatrixAligner

    @Override
    public String getScoreMatrixAsString() {
        boolean tempStoringScoreMatrix = isStoringScoreMatrix();
        if (scores == null) {
            setStoringScoreMatrix(true);
            align();
            if (scores == null) {
                return null;
            }
        }
        StringBuilder s = new StringBuilder();
        CompoundSet<C> compoundSet = query.getCompoundSet();
        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
                Math.max(Math.max(Short.toString(min).length(), Short.toString(max).length()), lengthCompound) + 1;
        String padCompound = "%" + Integer.toString(lengthCompound) + "s",
                padRest = "%" + Integer.toString(lengthRest);
        s.append(String.format(padCompound, ""));
        s.append(String.format(padRest + "s", ""));
        for (C col : target.getAsList()) {
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(String.format("%n"));
        for (int row = 0; row <= query.getLength(); row++) {
            s.append(String.format(padCompound, (row == 0) ? "" :
                    compoundSet.getStringForCompound(query.getCompoundAt(row))));
            for (int col = 0; col <= target.getLength(); col++) {
                s.append(String.format(padRest + "d", getScoreMatrixAt(row, col)));
            }
            s.append(String.format("%n"));
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return s.toString();
    }

    // method for PairwiseSequenceScorer

    @Override
    public SequencePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // helper methods

    // prepares for alignment; returns true if everything is set to run the alignment
    @Override
    protected boolean alignReady() {
        reset();
        if (query != null && target != null && getGapPenalty() != null && getSubstitutionMatrix() != null &&
                query.getCompoundSet().equals(target.getCompoundSet())) {
            scores = new short[query.getLength() + 1][target.getLength() + 1];
            return true;
        }
        return false;
    }

    // scores alignment of two columns
    @Override
    protected short alignScoreColumns(int queryColumn, int targetColumn) {
        return getSubstitutionMatrix().getValue(query.getCompoundAt(queryColumn), target.getCompoundAt(targetColumn));
    }

    // resets output fields
    @Override
    protected void reset() {
        if (query != null && target != null && getGapPenalty() != null && getSubstitutionMatrix() != null &&
                query.getCompoundSet().equals(target.getCompoundSet())) {
            int maxq = 0, maxt = 0;
            for (C c : query) {
                maxq += getSubstitutionMatrix().getValue(c, c);
            }
            for (C c : target) {
                maxt += getSubstitutionMatrix().getValue(c, c);
            }
            max = (short) Math.max(maxq, maxt);
            score = min = (short) (2 * getGapPenalty().getOpenPenalty() + (query.getLength() + target.getLength()) *
                    getGapPenalty().getExtensionPenalty());
        }
        scores = null;
        pair = null;
        time = -1;
    }

}
