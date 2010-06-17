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

import java.util.Arrays;

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
        implements PairwiseSequenceAligner<S, C> {

    // input fields
    private S query, target;
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean storingScoreMatrix;

    // output fields
    private short max, min;
    protected short score;
    protected short[][] scores;
    protected SequencePair<S, C> pair;
    protected long time = -1;

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
        this.query = query;
        this.target = target;
        this.gapPenalty = gapPenalty;
        this.subMatrix = subMatrix;
    }

    /**
     * Returns the query {@link Sequence}.
     *
     * @return the first {@link Sequence} of the pair to align
     */
    public S getQuery() {
        return query;
    }

    /**
     * Returns the target {@link Sequence}.
     *
     * @return the second {@link Sequence} of the pair to align
     */
    public S getTarget() {
        return target;
    }

    /**
     * Returns the gap penalties.
     *
     * @return the gap penalties used during alignment
     */
    public GapPenalty getGapPenalty() {
        return gapPenalty;
    }

    /**
     * Returns the substitution matrix.
     *
     * @return the set of substitution scores used during alignment
     */
    public SubstitutionMatrix<C> getSubstitutionMatrix() {
        return subMatrix;
    }

    /**
     * Returns choice to cache the score matrix or to save memory by deleting score matrix after alignment.
     *
     * @return choice to cache the score matrix
     */
    public boolean isStoringScoreMatrix() {
        return storingScoreMatrix;
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

    /**
     * Sets the gap penalties.
     *
     * @param gapPenalty the gap penalties used during alignment
     */
    public void setGapPenalty(GapPenalty gapPenalty) {
        this.gapPenalty = gapPenalty;
        reset();
    }

    /**
     * Sets the substitution matrix.
     *
     * @param subMatrix the set of substitution scores used during alignment
     */
    public void setSubstitutionMatrix(SubstitutionMatrix<C> subMatrix) {
        this.subMatrix = subMatrix;
        reset();
    }

    /**
     * Sets choice to cache the score matrix or to save memory by deleting score matrix after alignment.
     *
     * @param storingScoreMatrix choice to cache the score matrix
     */
    public void setStoringScoreMatrix(boolean storingScoreMatrix) {
        this.storingScoreMatrix = storingScoreMatrix;
        if (!storingScoreMatrix) {
            scores = null;
        }
    }

    // methods for MatrixAligner

    @Override
    public short[][] getScoreMatrix() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        short[][] copy = scores;
        if (tempStoringScoreMatrix) {
            copy = new short[scores.length][scores[0].length];
            for (int i = 0; i < copy.length; i++) {
                copy[i] = Arrays.copyOf(scores[i], scores[i].length);
            }
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return copy;
    }

    @Override
    public String getScoreMatrixAsString() {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return null;
            }
        }
        StringBuilder s = new StringBuilder();
        CompoundSet<C> compoundSet = query.getCompoundSet();
        int lengthCompound = compoundSet.getMaxSingleCompoundStringLength(), lengthRest =
                Math.max(Math.max(Short.toString(min).length(), Short.toString(max).length()), lengthCompound) + 1;
        String newLine = System.getProperty("line.separator"),
                padCompound = "%" + Integer.toString(lengthCompound) + "s",
                padRest = "%" + Integer.toString(lengthRest);
        s.append(String.format(padCompound, ""));
        s.append(String.format(padRest + "s", ""));
        for (C col : target.getAsList()) {
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(newLine);
        for (int row = 0; row <= query.getLength(); row++) {
            s.append(String.format(padCompound, (row == 0) ? "" :
                    compoundSet.getStringForCompound(query.getCompoundAt(row))));
            for (int col = 0; col <= target.getLength(); col++) {
                s.append(String.format(padRest + "d", getScoreMatrixAt(row, col)));
            }
            s.append(newLine);
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return s.toString();
    }

    @Override
    public short getScoreMatrixAt(int queryIndex, int targetIndex) {
        boolean tempStoringScoreMatrix = storingScoreMatrix;
        if (scores == null) {
            storingScoreMatrix = true;
            align();
            if (scores == null) {
                return Short.MIN_VALUE;
            }
        }
        short score = scores[queryIndex][targetIndex];
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return score;
    }

    // methods for Aligner

    @Override
    public long getComputationTime() {
        if (pair == null) {
            align();
        }
        return time;
    }

    @Override
    public Profile<S, C> getProfile() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // methods for Scorer

    @Override
    public int getMaxScore() {
        if (pair == null) {
            align();
        }
        return max;
    }

    @Override
    public int getMinScore() {
        if (pair == null) {
            align();
        }
        return min;
    }

    @Override
    public int getScore() {
        if (pair == null) {
            align();
        }
        return score;
    }

    // method for PairwiseSequenceScorer

    @Override
    public SequencePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // helper method that performs alignment
    protected abstract void align();

    // helper method that resets output fields
    protected void reset() {
        if (query != null && target != null && gapPenalty != null && subMatrix != null) {
            int subLength = Math.min(query.getLength(), target.getLength()), maxLength = query.getLength()
                    + target.getLength(), penalties = gapPenalty.getOpenPenalty() + gapPenalty.getExtensionPenalty();
            max = (short) (subLength * subMatrix.getMaxValue());
            score = min = (short) Math.min(subLength * subMatrix.getMinValue() + (maxLength - subLength) * penalties,
                    maxLength * penalties);
        }
        scores = null;
        pair = null;
        time = -1;
    }

}
