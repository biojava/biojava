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
 * Created on June 24, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment.template;

import java.util.Arrays;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

public abstract class AbstractProfileProfileAligner<S extends Sequence<C>, C extends Compound>
        implements ProfileProfileAligner<S, C> {

    // input fields
    private Profile<S, C> query, target;
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private boolean storingScoreMatrix;

    // output fields
    private short max, min;
    protected short score;
    protected short[][] scores;
    protected ProfilePair<S, C> pair;
    protected long time = -1;

    /**
     * Before running a profile-profile alignment, data must be sent in via calls to
     * {@link #setQuery(Profile)}, {@link #setTarget(Profile)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    protected AbstractProfileProfileAligner() {
    }

    /**
     * Prepares for a profile-profile alignment.
     *
     * @param query the first {@link Profile} of the pair to align
     * @param target the second {@link Profile} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractProfileProfileAligner(Profile<S, C> query, Profile<S, C> target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        this.query = query;
        this.target = target;
        this.gapPenalty = gapPenalty;
        this.subMatrix = subMatrix;
    }

    /**
     * Returns the query {@link Profile}.
     *
     * @return the first {@link Profile} of the pair to align
     */
    public Profile<S, C> getQuery() {
        return query;
    }

    /**
     * Returns the target {@link Profile}.
     *
     * @return the second {@link Profile} of the pair to align
     */
    public Profile<S, C> getTarget() {
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
     * Sets the query {@link Profile}.
     *
     * @param query the first {@link Profile} of the pair to align
     */
    public void setQuery(Profile<S, C> query) {
        this.query = query;
        reset();
    }

    /**
     * Sets the target {@link Profile}.
     *
     * @param target the second {@link Profile} of the pair to align
     */
    public void setTarget(Profile<S, C> target) {
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
        for (C col : target.getAlignedSequence(0).getAsList()) { // TODO print consensus sequences
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(newLine);
        for (int row = 0; row <= query.getLength(); row++) {
            s.append(String.format(padCompound, (row == 0) ? "" :
                    compoundSet.getStringForCompound(query.getAlignedSequence(0).getCompoundAt(row))));
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

    // method for ProfileProfileScorer

    @Override
    public ProfilePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // helper method that performs alignment
    protected abstract void align();

    // helper method that resets output fields; TODO better bounds for max and min
    protected void reset() {
        if (query != null && target != null && gapPenalty != null && subMatrix != null) {
            int subLength = Math.min(query.getLength(), target.getLength()), maxLength = query.getLength()
                    + target.getLength(), penalties = gapPenalty.getOpenPenalty() + gapPenalty.getExtensionPenalty(),
                    sumSize = query.getSize() * target.getSize();
            max = (short) (subLength * subMatrix.getMaxValue() * sumSize);
            score = min = (short) (Math.min(subLength * subMatrix.getMinValue() + (maxLength - subLength) * penalties,
                    maxLength * penalties) * sumSize);
        }
        scores = null;
        pair = null;
        time = -1;
    }

}
