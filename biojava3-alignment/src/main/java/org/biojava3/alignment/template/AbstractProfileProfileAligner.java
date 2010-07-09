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

import java.util.List;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements common code for an {@link Aligner} for a pair of {@link Profile}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the pair of alignment {@link Profile}s is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractProfileProfileAligner<S extends Sequence<C>, C extends Compound>
        extends AbstractMatrixAligner<S, C> implements ProfileProfileAligner<S, C> {

    // additional input fields
    private Profile<S, C> query, target;

    // additional output field
    protected ProfilePair<S, C> pair;

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
        super(gapPenalty, subMatrix);
        this.query = query;
        this.target = target;
        reset();
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
        for (C col : target.getAlignedSequence(1).getAsList()) { // TODO print consensus sequences
            s.append(String.format(padRest + "s", compoundSet.getStringForCompound(col)));
        }
        s.append(String.format("%n"));
        for (int row = 0; row <= query.getLength(); row++) {
            s.append(String.format(padCompound, (row == 0) ? "" :
                    compoundSet.getStringForCompound(query.getAlignedSequence(1).getCompoundAt(row))));
            for (int col = 0; col <= target.getLength(); col++) {
                s.append(String.format(padRest + "d", getScoreMatrixAt(row, col)));
            }
            s.append(String.format("%n"));
        }
        setStoringScoreMatrix(tempStoringScoreMatrix);
        return s.toString();
    }

    // method for ProfileProfileScorer

    @Override
    public ProfilePair<S, C> getPair() {
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

    // scores alignment of two columns; TODO add caching of column compounds
    @Override
    protected short alignScoreColumns(int queryColumn, int targetColumn) {
        List<C> cslist = query.getCompoundSet().getAllCompounds();
        float[] qfrac = query.getCompoundWeightsAt(queryColumn, cslist),
                tfrac = target.getCompoundWeightsAt(targetColumn, cslist);
        float score = 0.0f;
        for (int q = 0; q < qfrac.length; q++) {
            if (qfrac[q] > 0.0f) {
                for (int t = 0; t < tfrac.length; t++) {
                    if (tfrac[t] > 0.0f) {
                        score += qfrac[q]*tfrac[t]*getSubstitutionMatrix().getValue(cslist.get(q), cslist.get(t));
                    }
                }
            }
        }
        return (short) Math.round(score);
    }

    // resets output fields; TODO better bounds for max and min
    @Override
    protected void reset() {
        GapPenalty gapPenalty = getGapPenalty();
        SubstitutionMatrix<C> subMatrix = getSubstitutionMatrix();
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
