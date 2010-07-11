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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

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

    // concurrent execution fields
    private Future<ProfilePair<S, C>> queryFuture, targetFuture;

    // cached fields
    private List<C> cslist;
    private float[][] qfrac, tfrac;

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
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align, still to be calculated
     * @param target the second {@link Profile} of the pair to align, still to be calculated
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractProfileProfileAligner(Future<ProfilePair<S, C>> query, Future<ProfilePair<S, C>> target,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(gapPenalty, subMatrix);
        queryFuture = query;
        targetFuture = target;
        reset();
    }

    /**
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align
     * @param target the second {@link Profile} of the pair to align, still to be calculated
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractProfileProfileAligner(Profile<S, C> query, Future<ProfilePair<S, C>> target,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(gapPenalty, subMatrix);
        this.query = query;
        targetFuture = target;
        reset();
    }

    /**
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align, still to be calculated
     * @param target the second {@link Profile} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    protected AbstractProfileProfileAligner(Future<ProfilePair<S, C>> query, Profile<S, C> target,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(gapPenalty, subMatrix);
        queryFuture = query;
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
        queryFuture = null;
        reset();
    }

    /**
     * Sets the target {@link Profile}.
     *
     * @param target the second {@link Profile} of the pair to align
     */
    public void setTarget(Profile<S, C> target) {
        this.target = target;
        targetFuture = null;
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
        // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
        try {
            if (query == null && queryFuture != null) {
                query = queryFuture.get();
            }
            if (target == null && targetFuture != null) {
                target = targetFuture.get();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }
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
        return alignScoreVectors(qfrac[queryColumn - 1], tfrac[targetColumn - 1]);
    }

    // scores alignment of two column vectors
    private short alignScoreVectors(float[] qv, float[] tv) {
        float score = 0.0f;
        for (int q = 0; q < qv.length; q++) {
            if (qv[q] > 0.0f) {
                for (int t = 0; t < tv.length; t++) {
                    if (tv[t] > 0.0f) {
                        score += qv[q]*tv[t]*getSubstitutionMatrix().getValue(cslist.get(q), cslist.get(t));
                    }
                }
            }
        }
        return (short) Math.round(score);
    }

    // resets output fields; caches profile vectors
    @Override
    protected void reset() {
        if (query != null && target != null && getGapPenalty() != null && getSubstitutionMatrix() != null &&
                query.getCompoundSet().equals(target.getCompoundSet())) {
            int maxq = 0, maxt = 0;
            cslist = query.getCompoundSet().getAllCompounds();
            qfrac = new float[query.getLength()][];
            for (int i = 0; i < qfrac.length; i++) {
                qfrac[i] = query.getCompoundWeightsAt(i + 1, cslist);
                maxq += alignScoreVectors(qfrac[i], qfrac[i]);
            }
            tfrac = new float[target.getLength()][];
            for (int i = 0; i < tfrac.length; i++) {
                tfrac[i] = target.getCompoundWeightsAt(i + 1, cslist);
                maxt += alignScoreVectors(tfrac[i], tfrac[i]);
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
