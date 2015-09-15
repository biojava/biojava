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

package org.biojava.nbio.alignment.template;

import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.alignment.template.GapPenalty.Type;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Implements common code for an {@link Aligner} for a pair of {@link Profile}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the pair of alignment {@link Profile}s is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public abstract class AbstractProfileProfileAligner<S extends Sequence<C>, C extends Compound>
        extends AbstractMatrixAligner<S, C> implements ProfileProfileAligner<S, C> {

	private final static Logger logger = LoggerFactory.getLogger(AbstractProfileProfileAligner.class);

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

    // method for ProfileProfileAligner

    @Override
    public ProfilePair<S, C> getPair() {
        if (pair == null) {
            align();
        }
        return pair;
    }

    // methods for ProfileProfileScorer

    @Override
    public Profile<S, C> getQuery() {
        return query;
    }

    @Override
    public Profile<S, C> getTarget() {
        return target;
    }

    // methods for AbstractMatrixAligner

    @Override
    protected CompoundSet<C> getCompoundSet() {
        return (query == null) ? null : query.getCompoundSet();
    }

    @Override
    protected List<C> getCompoundsOfQuery() {
        // TODO replace with consensus sequence
        return (query == null) ? new ArrayList<C>() : query.getAlignedSequence(1).getAsList();
    }

    @Override
    protected List<C> getCompoundsOfTarget() {
        // TODO replace with consensus sequence
        return (target == null) ? new ArrayList<C>() : target.getAlignedSequence(1).getAsList();
    }

    @Override
    protected int[] getScoreMatrixDimensions() {
        return new int[] { query.getLength() + 1, target.getLength() + 1, (getGapPenalty().getType() == Type.LINEAR) ?
                1 : 3 };
    }

    @Override
    protected int getSubstitutionScore(int queryColumn, int targetColumn) {
        return getSubstitutionScore(qfrac[queryColumn - 1], tfrac[targetColumn - 1]);
    }

    @Override
    protected boolean isReady() {
        // TODO when added to ConcurrencyTools, log completions and exceptions instead of printing stack traces
        try {
            if (query == null && queryFuture != null) {
                query = queryFuture.get();
            }
            if (target == null && targetFuture != null) {
                target = targetFuture.get();
            }
            reset();
        } catch (InterruptedException e) {
            logger.error("Interrupted Exception: ", e);
        } catch (ExecutionException e) {
            logger.error("Execution Exception: ", e);
        }
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
            cslist = query.getCompoundSet().getAllCompounds();
            qfrac = new float[query.getLength()][];
            for (int i = 0; i < qfrac.length; i++) {
                qfrac[i] = query.getCompoundWeightsAt(i + 1, cslist);
                maxq += getSubstitutionScore(qfrac[i], qfrac[i]);
            }
            tfrac = new float[target.getLength()][];
            for (int i = 0; i < tfrac.length; i++) {
                tfrac[i] = target.getCompoundWeightsAt(i + 1, cslist);
                maxt += getSubstitutionScore(tfrac[i], tfrac[i]);
            }
            max = (int) Math.max(maxq, maxt);
            score = min = isLocal() ? 0 : (int) (2 * getGapPenalty().getOpenPenalty() + (query.getLength() +
                    target.getLength()) * getGapPenalty().getExtensionPenalty());
        }
    }

    // helper method that scores alignment of two column vectors
    private int getSubstitutionScore(float[] qv, float[] tv) {
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
        return (int) Math.round(score);
    }

}
