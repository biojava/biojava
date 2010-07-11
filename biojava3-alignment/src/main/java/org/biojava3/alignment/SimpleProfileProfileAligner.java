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
 * Created on June 30, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.util.List;
import java.util.concurrent.Future;

import org.biojava3.alignment.template.AbstractProfileProfileAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.Aligner;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.ProfilePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements a simple (naive) {@link Aligner} for a pair of {@link Profile}s.  This is basically an extension of the
 * {@link NeedlemanWunsch} pairwise sequence aligner to pairwise profile alignment using a sum-of-pairs score.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the pair of alignment {@link Profile}s is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SimpleProfileProfileAligner<S extends Sequence<C>, C extends Compound>
        extends AbstractProfileProfileAligner<S, C> {

    /**
     * Before running a profile-profile alignment, data must be sent in via calls to
     * {@link #setQuery(Profile)}, {@link #setTarget(Profile)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public SimpleProfileProfileAligner() {
    }

    /**
     * Prepares for a profile-profile alignment.
     *
     * @param query the first {@link Profile} of the pair to align
     * @param target the second {@link Profile} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SimpleProfileProfileAligner(Profile<S, C> query, Profile<S, C> target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
    }

    /**
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align, still to be calculated
     * @param target the second {@link Profile} of the pair to align, still to be calculated
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SimpleProfileProfileAligner(Future<ProfilePair<S, C>> query, Future<ProfilePair<S, C>> target,
            GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
    }

    /**
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align
     * @param target the second {@link Profile} of the pair to align, still to be calculated
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SimpleProfileProfileAligner(Profile<S, C> query, Future<ProfilePair<S, C>> target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
    }

    /**
     * Prepares for a profile-profile alignment run concurrently.
     *
     * @param query the first {@link Profile} of the pair to align, still to be calculated
     * @param target the second {@link Profile} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SimpleProfileProfileAligner(Future<ProfilePair<S, C>> query, Profile<S, C> target, GapPenalty gapPenalty,
            SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix);
    }

    // helper methods

    // scores with linear gap penalty; saves memory by skipping allocation of separate matching and gap arrays
    @Override
    protected void alignScoreLinear() {
        for (int x = 1; x < scores.length; x++) {
            scores[x][0] = (short) (scores[x - 1][0] + getGapPenalty().getExtensionPenalty());
        }
        for (int y = 1; y < scores[0].length; y++) {
            scores[0][y] = (short) (scores[0][y - 1] + getGapPenalty().getExtensionPenalty());
        }
        for (int x = 1; x < scores.length; x++) {
            for (int y = 1; y < scores[0].length; y++) {
                scores[x][y] = (short) Math.max(Math.max(scores[x - 1][y] + getGapPenalty().getExtensionPenalty(),
                        scores[x][y - 1] + getGapPenalty().getExtensionPenalty()), scores[x - 1][y - 1] +
                        alignScoreColumns(x, y));
            }
        }
    }

    // traces back through score matrix; chooses highroad alignment
    @Override
    protected void alignTracebackLinear(List<Step> sx, List<Step> sy) {
        int x = scores.length - 1, y = scores[0].length - 1;
        while (x > 0 || y > 0) {
            if (x == 0) {
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
            } else if (y == 0 || scores[x][y] == scores[x - 1][y] + getGapPenalty().getExtensionPenalty()) {
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
            } else if (scores[x][y] == scores[x - 1][y - 1] + alignScoreColumns(x, y)) {
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.COMPOUND);
                x--;
                y--;
            } else {
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
            }
        }
    }

    // scores with affine gap penalty
    @Override
    protected void alignScoreAffine(short[][] ix, short[][] iy) {
        GapPenalty gapPenalty = getGapPenalty();
        short min = (short) (Short.MIN_VALUE - gapPenalty.getOpenPenalty() - gapPenalty.getExtensionPenalty());
        ix[0][0] = iy[0][0] = gapPenalty.getOpenPenalty();
        for (int x = 1; x < scores.length; x++) {
            scores[x][0] = iy[x][0] = min;
            ix[x][0] = (short) (ix[x - 1][0] + gapPenalty.getExtensionPenalty());
        }
        for (int y = 1; y < scores[0].length; y++) {
            scores[0][y] = ix[0][y] = min;
            iy[0][y] = (short) (iy[0][y - 1] + gapPenalty.getExtensionPenalty());
        }
        for (int x = 1; x < scores.length; x++) {
            for (int y = 1; y < scores[0].length; y++) {
                scores[x][y] = (short) (Math.max(Math.max(scores[x - 1][y - 1], ix[x - 1][y - 1]), iy[x - 1][y - 1]) +
                        alignScoreColumns(x, y));
                ix[x][y] = (short) (Math.max(scores[x - 1][y] + gapPenalty.getOpenPenalty(), ix[x - 1][y]) +
                        gapPenalty.getExtensionPenalty());
                iy[x][y] = (short) (Math.max(scores[x][y - 1] + gapPenalty.getOpenPenalty(), iy[x][y - 1]) +
                        gapPenalty.getExtensionPenalty());
            }
        }
    }

    // traces back through score matrices; chooses highroad alignment
    @Override
    protected void alignTracebackAffine(List<Step> sx, List<Step> sy, short[][] ix, short[][] iy) {
        int x = scores.length - 1, y = scores[0].length - 1;
        int max = Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
        Last last = (max == ix[x][y]) ? Last.IX : ((max == scores[x][y]) ? Last.M : Last.IY);
        while (x > 0 || y > 0) {
            switch (last) {
            case IX:
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.GAP);
                x--;
                last = (scores[x][y] + getGapPenalty().getOpenPenalty() > ix[x][y]) ? Last.M : Last.IX;
                break;
            case M:
                sx.add(0, Step.COMPOUND);
                sy.add(0, Step.COMPOUND);
                x--;
                y--;
                max = Math.max(Math.max(scores[x][y], ix[x][y]), iy[x][y]);
                last = (max == ix[x][y]) ? Last.IX : ((max == scores[x][y]) ? Last.M : Last.IY);
                break;
            case IY:
                sx.add(0, Step.GAP);
                sy.add(0, Step.COMPOUND);
                y--;
                last = (scores[x][y] + getGapPenalty().getOpenPenalty() >= iy[x][y]) ? Last.M : Last.IY;
            }
        }
    }

    // sets output fields
    @Override
    protected void alignSetOutputs(List<Step> sx, List<Step> sy) {
        score = scores[scores.length - 1][scores[0].length - 1];
        profile = pair = new SimpleProfilePair<S, C>(getQuery(), getTarget(), sx, sy);
    }

}
