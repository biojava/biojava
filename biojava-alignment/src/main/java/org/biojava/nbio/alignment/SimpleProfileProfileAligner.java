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

package org.biojava.nbio.alignment;

import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.alignment.SimpleProfilePair;
import org.biojava.nbio.alignment.template.*;
import org.biojava.nbio.core.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.List;
import java.util.concurrent.Future;

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

    // method for AbstractMatrixAligner

    @Override
    protected void setProfile(List<Step> sx, List<Step> sy) {
        profile = pair = new SimpleProfilePair<S, C>(getQuery(), getTarget(), sx, sy);
    }

}
