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

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.List;

/**
 * Smith and Waterman defined an algorithm for pairwise local sequence alignments (best match of sections from each
 * {@link Sequence}).  This class performs such local sequence comparisons efficiently by dynamic programming.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SmithWaterman<S extends Sequence<C>, C extends Compound> extends AbstractPairwiseSequenceAligner<S, C> {

    /**
     * Before running a pairwise local sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public SmithWaterman() {
        super(null, null, null, null, true);
    }

    /**
     * Prepares for a pairwise local sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public SmithWaterman(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix, true);
    }

    // method for AbstractMatrixAligner

    @Override
    protected void setProfile(List<Step> sx, List<Step> sy) {
        profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, xyStart[0],
                getQuery().getLength() - xyMax[0], sy, xyStart[1], getTarget().getLength() - xyMax[1]);
    }

}
