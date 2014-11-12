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

package org.biojava3.alignment;

import org.biojava3.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

import java.util.List;

/**
 * Needleman and Wunsch defined an algorithm for pairwise global sequence alignments (best match of sections from each
 * {@link Sequence}).  This class performs such global sequence comparisons efficiently by dynamic programming, using
 * memory proportional to the length of the first sequence times the length of the second.
 *
 * @author Douglas Myers-Turnbull
 * @param <S> each {@link Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 * @see org.biojava3.alignment.NeedlemanWunsch For linear-space global pairwise alignment
 */
public class NeedlemanWunschQuadratic<S extends Sequence<C>, C extends Compound> extends
        AbstractPairwiseSequenceAligner<S, C> {

    /**
     * Before running a pairwise global sequence alignment, data must be sent in via calls to
     * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(GapPenalty)}, and
     * {@link #setSubstitutionMatrix(SubstitutionMatrix)}.
     */
    public NeedlemanWunschQuadratic() {
        super(null, null, null, null, false);
    }

    /**
     * Prepares for a pairwise global sequence alignment.
     *
     * @param query the first {@link Sequence} of the pair to align
     * @param target the second {@link Sequence} of the pair to align
     * @param gapPenalty the gap penalties used during alignment
     * @param subMatrix the set of substitution scores used during alignment
     */
    public NeedlemanWunschQuadratic(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        super(query, target, gapPenalty, subMatrix, false);
    }

    // method for AbstractMatrixAligner

    @Override
    protected void setProfile(List<Step> sx, List<Step> sy) {
        profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, xyStart[0],
                getQuery().getLength() - xyMax[0], sy, xyStart[1], getTarget().getLength() - xyMax[1]);
    }

}
