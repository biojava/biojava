package org.biojava3.alignment;

import org.biojava3.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

import java.util.List;

/**
 * A pairwise sequence alignment that penalizes end gaps in the target but not the query.
 *
 * @author dmyersturnbull
 * @param <S> each {@link org.biojava3.core.sequence.template.Sequence} of the alignment pair is of type S
 * @param <C> each element of an {@link org.biojava3.alignment.template.AlignedSequence} is a {@link org.biojava3.core.sequence.template.Compound} of type C
 * @see org.biojava3.alignment.NeedlemanWunsch
 * @see org.biojava3.alignment.SmithWaterman
 * @see org.biojava3.alignment.AlignmentMode
 */
public class SemiglobalAlignment<S extends Sequence<C>, C extends Compound> extends AbstractPairwiseSequenceAligner<S, C> {

	/**
	 * Before running a pairwise semiglobal sequence alignment, data must be sent in via calls to
	 * {@link #setQuery(Sequence)}, {@link #setTarget(Sequence)}, {@link #setGapPenalty(org.biojava3.alignment.template.GapPenalty)}, and
	 * {@link #setSubstitutionMatrix(org.biojava3.alignment.template.SubstitutionMatrix)}.
	 */
	public SemiglobalAlignment() {
		super(null, null, null, null, AlignmentMode.LOCAL);
	}

	/**
	 * Prepares for a pairwise semiglobal sequence alignment.
	 *
	 * @param query the first {@link Sequence} of the pair to align
	 * @param target the second {@link Sequence} of the pair to align
	 * @param gapPenalty the gap penalties used during alignment
	 * @param subMatrix the set of substitution scores used during alignment
	 */
	public SemiglobalAlignment(S query, S target, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		super(query, target, gapPenalty, subMatrix, AlignmentMode.SEMIGLOBAL);
	}

	@Override
	protected void setProfile(List<AlignedSequence.Step> sx, List<AlignedSequence.Step> sy) {
		profile = pair = new SimpleSequencePair<S, C>(getQuery(), getTarget(), sx, xyStart[0],
		                                              getQuery().getLength() - xyMax[0], sy, xyStart[1], getTarget().getLength() - xyMax[1]); // TODO Verify
	}

}
