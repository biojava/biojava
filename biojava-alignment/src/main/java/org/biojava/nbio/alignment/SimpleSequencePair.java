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
 * Created on June 14, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.AlignedSequence.Step;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.List;

/**
 * Implements a data structure for the results of pairwise sequence alignment.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public class SimpleSequencePair<S extends Sequence<C>, C extends Compound> extends SimpleProfile<S, C>
implements SequencePair<S, C> {

	private int identicals = -1, similars = -1;

	/**
	 * Creates a pair profile for the given already aligned sequences.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @throws IllegalArgumentException if sequences differ in size
	 */
	public SimpleSequencePair(AlignedSequence<S, C> query, AlignedSequence<S, C> target) {
		super(query, target);
	}

	/**
	 * Creates a pair profile for the given sequences with a global alignment.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @param sx lists whether the query sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @param sy lists whether the target sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @throws IllegalArgumentException if alignments differ in size or given sequences do not fit in alignments
	 */
	public SimpleSequencePair(S query, S target, List<Step> sx, List<Step> sy) {
		this(query, target, sx, 0, 0, sy, 0, 0);
	}

	/**
	 * Creates a pair profile for the given sequences with a local alignment.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @param sx lists whether the query sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @param xb number of {@link Compound}s skipped in the query sequence before the aligned region
	 * @param xa number of {@link Compound}s skipped in the query sequence after the aligned region
	 * @param sy lists whether the target sequence aligns a {@link Compound} or gap at each index of the alignment
	 * @param yb number of {@link Compound}s skipped in the target sequence before the aligned region
	 * @param ya number of {@link Compound}s skipped in the target sequence after the aligned region
	 * @throws IllegalArgumentException if alignments differ in size or given sequences do not fit in alignments
	 */
	public SimpleSequencePair(S query, S target, List<Step> sx, int xb, int xa, List<Step> sy, int yb, int ya) {
		super(query, target, sx, xb, xa, sy, yb, ya);
	}

	@Override
	public C getCompoundInQueryAt(int alignmentIndex) {
		return getAlignedSequence(1).getCompoundAt(alignmentIndex);
	}

	@Override
	public C getCompoundInTargetAt(int alignmentIndex) {
		return getAlignedSequence(2).getCompoundAt(alignmentIndex);
	}

	@Override
	public int getIndexInQueryAt(int alignmentIndex) {
		return getAlignedSequence(1).getSequenceIndexAt(alignmentIndex);
	}

	@Override
	public int getIndexInQueryForTargetAt(int targetIndex) {
		return getAlignedSequence(1).getSequenceIndexAt(getAlignedSequence(2).getAlignmentIndexAt(targetIndex));
	}

	@Override
	public int getIndexInTargetAt(int alignmentIndex) {
		return getAlignedSequence(2).getSequenceIndexAt(alignmentIndex);
	}

	@Override
	public int getIndexInTargetForQueryAt(int queryIndex) {
		return getAlignedSequence(2).getSequenceIndexAt(getAlignedSequence(1).getAlignmentIndexAt(queryIndex));
	}

	@Override
	public int getNumIdenticals() {
		if (identicals == -1) {
			identicals = 0;
			for (int i = 1; i <= getLength(); i++) {
				if (getCompoundInQueryAt(i).equalsIgnoreCase(getCompoundInTargetAt(i))) {
					identicals++;
				}
			}
			getQuery().clearCache();
			getTarget().clearCache();
		}
		return identicals;
	}

	@Override
	public int getNumSimilars() {
		if (similars == -1) {
			similars = 0;
			for (int i = 1; i <= getLength(); i++) {

				C c1 = getCompoundInQueryAt(i);
				C c2 = getCompoundInTargetAt(i);

				if ( c1 instanceof AminoAcidCompound && c2 instanceof AminoAcidCompound) {
					short value = matrix.getValue((AminoAcidCompound)c1, (AminoAcidCompound)c2);
					if ( value > 0)
						similars++;
				} else {

					if (getCompoundSet().compoundsEquivalent(c1,c2)) {
						similars++;
					}
				}
			}
			getQuery().clearCache();
			getTarget().clearCache();
		}
		return similars;
	}

	@Override
	public AlignedSequence<S, C> getQuery() {
		return getAlignedSequence(1);
	}

	@Override
	public AlignedSequence<S, C> getTarget() {
		return getAlignedSequence(2);
	}

}
