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
 */

package org.biojava.nbio.alignment;

/**
 * Summary statistics of a pairwise sequence alignment, computed by {@link PairwiseStatisticsAligner} without
 * building the alignment itself.
 * <p>
 * Positions are 1-based indices into the original sequences, as in the rest of BioJava (for example
 * {@link org.biojava.nbio.core.alignment.template.SequencePair#getIndexInQueryAt(int)}): the start is the first
 * residue of the sequence in the alignment and the end the last one, both inclusive. If the alignment is empty
 * (a local alignment where no pair of residues scores above zero), the length, the number of identical positions
 * and all four positions are 0.
 * <p>
 * Instances are immutable.
 */
public final class PairwiseAlignmentStatistics {

	private final int score, numIdenticals, length, queryStart, queryEnd, targetStart, targetEnd;

	/**
	 * @param score the score of the alignment
	 * @param numIdenticals the number of columns of the alignment with identical compounds
	 * @param length the number of columns of the alignment, including those with a gap
	 * @param queryStart the first aligned position of the query, 1-based, or 0 for an empty alignment
	 * @param queryEnd the last aligned position of the query, 1-based, or 0 for an empty alignment
	 * @param targetStart the first aligned position of the target, 1-based, or 0 for an empty alignment
	 * @param targetEnd the last aligned position of the target, 1-based, or 0 for an empty alignment
	 */
	public PairwiseAlignmentStatistics(int score, int numIdenticals, int length, int queryStart, int queryEnd,
			int targetStart, int targetEnd) {
		this.score = score;
		this.numIdenticals = numIdenticals;
		this.length = length;
		this.queryStart = queryStart;
		this.queryEnd = queryEnd;
		this.targetStart = targetStart;
		this.targetEnd = targetEnd;
	}

	/**
	 * @return the score of the alignment, as {@link SmithWaterman#getScore()} or {@link NeedlemanWunsch#getScore()}
	 *         give it
	 */
	public int getScore() {
		return score;
	}

	/**
	 * @return the number of columns of the alignment whose two compounds are identical, ignoring case, as
	 *         {@link org.biojava.nbio.core.alignment.template.SequencePair#getNumIdenticals()} counts them
	 */
	public int getNumIdenticals() {
		return numIdenticals;
	}

	/**
	 * @return the number of columns of the alignment, including those with a gap in either sequence
	 */
	public int getLength() {
		return length;
	}

	/**
	 * @return the first position of the query in the alignment, 1-based, or 0 if the alignment is empty
	 */
	public int getQueryStart() {
		return queryStart;
	}

	/**
	 * @return the last position of the query in the alignment, 1-based and inclusive, or 0 if the alignment is
	 *         empty
	 */
	public int getQueryEnd() {
		return queryEnd;
	}

	/**
	 * @return the first position of the target in the alignment, 1-based, or 0 if the alignment is empty
	 */
	public int getTargetStart() {
		return targetStart;
	}

	/**
	 * @return the last position of the target in the alignment, 1-based and inclusive, or 0 if the alignment is
	 *         empty
	 */
	public int getTargetEnd() {
		return targetEnd;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (!(o instanceof PairwiseAlignmentStatistics)) {
			return false;
		}
		PairwiseAlignmentStatistics s = (PairwiseAlignmentStatistics) o;
		return score == s.score && numIdenticals == s.numIdenticals && length == s.length
				&& queryStart == s.queryStart && queryEnd == s.queryEnd
				&& targetStart == s.targetStart && targetEnd == s.targetEnd;
	}

	@Override
	public int hashCode() {
		int h = score;
		h = 31 * h + numIdenticals;
		h = 31 * h + length;
		h = 31 * h + queryStart;
		h = 31 * h + queryEnd;
		h = 31 * h + targetStart;
		return 31 * h + targetEnd;
	}

	@Override
	public String toString() {
		return "PairwiseAlignmentStatistics[score=" + score + ", numIdenticals=" + numIdenticals + ", length="
				+ length + ", query=" + queryStart + "-" + queryEnd + ", target=" + targetStart + "-" + targetEnd
				+ "]";
	}
}
