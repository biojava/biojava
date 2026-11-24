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

import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Aligns pairs of sequences, locally like {@link SmithWaterman} or globally like {@link NeedlemanWunsch}, and
 * returns only the {@link PairwiseAlignmentStatistics statistics} of the alignment: its score, number of identical
 * positions, length and where it starts and ends. It does not build the alignment, so it needs no traceback
 * matrix: memory is linear in the length of the target and time is proportional to the product of the lengths.
 * This suits all-against-all comparisons of many sequences, for example to cluster them.
 * <p>
 * <b>Results.</b> The score is the same as the one of {@link SmithWaterman} or {@link NeedlemanWunsch} with the same
 * gap penalty and substitution matrix: the same recurrence is computed, with the same integer arithmetic. The
 * statistics are carried forward through the dynamic programming along the same choices of predecessor that those
 * aligners record as traceback pointers, so they are those of the very alignment that their traceback returns,
 * also where several alignments have the optimal score. In particular, for a local alignment:
 * <ul>
 * <li>the alignment ends at the first cell with the best score, scanning the query positions in order and, for
 * each of them, the target positions in order;</li>
 * <li>where moves tie, the alignment prefers a gap in the target (deletion), then a substitution, then a gap in
 * the query (insertion), as the traceback pointers of those aligners do;</li>
 * <li>it starts and ends with an aligned pair of compounds, never with a gap.</li>
 * </ul>
 * <b>Gap penalties.</b> The penalties are those of the {@link GapPenalty}, taken when this aligner is created:
 * {@link GapPenalty#getOpenPenalty()} and {@link GapPenalty#getExtensionPenalty()}, both zero or negative. A gap
 * of length n adds open + n &times; extension to the score: with {@code new SimpleGapPenalty(11, 1)}, a gap of 1
 * costs 12 and a gap of 3 costs 14. With an open penalty of 0 ({@link GapPenalty.Type#LINEAR}) a gap of length n
 * costs n &times; extension, and as in the other aligners a gap in one sequence may then directly follow a gap in
 * the other; with a nonzero open penalty it may not.
 * <p>
 * <b>Compounds.</b> Compounds are identical as {@link Compound#equalsIgnoreCase(Compound)} says, which is how
 * {@link org.biojava.nbio.core.alignment.template.SequencePair#getNumIdenticals()} counts them. Every compound of
 * an aligned sequence must be in the {@link org.biojava.nbio.core.sequence.template.CompoundSet} of the
 * substitution matrix: when this aligner is created, it looks up the scores of all pairs of compounds of that set
 * once.
 * <p>
 * <b>Many pairs.</b> To align a sequence against many others, {@link #encode(Sequence) encode} each sequence once
 * and align the encoded sequences:
 * <pre>
 * PairwiseStatisticsAligner&lt;AminoAcidCompound&gt; aligner = new PairwiseStatisticsAligner&lt;&gt;(
 *         new SimpleGapPenalty(11, 1), SubstitutionMatrixHelper.getBlosum62(), true);
 * List&lt;PairwiseStatisticsAligner.EncodedSequence&lt;AminoAcidCompound&gt;&gt; encoded = new ArrayList&lt;&gt;();
 * for (ProteinSequence s : sequences) {
 *     encoded.add(aligner.encode(s));
 * }
 * PairwiseAlignmentStatistics stats = aligner.align(encoded.get(i), encoded.get(j));
 * </pre>
 * <b>Thread safety.</b> An aligner is immutable once created, and so are the sequences it encodes: one aligner
 * and its encoded sequences can be shared by any number of threads, each alignment allocating only its own rows.
 *
 * @param <C> each element of a sequence is a {@link Compound} of type C
 */
public final class PairwiseStatisticsAligner<C extends Compound> {

	/**
	 * A sequence encoded by a {@link PairwiseStatisticsAligner} for aligning it with that aligner. It is immutable.
	 *
	 * @param <C> each element of the sequence is a {@link Compound} of type C
	 */
	public static final class EncodedSequence<C extends Compound> {

		private final PairwiseStatisticsAligner<C> aligner;
		private final int[] codes;

		private EncodedSequence(PairwiseStatisticsAligner<C> aligner, int[] codes) {
			this.aligner = aligner;
			this.codes = codes;
		}

		/**
		 * @return the number of compounds of the sequence
		 */
		public int getLength() {
			return codes.length;
		}
	}

	// the statistics of an alignment are carried in a long: identicals in the high 32 bits, length in the low ones
	private static final long ONE_IDENTICAL = 1L << 32;

	// the start of a local alignment is carried in a long: query position in the high 32 bits, target position in
	// the low ones

	private final boolean local, linear;
	private final int gop, gep;
	private final Map<C, Integer> codes;
	// substitution score of each pair of compound codes
	private final int[][] scores;
	// what aligning each pair of compound codes adds to the statistics: one column, and one identical if they are
	private final long[][] steps;

	/**
	 * Creates an aligner. Later changes to the gap penalty object do not affect it.
	 *
	 * @param gapPenalty the gap penalties used during alignment; both must be zero or negative, as
	 *        {@link SimpleGapPenalty} makes them
	 * @param subMatrix the set of substitution scores used during alignment
	 * @param local true for a local alignment, as {@link SmithWaterman} computes, false for a global one, as
	 *        {@link NeedlemanWunsch} computes
	 * @throws IllegalArgumentException if a gap penalty is positive
	 */
	public PairwiseStatisticsAligner(GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix, boolean local) {
		Objects.requireNonNull(gapPenalty, "gapPenalty");
		Objects.requireNonNull(subMatrix, "subMatrix");
		this.local = local;
		this.gop = gapPenalty.getOpenPenalty();
		this.gep = gapPenalty.getExtensionPenalty();
		if (gop > 0 || gep > 0) {
			throw new IllegalArgumentException("Gap penalties must be zero or negative, but they are " + gop
					+ " to open and " + gep + " to extend");
		}
		this.linear = gapPenalty.getType() == GapPenalty.Type.LINEAR;

		Map<C, Integer> codes = new HashMap<>();
		List<C> compounds = new ArrayList<>();
		for (C c : subMatrix.getCompoundSet().getAllCompounds()) {
			if (!codes.containsKey(c)) {
				codes.put(c, compounds.size());
				compounds.add(c);
			}
		}
		int k = compounds.size();
		scores = new int[k][k];
		steps = new long[k][k];
		for (int i = 0; i < k; i++) {
			C from = compounds.get(i);
			for (int j = 0; j < k; j++) {
				C to = compounds.get(j);
				scores[i][j] = subMatrix.getValue(from, to);
				steps[i][j] = from.equalsIgnoreCase(to) ? ONE_IDENTICAL + 1 : 1;
			}
		}
		this.codes = Collections.unmodifiableMap(codes);
	}

	/**
	 * @return true if this aligner computes local alignments, false if global ones
	 */
	public boolean isLocal() {
		return local;
	}

	/**
	 * Encodes a sequence for aligning it with this aligner. Encoding a sequence once and aligning it many times
	 * avoids looking up its compounds for every pair.
	 *
	 * @param sequence the sequence to encode
	 * @return the encoded sequence, which only this aligner can align
	 * @throws IllegalArgumentException if a compound of the sequence is not in the compound set of the
	 *         substitution matrix
	 */
	public EncodedSequence<C> encode(Sequence<C> sequence) {
		int[] encoded = new int[sequence.getLength()];
		int i = 0;
		for (C c : sequence) {
			Integer code = codes.get(c);
			if (code == null) {
				throw new IllegalArgumentException("Compound " + c + " at position " + (i + 1)
						+ " is not in the compound set of the substitution matrix");
			}
			encoded[i++] = code;
		}
		return new EncodedSequence<>(this, encoded);
	}

	/**
	 * Aligns two sequences. To align a sequence with many others, {@link #encode(Sequence) encode} it once and use
	 * {@link #align(EncodedSequence, EncodedSequence)} instead.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @return the statistics of the alignment
	 * @throws IllegalArgumentException if a compound of either sequence is not in the compound set of the
	 *         substitution matrix
	 */
	public PairwiseAlignmentStatistics align(Sequence<C> query, Sequence<C> target) {
		return align(encode(query), encode(target));
	}

	/**
	 * Aligns two sequences encoded by this aligner.
	 *
	 * @param query the first sequence of the pair
	 * @param target the second sequence of the pair
	 * @return the statistics of the alignment
	 * @throws IllegalArgumentException if a sequence was encoded by another aligner
	 */
	public PairwiseAlignmentStatistics align(EncodedSequence<C> query, EncodedSequence<C> target) {
		if (query.aligner != this || target.aligner != this) {
			throw new IllegalArgumentException("The sequences must be encoded by this aligner");
		}
		int[] q = query.codes, t = target.codes;
		if (local) {
			return linear ? alignLocalLinear(q, t) : alignLocalAffine(q, t);
		}
		return linear ? alignGlobalLinear(q, t) : alignGlobalAffine(q, t);
	}

	/*
	 * The recurrences below are those of AlignerHelper, with the same comparisons in the same order. With the
	 * affine gap penalties there are three states per cell: the alignment ending there ends with a substitution, a
	 * deletion (gap in the target) or an insertion (gap in the query). A substitution follows the best state of
	 * the cell on the diagonal, which ties prefer in the order deletion, substitution, insertion, as the traceback
	 * pointers of AlignerHelper.setScorePoint do.
	 *
	 * Only one row is kept for each state. While a row is computed, the arrays hold the current row up to the
	 * previous column and the previous row from the current column on. The cell on the diagonal and the cell to
	 * the left are kept in local variables.
	 */

	private PairwiseAlignmentStatistics alignLocalAffine(int[] q, int[] t) {
		final int n = q.length, m = t.length, gop = this.gop, gep = this.gep;
		int[] sub = new int[m + 1], del = new int[m + 1], best = new int[m + 1];
		long[] subStats = new long[m + 1], delStats = new long[m + 1], bestStats = new long[m + 1];
		long[] subStart = new long[m + 1], delStart = new long[m + 1], bestStart = new long[m + 1];
		int maxScore = 0, maxX = 0, maxY = 0;
		long maxStats = 0, maxStart = 0;

		// row 0 and column 0 are all 0: a local alignment can start anywhere
		for (int x = 1; x <= n; x++) {
			int[] scoreRow = scores[q[x - 1]];
			long[] stepRow = steps[q[x - 1]];
			long rowStart = (long) x << 32;
			int diag = 0, left = 0, ins = 0;
			long diagStats = 0, diagStart = 0, leftStats = 0, leftStart = 0, insStats = 0, insStart = 0;
			for (int y = 1; y <= m; y++) {
				int ty = t[y - 1];

				// substitution, after the best state on the diagonal, or starting the alignment here
				int s = diag + scoreRow[ty];
				long sStats, sStart;
				if (diag > 0) {
					sStats = diagStats + stepRow[ty];
					sStart = diagStart;
				} else {
					sStats = stepRow[ty];
					sStart = rowStart | y;
				}

				// deletion: extend one, or open one after a substitution
				int up = sub[y], upDel = del[y];
				int d;
				long dStats, dStart;
				if (upDel >= up + gop) {
					d = upDel + gep;
					dStats = delStats[y] + 1;
					dStart = delStart[y];
				} else {
					d = up + gop + gep;
					dStats = subStats[y] + 1;
					dStart = subStart[y];
				}

				// insertion: open one after a substitution, or extend one
				int i;
				long iStats, iStart;
				if (left + gop >= ins) {
					i = left + gop + gep;
					iStats = leftStats + 1;
					iStart = leftStart;
				} else {
					i = ins + gep;
					iStats = insStats + 1;
					iStart = insStart;
				}

				// nothing scoring zero or less is extended: a positive state always follows a positive one
				if (s < 0) {
					s = 0;
				}
				if (d < 0) {
					d = 0;
				}
				if (i < 0) {
					i = 0;
				}

				int b;
				long bStats, bStart;
				if (d >= s && d >= i) {
					b = d;
					bStats = dStats;
					bStart = dStart;
				} else if (s >= i) {
					b = s;
					bStats = sStats;
					bStart = sStart;
				} else {
					b = i;
					bStats = iStats;
					bStart = iStart;
				}

				diag = best[y];
				diagStats = bestStats[y];
				diagStart = bestStart[y];
				sub[y] = s;
				subStats[y] = sStats;
				subStart[y] = sStart;
				del[y] = d;
				delStats[y] = dStats;
				delStart[y] = dStart;
				best[y] = b;
				bestStats[y] = bStats;
				bestStart[y] = bStart;
				left = s;
				leftStats = sStats;
				leftStart = sStart;
				ins = i;
				insStats = iStats;
				insStart = iStart;

				if (s > maxScore) {
					maxScore = s;
					maxStats = sStats;
					maxStart = sStart;
					maxX = x;
					maxY = y;
				}
			}
		}
		return localStatistics(maxScore, maxStats, maxStart, maxX, maxY);
	}

	private PairwiseAlignmentStatistics alignLocalLinear(int[] q, int[] t) {
		final int n = q.length, m = t.length, gep = this.gep;
		int[] score = new int[m + 1];
		long[] stats = new long[m + 1], start = new long[m + 1];
		int maxScore = 0, maxX = 0, maxY = 0;
		long maxStats = 0, maxStart = 0;

		for (int x = 1; x <= n; x++) {
			int[] scoreRow = scores[q[x - 1]];
			long[] stepRow = steps[q[x - 1]];
			long rowStart = (long) x << 32;
			int diag = 0, left = 0;
			long diagStats = 0, diagStart = 0, leftStats = 0, leftStart = 0;
			for (int y = 1; y <= m; y++) {
				int ty = t[y - 1];
				int up = score[y];
				int d = up + gep, i = left + gep, s = diag + scoreRow[ty];
				int v;
				long vStats, vStart;
				if (d >= s && d >= i) {
					v = d;
					vStats = stats[y] + 1;
					vStart = start[y];
				} else if (s >= i) {
					v = s;
					if (diag > 0) {
						vStats = diagStats + stepRow[ty];
						vStart = diagStart;
					} else {
						vStats = stepRow[ty];
						vStart = rowStart | y;
					}
				} else {
					v = i;
					vStats = leftStats + 1;
					vStart = leftStart;
				}
				if (v < 0) {
					v = 0;
				}

				diag = up;
				diagStats = stats[y];
				diagStart = start[y];
				score[y] = v;
				stats[y] = vStats;
				start[y] = vStart;
				left = v;
				leftStats = vStats;
				leftStart = vStart;

				if (v > maxScore) {
					maxScore = v;
					maxStats = vStats;
					maxStart = vStart;
					maxX = x;
					maxY = y;
				}
			}
		}
		return localStatistics(maxScore, maxStats, maxStart, maxX, maxY);
	}

	private static PairwiseAlignmentStatistics localStatistics(int score, long stats, long start, int endX, int endY) {
		if (score <= 0) {
			return new PairwiseAlignmentStatistics(0, 0, 0, 0, 0, 0, 0);
		}
		return new PairwiseAlignmentStatistics(score, (int) (stats >>> 32), (int) stats, (int) (start >>> 32), endX,
				(int) start, endY);
	}

	private PairwiseAlignmentStatistics alignGlobalAffine(int[] q, int[] t) {
		final int n = q.length, m = t.length, gop = this.gop, gep = this.gep;
		// as in AlignerHelper: low enough to never be chosen, high enough not to overflow when penalties are added
		final int min = Integer.MIN_VALUE - gop - gep;
		int[] sub = new int[m + 1], del = new int[m + 1], best = new int[m + 1];
		long[] subStats = new long[m + 1], delStats = new long[m + 1], bestStats = new long[m + 1];

		// row 0: only insertions lead there from the origin
		sub[0] = 0;
		del[0] = gop;
		int ins = gop;
		long insStats = 0;
		setBest(best, bestStats, 0, sub[0], 0, del[0], 0, ins, 0);
		for (int y = 1; y <= m; y++) {
			sub[y] = del[y] = min;
			ins += gep;
			insStats += 1;
			setBest(best, bestStats, y, min, 0, min, 0, ins, insStats);
		}

		for (int x = 1; x <= n; x++) {
			int[] scoreRow = scores[q[x - 1]];
			long[] stepRow = steps[q[x - 1]];

			// column 0: only deletions lead there from the origin
			int diag = best[0];
			long diagStats = bestStats[0];
			sub[0] = min;
			del[0] += gep;
			delStats[0] += 1;
			ins = min;
			insStats = 0;
			setBest(best, bestStats, 0, min, 0, del[0], delStats[0], min, 0);
			int left = min;
			long leftStats = 0;

			for (int y = 1; y <= m; y++) {
				int ty = t[y - 1];

				int s = diag + scoreRow[ty];
				long sStats = diagStats + stepRow[ty];

				int up = sub[y], upDel = del[y];
				int d;
				long dStats;
				if (upDel >= up + gop) {
					d = upDel + gep;
					dStats = delStats[y] + 1;
				} else {
					d = up + gop + gep;
					dStats = subStats[y] + 1;
				}

				int i;
				long iStats;
				if (left + gop >= ins) {
					i = left + gop + gep;
					iStats = leftStats + 1;
				} else {
					i = ins + gep;
					iStats = insStats + 1;
				}

				int b;
				long bStats;
				if (d >= s && d >= i) {
					b = d;
					bStats = dStats;
				} else if (s >= i) {
					b = s;
					bStats = sStats;
				} else {
					b = i;
					bStats = iStats;
				}

				diag = best[y];
				diagStats = bestStats[y];
				sub[y] = s;
				subStats[y] = sStats;
				del[y] = d;
				delStats[y] = dStats;
				best[y] = b;
				bestStats[y] = bStats;
				left = s;
				leftStats = sStats;
				ins = i;
				insStats = iStats;
			}
		}

		// the state the traceback of AlignerHelper.setSteps starts from, and the score of AbstractMatrixAligner
		int s = sub[m], d = del[m];
		long stats;
		if (d > s && d > ins) {
			stats = delStats[m];
		} else if (s > ins) {
			stats = subStats[m];
		} else {
			stats = insStats;
		}
		return globalStatistics(Math.max(s, Math.max(d, ins)), stats, n, m);
	}

	// the best of the three states of a cell, preferring a deletion, then a substitution, then an insertion
	private static void setBest(int[] best, long[] bestStats, int y, int s, long sStats, int d, long dStats, int i,
			long iStats) {
		if (d >= s && d >= i) {
			best[y] = d;
			bestStats[y] = dStats;
		} else if (s >= i) {
			best[y] = s;
			bestStats[y] = sStats;
		} else {
			best[y] = i;
			bestStats[y] = iStats;
		}
	}

	private PairwiseAlignmentStatistics alignGlobalLinear(int[] q, int[] t) {
		final int n = q.length, m = t.length, gep = this.gep;
		int[] score = new int[m + 1];
		long[] stats = new long[m + 1];
		for (int y = 1; y <= m; y++) {
			score[y] = score[y - 1] + gep;
			stats[y] = y;
		}

		for (int x = 1; x <= n; x++) {
			int[] scoreRow = scores[q[x - 1]];
			long[] stepRow = steps[q[x - 1]];
			int diag = score[0];
			long diagStats = stats[0];
			score[0] += gep;
			stats[0] += 1;
			int left = score[0];
			long leftStats = stats[0];
			for (int y = 1; y <= m; y++) {
				int ty = t[y - 1];
				int up = score[y];
				int d = up + gep, i = left + gep, s = diag + scoreRow[ty];
				int v;
				long vStats;
				if (d >= s && d >= i) {
					v = d;
					vStats = stats[y] + 1;
				} else if (s >= i) {
					v = s;
					vStats = diagStats + stepRow[ty];
				} else {
					v = i;
					vStats = leftStats + 1;
				}
				diag = up;
				diagStats = stats[y];
				score[y] = v;
				stats[y] = vStats;
				left = v;
				leftStats = vStats;
			}
		}
		return globalStatistics(score[m], stats[m], n, m);
	}

	private static PairwiseAlignmentStatistics globalStatistics(int score, long stats, int n, int m) {
		return new PairwiseAlignmentStatistics(score, (int) (stats >>> 32), (int) stats, n > 0 ? 1 : 0, n,
				m > 0 ? 1 : 0, m);
	}
}
