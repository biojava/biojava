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

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.template.AbstractScorer;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * The scores of all sequence pairs for {@link PairwiseSequenceScorerType#GLOBAL GLOBAL} and
 * {@link PairwiseSequenceScorerType#LOCAL LOCAL} scoring, without building alignments.
 * <p>
 * A score needs only the last row of the dynamic programming matrix, so each pair is scored in linear memory on
 * primitive arrays, with no traceback. The recurrences are those of
 * {@link org.biojava.nbio.alignment.routines.AlignerHelper AlignerHelper} that {@link NeedlemanWunsch} and
 * {@link SmithWaterman} use, so the scores are exactly theirs, as are the maximum and minimum scores (and hence
 * the distances and similarities) of the scorers returned by {@link #getScorers}.
 * <p>
 * The pairs are scored concurrently on the shared thread pool of {@link ConcurrencyTools}.
 *
 * @author Michalis Papadimitriou
 */
final class AllPairsScores {

	private AllPairsScores() {
	}

	/**
	 * Whether the scores of {@code sequences} can be computed here: {@code GLOBAL} or {@code LOCAL} scoring, no
	 * empty sequence, and one compound set (as the aligners require to compute the maximum and minimum scores).
	 */
	static <S extends Sequence<C>, C extends Compound> boolean isSupported(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		if ((type != PairwiseSequenceScorerType.GLOBAL && type != PairwiseSequenceScorerType.LOCAL)
				|| gapPenalty == null || subMatrix == null || sequences.size() < 2) {
			return false;
		}
		for (S s : sequences) {
			if (s.getLength() == 0 || !s.getCompoundSet().equals(sequences.get(0).getCompoundSet())) {
				return false;
			}
		}
		return true;
	}

	/**
	 * The score of every pair (i, j), i &lt; j, in the order of {@link Alignments#getAllPairsScorers}.
	 */
	static <S extends Sequence<C>, C extends Compound> double[] getScores(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		int[] scores = computeScores(sequences, type == PairwiseSequenceScorerType.LOCAL, gapPenalty, subMatrix);
		double[] out = new double[scores.length];
		for (int k = 0; k < scores.length; k++) {
			out[k] = scores[k];
		}
		return out;
	}

	/**
	 * Scorers whose scores are already computed, in the order of {@link Alignments#getAllPairsScorers}, with the
	 * maximum and minimum scores that {@link NeedlemanWunsch} and {@link SmithWaterman} give the same pair.
	 */
	static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getScorers(
			List<S> sequences, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
			SubstitutionMatrix<C> subMatrix) {
		boolean local = type == PairwiseSequenceScorerType.LOCAL;
		int[] scores = computeScores(sequences, local, gapPenalty, subMatrix);
		// as AbstractPairwiseSequenceAligner.reset
		int[] selfScores = new int[sequences.size()];
		for (int i = 0; i < sequences.size(); i++) {
			for (C c : sequences.get(i)) {
				selfScores[i] += subMatrix.getValue(c, c);
			}
		}
		List<PairwiseSequenceScorer<S, C>> scorers = new ArrayList<>(scores.length);
		int k = 0;
		for (int i = 0; i < sequences.size(); i++) {
			for (int j = i + 1; j < sequences.size(); j++) {
				S query = sequences.get(i);
				S target = sequences.get(j);
				int max = Math.max(selfScores[i], selfScores[j]);
				int min = local ? 0 : (int) (2 * gapPenalty.getOpenPenalty()
						+ (query.getLength() + target.getLength()) * gapPenalty.getExtensionPenalty());
				scorers.add(new ComputedScorer<>(query, target, scores[k++], max, min));
			}
		}
		return scorers;
	}

	/** A {@link PairwiseSequenceScorer} holding a score computed by {@link AllPairsScores}. */
	static final class ComputedScorer<S extends Sequence<C>, C extends Compound> extends AbstractScorer
			implements PairwiseSequenceScorer<S, C> {

		private final S query;
		private final S target;
		private final int score;
		private final int max;
		private final int min;

		ComputedScorer(S query, S target, int score, int max, int min) {
			this.query = query;
			this.target = target;
			this.score = score;
			this.max = max;
			this.min = min;
		}

		@Override
		public S getQuery() {
			return query;
		}

		@Override
		public S getTarget() {
			return target;
		}

		@Override
		public double getMaxScore() {
			return max;
		}

		@Override
		public double getMinScore() {
			return min;
		}

		@Override
		public double getScore() {
			return score;
		}
	}

	private static <S extends Sequence<C>, C extends Compound> int[] computeScores(List<S> sequences, boolean local,
			GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		// The sequences as codes of their distinct compounds, and the matrix's own score for every pair of them,
		// so that compounds missing from the matrix score exactly as getValue scores them.
		Map<C, Integer> codes = new HashMap<>();
		List<C> alphabet = new ArrayList<>();
		int[][] seqs = new int[sequences.size()][];
		for (int i = 0; i < seqs.length; i++) {
			S seq = sequences.get(i);
			int[] enc = new int[seq.getLength()];
			int pos = 0;
			for (C c : seq) {
				Integer code = codes.get(c);
				if (code == null) {
					code = alphabet.size();
					codes.put(c, code);
					alphabet.add(c);
				}
				enc[pos++] = code;
			}
			seqs[i] = enc;
		}
		int k = alphabet.size();
		int[] subs = new int[k * k];
		for (int a = 0; a < k; a++) {
			for (int b = 0; b < k; b++) {
				subs[a * k + b] = subMatrix.getValue(alphabet.get(a), alphabet.get(b));
			}
		}
		// as AbstractMatrixAligner.align: the linear recurrence for LINEAR gap penalties, the affine one otherwise
		boolean linear = gapPenalty.getType() == GapPenalty.Type.LINEAR;
		int gop = gapPenalty.getOpenPenalty();
		int gep = gapPenalty.getExtensionPenalty();

		int n = seqs.length;
		int[] rowStart = new int[n];
		for (int i = 1; i < n; i++) {
			rowStart[i] = rowStart[i - 1] + (n - i);
		}
		int[] scores = new int[n * (n - 1) / 2];
		// one task per query sequence: the pairs (i, j) for every j > i
		List<Future<Void>> futures = new ArrayList<>(n);
		for (int row = 0; row < n - 1; row++) {
			final int i = row;
			futures.add(ConcurrencyTools.submit(() -> {
				for (int j = i + 1; j < n; j++) {
					scores[rowStart[i] + j - i - 1] = linear
							? scoreLinear(seqs[i], seqs[j], subs, k, gep, local)
							: scoreAffine(seqs[i], seqs[j], subs, k, gop, gep, local);
				}
				return null;
			}, String.format("Scoring the pairs of sequence %d of %d", i + 1, n)));
		}
		for (Future<Void> f : futures) {
			try {
				f.get();
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new IllegalStateException("Interrupted while scoring sequence pairs", e);
			} catch (ExecutionException e) {
				throw new IllegalStateException("Failed to score sequence pairs", e.getCause());
			}
		}
		return scores;
	}

	/**
	 * The affine-gap score of query {@code q} against target {@code t}, as the three-state
	 * {@code AlignerHelper.setScoreVector}: substitution, gap in the target and gap in the query. Only the
	 * previous row of each state is kept.
	 */
	static int scoreAffine(int[] q, int[] t, int[] subs, int k, int gop, int gep, boolean local) {
		// as AlignerHelper: the lowest score that can still be extended without overflow
		int min = Integer.MIN_VALUE - gop - gep;
		int n = t.length;
		int[] maxRow = new int[n + 1];
		int[] s0 = new int[n + 1];
		int[] s1 = new int[n + 1];
		if (!local) {
			s1[0] = gop;
			int s2 = gop;
			for (int y = 1; y <= n; y++) {
				s2 += gep;
				s0[y] = min;
				s1[y] = min;
				maxRow[y] = s2;
			}
		}
		int best = 0;
		for (int x = 1; x <= q.length; x++) {
			int row = q[x - 1] * k;
			int diag = maxRow[0];
			int left0;
			int left2;
			if (local) {
				left0 = 0;
				left2 = 0;
				maxRow[0] = 0;
				s0[0] = 0;
				s1[0] = 0;
			} else {
				int v = s1[0] + gep;
				left0 = min;
				left2 = min;
				s0[0] = min;
				s1[0] = v;
				maxRow[0] = Math.max(v, min);
			}
			for (int y = 1; y <= n; y++) {
				int up0 = s0[y];
				int up1 = s1[y];
				int upMax = maxRow[y];
				int a = diag + subs[row + t[y - 1]];
				int b = up1 >= up0 + gop ? up1 + gep : up0 + gop + gep;
				int c = left0 + gop >= left2 ? left0 + gop + gep : left2 + gep;
				if (local) {
					a = Math.max(a, 0);
					b = Math.max(b, 0);
					c = Math.max(c, 0);
					best = Math.max(best, a);
				}
				s0[y] = a;
				s1[y] = b;
				maxRow[y] = Math.max(Math.max(a, b), c);
				diag = upMax;
				left0 = a;
				left2 = c;
			}
		}
		return local ? best : maxRow[n];
	}

	/** The linear-gap score of query {@code q} against target {@code t}, as the one-state {@code setScoreVector}. */
	static int scoreLinear(int[] q, int[] t, int[] subs, int k, int gep, boolean local) {
		int n = t.length;
		int[] h = new int[n + 1];
		if (!local) {
			for (int y = 1; y <= n; y++) {
				h[y] = h[y - 1] + gep;
			}
		}
		int best = 0;
		for (int x = 1; x <= q.length; x++) {
			int row = q[x - 1] * k;
			int diag = h[0];
			if (!local) {
				h[0] = h[0] + gep;
			}
			int left = h[0];
			for (int y = 1; y <= n; y++) {
				int up = h[y];
				int deletion = up + gep;
				int insertion = left + gep;
				int substitution = diag + subs[row + t[y - 1]];
				int v = deletion >= substitution && deletion >= insertion ? deletion
						: (substitution >= insertion ? substitution : insertion);
				if (local) {
					v = Math.max(v, 0);
					best = Math.max(best, v);
				}
				h[y] = v;
				diag = up;
				left = v;
			}
		}
		return local ? best : h[n];
	}
}
