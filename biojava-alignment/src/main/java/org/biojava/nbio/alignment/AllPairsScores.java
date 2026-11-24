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
import org.biojava.nbio.alignment.PairwiseStatisticsAligner.EncodedSequence;
import org.biojava.nbio.alignment.template.AbstractScorer;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.ConcurrencyTools;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * The scores of all sequence pairs for {@link PairwiseSequenceScorerType#GLOBAL GLOBAL},
 * {@link PairwiseSequenceScorerType#LOCAL LOCAL}, {@link PairwiseSequenceScorerType#GLOBAL_IDENTITIES
 * GLOBAL_IDENTITIES} and {@link PairwiseSequenceScorerType#LOCAL_IDENTITIES LOCAL_IDENTITIES} scoring, without
 * building alignments.
 * <p>
 * Each pair is aligned by a {@link PairwiseStatisticsAligner}, in linear memory and with no traceback. Its score is
 * that of {@link NeedlemanWunsch} or {@link SmithWaterman}, and its number of identical positions and length are
 * those of the alignment they return, so the scorers returned by {@link #getScorers} have the same scores, maximum
 * and minimum scores (and hence distances and similarities) as the {@link PairwiseSequenceScorer}s of
 * {@link Alignments#getAllPairsScorers}.
 * <p>
 * Each sequence is encoded once. The pairs are scored concurrently on the shared thread pool of
 * {@link ConcurrencyTools}, one task per query sequence.
 *
 * @author Michalis Papadimitriou
 */
final class AllPairsScores {

	private AllPairsScores() {
	}

	/**
	 * Whether the scores of {@code sequences} can be computed here: {@code GLOBAL}, {@code LOCAL},
	 * {@code GLOBAL_IDENTITIES} or {@code LOCAL_IDENTITIES} scoring, gap penalties that are zero or negative, no
	 * empty sequence, one compound set (as the aligners require to compute the maximum and minimum scores), and
	 * every compound in the compound set of the substitution matrix (as {@link PairwiseStatisticsAligner} requires).
	 */
	static <S extends Sequence<C>, C extends Compound> boolean isSupported(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		if (!isSupported(type) || gapPenalty == null || subMatrix == null || sequences.size() < 2
				|| gapPenalty.getOpenPenalty() > 0 || gapPenalty.getExtensionPenalty() > 0) {
			return false;
		}
		CompoundSet<C> compoundSet = sequences.get(0).getCompoundSet();
		Set<C> matrixCompounds = new HashSet<>(subMatrix.getCompoundSet().getAllCompounds());
		for (S s : sequences) {
			if (s.getLength() == 0 || !s.getCompoundSet().equals(compoundSet)) {
				return false;
			}
			for (C c : s) {
				if (!matrixCompounds.contains(c)) {
					return false;
				}
			}
		}
		return true;
	}

	private static boolean isSupported(PairwiseSequenceScorerType type) {
		return type == PairwiseSequenceScorerType.GLOBAL || type == PairwiseSequenceScorerType.LOCAL
				|| type == PairwiseSequenceScorerType.GLOBAL_IDENTITIES
				|| type == PairwiseSequenceScorerType.LOCAL_IDENTITIES;
	}

	private static boolean isLocal(PairwiseSequenceScorerType type) {
		return type == PairwiseSequenceScorerType.LOCAL || type == PairwiseSequenceScorerType.LOCAL_IDENTITIES;
	}

	private static boolean isIdentities(PairwiseSequenceScorerType type) {
		return type == PairwiseSequenceScorerType.GLOBAL_IDENTITIES
				|| type == PairwiseSequenceScorerType.LOCAL_IDENTITIES;
	}

	/**
	 * The score of every pair (i, j), i &lt; j, in the order of {@link Alignments#getAllPairsScorers}: the alignment
	 * score, or for the {@code *_IDENTITIES} types the number of identical positions.
	 */
	static <S extends Sequence<C>, C extends Compound> double[] getScores(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		PairwiseAlignmentStatistics[] stats = computeStatistics(sequences, isLocal(type), gapPenalty, subMatrix);
		boolean identities = isIdentities(type);
		double[] out = new double[stats.length];
		for (int k = 0; k < stats.length; k++) {
			out[k] = identities ? stats[k].getNumIdenticals() : stats[k].getScore();
		}
		return out;
	}

	/**
	 * Scorers whose scores are already computed, in the order of {@link Alignments#getAllPairsScorers}, with the
	 * scores, maximum and minimum scores of the scorers that it would create for the same pair: those of
	 * {@link NeedlemanWunsch} and {@link SmithWaterman}, or of {@link FractionalIdentityScorer} (the number of
	 * identical positions, out of the alignment length) for the {@code *_IDENTITIES} types.
	 */
	static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getScorers(
			List<S> sequences, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
			SubstitutionMatrix<C> subMatrix) {
		boolean local = isLocal(type);
		PairwiseAlignmentStatistics[] stats = computeStatistics(sequences, local, gapPenalty, subMatrix);
		List<PairwiseSequenceScorer<S, C>> scorers = new ArrayList<>(stats.length);
		if (isIdentities(type)) {
			// as FractionalIdentityScorer
			int k = 0;
			for (int i = 0; i < sequences.size(); i++) {
				for (int j = i + 1; j < sequences.size(); j++) {
					PairwiseAlignmentStatistics s = stats[k++];
					scorers.add(new ComputedScorer<>(sequences.get(i), sequences.get(j), s.getNumIdenticals(),
							s.getLength(), 0));
				}
			}
			return scorers;
		}
		// as AbstractPairwiseSequenceAligner.reset
		int[] selfScores = new int[sequences.size()];
		for (int i = 0; i < sequences.size(); i++) {
			for (C c : sequences.get(i)) {
				selfScores[i] += subMatrix.getValue(c, c);
			}
		}
		int k = 0;
		for (int i = 0; i < sequences.size(); i++) {
			for (int j = i + 1; j < sequences.size(); j++) {
				S query = sequences.get(i);
				S target = sequences.get(j);
				int max = Math.max(selfScores[i], selfScores[j]);
				int min = local ? 0 : (int) (2 * gapPenalty.getOpenPenalty()
						+ (query.getLength() + target.getLength()) * gapPenalty.getExtensionPenalty());
				scorers.add(new ComputedScorer<>(query, target, stats[k++].getScore(), max, min));
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

	/** The statistics of every pair (i, j), i &lt; j, in the order of {@link Alignments#getAllPairsScorers}. */
	private static <S extends Sequence<C>, C extends Compound> PairwiseAlignmentStatistics[] computeStatistics(
			List<S> sequences, boolean local, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		PairwiseStatisticsAligner<C> aligner = new PairwiseStatisticsAligner<>(gapPenalty, subMatrix, local);
		int n = sequences.size();
		List<EncodedSequence<C>> encoded = new ArrayList<>(n);
		for (S s : sequences) {
			encoded.add(aligner.encode(s));
		}
		int[] rowStart = new int[n];
		for (int i = 1; i < n; i++) {
			rowStart[i] = rowStart[i - 1] + (n - i);
		}
		PairwiseAlignmentStatistics[] stats = new PairwiseAlignmentStatistics[n * (n - 1) / 2];
		// one task per query sequence: the pairs (i, j) for every j > i
		List<Future<Void>> futures = new ArrayList<>(n);
		for (int row = 0; row < n - 1; row++) {
			final int i = row;
			futures.add(ConcurrencyTools.submit(() -> {
				for (int j = i + 1; j < n; j++) {
					stats[rowStart[i] + j - i - 1] = aligner.align(encoded.get(i), encoded.get(j));
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
		return stats;
	}
}
