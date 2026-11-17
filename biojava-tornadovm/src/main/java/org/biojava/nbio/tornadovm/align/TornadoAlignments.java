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
package org.biojava.nbio.tornadovm.align;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.GuideTree;
import org.biojava.nbio.alignment.template.AbstractScorer;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.tornadovm.TornadoSupport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;

/**
 * GPU versions of the all-pairs scoring entry points of {@link Alignments}.
 * <p>
 * Supports the {@link PairwiseSequenceScorerType#GLOBAL GLOBAL} (Needleman-Wunsch) and
 * {@link PairwiseSequenceScorerType#LOCAL LOCAL} (Smith-Waterman) score types with affine, constant
 * and linear gap penalties; the scores are exactly BioJava's. Large inputs run on the GPU; small inputs, and runs without a TornadoVM runtime,
 * run the same score-only algorithm on all CPU cores ({@link CpuAlignmentScores}). Other score types (identities and
 * similarities need a traceback) use {@link Alignments}.
 */
public final class TornadoAlignments {

	private static final Logger logger = LoggerFactory.getLogger(TornadoAlignments.class);

	/** Below this many DP cells in total the CPU implementation is used. */
	public static final long DEFAULT_MIN_GPU_CELLS = 10_000_000L;

	/** Device memory budget for the DP row buffers of one batch, in bytes. */
	static final long BUFFER_BUDGET = Long.getLong("biojava.tornado.sw.bufferBytes", 4L << 30);

	/** The shared GPU engines (affine and linear gaps); executions are serialised on the device. */
	private static final Engine AFFINE = new Engine(false);
	private static final Engine LINEAR = new Engine(true);

	private TornadoAlignments() {
	}

	/** @return true if the given scorer type is computed on the GPU */
	public static boolean isSupported(PairwiseSequenceScorerType type, GapPenalty gapPenalty) {
		return type == PairwiseSequenceScorerType.GLOBAL || type == PairwiseSequenceScorerType.LOCAL;
	}

	/**
	 * Drop-in for {@link Alignments#getAllPairsScores(List, PairwiseSequenceScorerType, GapPenalty,
	 * SubstitutionMatrix)}: the score of every pair (i, j), i &lt; j, in the same order.
	 */
	public static <S extends Sequence<C>, C extends Compound> double[] getAllPairsScores(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		if (!isSupported(type, gapPenalty)) {
			return Alignments.getAllPairsScores(sequences, type, gapPenalty, subMatrix);
		}
		int[] scores = computeScores(sequences, type == PairwiseSequenceScorerType.LOCAL, gapPenalty, subMatrix);
		double[] out = new double[scores.length];
		for (int k = 0; k < scores.length; k++) {
			out[k] = scores[k];
		}
		return out;
	}

	/**
	 * Drop-in for {@link Alignments#getAllPairsScorers} followed by {@link Alignments#runPairwiseScorers}: scorers
	 * whose scores are already computed, with the same {@link PairwiseSequenceScorer#getMaxScore() max} and
	 * {@link PairwiseSequenceScorer#getMinScore() min} (hence distances) as BioJava's aligners. They can be passed
	 * to {@link GuideTree#GuideTree(List, List)}.
	 */
	public static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getAllPairsScorers(
			List<S> sequences, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
			SubstitutionMatrix<C> subMatrix) {
		if (!isSupported(type, gapPenalty)) {
			List<PairwiseSequenceScorer<S, C>> scorers = Alignments.getAllPairsScorers(sequences, type, gapPenalty,
					subMatrix);
			Alignments.runPairwiseScorers(scorers);
			return scorers;
		}
		boolean local = type == PairwiseSequenceScorerType.LOCAL;
		int[] scores = computeScores(sequences, local, gapPenalty, subMatrix);
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
				// as AbstractPairwiseSequenceAligner.reset
				int max = Math.max(selfScores[i], selfScores[j]);
				int min = local ? 0 : (int) (2 * gapPenalty.getOpenPenalty()
						+ (query.getLength() + target.getLength()) * gapPenalty.getExtensionPenalty());
				scorers.add(new PrecomputedScorer<>(query, target, scores[k++], max, min));
			}
		}
		return scorers;
	}

	/** A {@link PairwiseSequenceScorer} holding an already computed score. */
	public static final class PrecomputedScorer<S extends Sequence<C>, C extends Compound> extends AbstractScorer
			implements PairwiseSequenceScorer<S, C> {
		private final S query;
		private final S target;
		private final int score;
		private final int max;
		private final int min;

		PrecomputedScorer(S query, S target, int score, int max, int min) {
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

	private static <S extends Sequence<C>, C extends Compound> long totalCells(List<S> sequences) {
		long sum = 0, sumSq = 0;
		for (S s : sequences) {
			sum += s.getLength();
			sumSq += (long) s.getLength() * s.getLength();
		}
		return (sum * sum - sumSq) / 2;
	}

	/**
	 * Scores all pairs: on the GPU for large inputs, else with {@link CpuAlignmentScores}. Same result either way.
	 */
	static <S extends Sequence<C>, C extends Compound> int[] computeScores(List<S> sequences, boolean local,
			GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		Encoded enc = new Encoded(sequences, subMatrix);
		boolean linear = gapPenalty.getType() == GapPenalty.Type.LINEAR;
		// BioJava's linear aligners ignore the open penalty (it is 0 for a LINEAR SimpleGapPenalty)
		int gop = linear ? 0 : gapPenalty.getOpenPenalty();
		int gep = gapPenalty.getExtensionPenalty();
		if (enc.nPairs() > 0 && TornadoSupport.useGpu("alignment", totalCells(sequences), DEFAULT_MIN_GPU_CELLS)) {
			try {
				synchronized (AFFINE) {
					return (linear ? LINEAR : AFFINE).run(enc, gop, gep, local);
				}
			} catch (RuntimeException | Error e) {
				TornadoSupport.disable("alignment", e);
			}
		}
		return CpuAlignmentScores.allPairs(enc.seqs, enc.subs, enc.k, gop, gep, local);
	}

	/**
	 * Sequences as alphabet codes. The alphabet is every distinct compound, scored with the matrix's own
	 * {@code getValue}, so that unknown compounds get exactly the value BioJava gives them.
	 */
	static final class Encoded {
		final int[][] seqs;
		final int[] subs;
		final int k;
		final int maxLength;

		<S extends Sequence<C>, C extends Compound> Encoded(List<S> sequences, SubstitutionMatrix<C> subMatrix) {
			Map<C, Integer> codes = new HashMap<>();
			List<C> alphabet = new ArrayList<>();
			seqs = new int[sequences.size()][];
			int longest = 0;
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
				longest = Math.max(longest, enc.length);
			}
			maxLength = longest;
			k = alphabet.size();
			subs = new int[Math.max(1, k * k)];
			for (int a = 0; a < k; a++) {
				for (int b = 0; b < k; b++) {
					subs[a * k + b] = subMatrix.getValue(alphabet.get(a), alphabet.get(b));
				}
			}
		}

		int nPairs() {
			return seqs.length * (seqs.length - 1) / 2;
		}
	}

	/**
	 * A persistent execution plan with buffers sized for a capacity (residues, sequences, alphabet, sequence
	 * length), compiled once and re-executed for every batch of every call that fits; rebuilt larger when a call
	 * does not fit.
	 */
	private static final class Engine {
		private final boolean linear;
		private int residueCapacity;
		private int seqCapacity;
		private int subsCapacity;
		private int lengthCapacity;
		private int batch;
		private IntArray residues;
		private IntArray seqStart;
		private IntArray subs;
		private IntArray pairQuery;
		private IntArray pairTarget;
		private IntArray params;
		private IntArray scores;
		private WorkerGrid worker;
		private GridScheduler gridScheduler;
		private TornadoExecutionPlan plan;

		Engine(boolean linear) {
			this.linear = linear;
		}

		int[] run(Encoded enc, int gop, int gep, boolean local) {
			int nSeq = enc.seqs.length;
			int nPairs = enc.nPairs();
			int totalResidues = 0;
			for (int[] s : enc.seqs) {
				totalResidues += s.length;
			}
			if (plan == null || totalResidues > residueCapacity || nSeq + 1 > seqCapacity
					|| enc.subs.length > subsCapacity || enc.maxLength + 1 > lengthCapacity) {
				build(totalResidues, nSeq + 1, enc.subs.length, enc.maxLength + 1);
			}
			int pos = 0;
			for (int i = 0; i < nSeq; i++) {
				seqStart.set(i, pos);
				for (int r : enc.seqs[i]) {
					residues.set(pos++, r);
				}
			}
			seqStart.set(nSeq, pos);
			for (int i = 0; i < enc.subs.length; i++) {
				subs.set(i, enc.subs[i]);
			}
			params.set(SwKernels.P_ALPHABET, enc.k);
			params.set(SwKernels.P_GOP, gop);
			params.set(SwKernels.P_GEP, gep);
			params.set(SwKernels.P_LOCAL, local ? 1 : 0);
			params.set(SwKernels.P_STRIDE, batch);

			// pairs in BioJava's order, then sorted by decreasing DP size so that a warp does similar work
			int[] queryOf = new int[nPairs];
			int[] targetOf = new int[nPairs];
			long[] keyed = new long[nPairs];
			int idx = 0;
			for (int i = 0; i < nSeq; i++) {
				for (int j = i + 1; j < nSeq; j++) {
					queryOf[idx] = i;
					targetOf[idx] = j;
					long cells = (long) enc.seqs[i].length * enc.seqs[j].length;
					// high bits: DP size (descending), low 32 bits: pair index
					keyed[idx] = (Integer.MAX_VALUE - Math.min(cells, Integer.MAX_VALUE)) << 32 | idx;
					idx++;
				}
			}
			Arrays.sort(keyed);
			int[] order = new int[nPairs];
			for (int p = 0; p < nPairs; p++) {
				order[p] = (int) keyed[p];
			}

			int[] result = new int[nPairs];
			for (int from = 0; from < nPairs; from += batch) {
				int size = Math.min(batch, nPairs - from);
				for (int p = 0; p < size; p++) {
					pairQuery.set(p, queryOf[order[from + p]]);
					pairTarget.set(p, targetOf[order[from + p]]);
				}
				params.set(SwKernels.P_PAIRS, size);
				worker.setGlobalWork(roundUp(size), 1, 1);
				long t0 = System.nanoTime();
				plan.withGridScheduler(gridScheduler).execute();
				if (Boolean.getBoolean("biojava.tornado.trace")) {
					logger.info("sw: batch of {} pairs in {} ms", size, (System.nanoTime() - t0) / 1e6);
				}
				for (int p = 0; p < size; p++) {
					result[order[from + p]] = scores.get(p);
				}
			}
			return result;
		}

		private void build(int minResidues, int minSeqs, int minSubs, int minLength) {
			if (plan != null) {
				try {
					plan.close();
				} catch (Exception e) {
					logger.warn("Could not release the previous alignment execution plan", e);
				}
			}
			residueCapacity = Math.max(residueCapacity, pow2(Math.max(minResidues, 1 << 16)));
			seqCapacity = Math.max(seqCapacity, pow2(Math.max(minSeqs, 1 << 12)));
			subsCapacity = Math.max(subsCapacity, pow2(Math.max(minSubs, 1024)));
			lengthCapacity = Math.max(lengthCapacity, (Math.max(minLength, 128) + 127) / 128 * 128);
			batch = (int) Math.max(SwKernels.GROUP_SIZE, BUFFER_BUDGET / (3L * 4 * lengthCapacity));
			batch = batch / SwKernels.GROUP_SIZE * SwKernels.GROUP_SIZE;

			residues = new IntArray(residueCapacity);
			seqStart = new IntArray(seqCapacity);
			subs = new IntArray(subsCapacity);
			pairQuery = new IntArray(batch);
			pairTarget = new IntArray(batch);
			params = new IntArray(6);
			scores = new IntArray(batch);
			IntArray bufM = new IntArray(lengthCapacity * batch);
			IntArray bufS0 = new IntArray(lengthCapacity * batch);
			IntArray bufS1 = new IntArray(lengthCapacity * batch);
			KernelContext ctx = new KernelContext();

			TaskGraph taskGraph = new TaskGraph("sw")
					.transferToDevice(DataTransferMode.EVERY_EXECUTION, residues, seqStart, subs, pairQuery,
							pairTarget, params)
					.task("scores", linear ? SwKernels::alignmentScoresLinear : SwKernels::alignmentScores, ctx,
							residues, seqStart, pairQuery, pairTarget, subs, params, bufM, bufS0, bufS1, scores)
					.transferToHost(DataTransferMode.EVERY_EXECUTION, scores);
			worker = new WorkerGrid1D(batch);
			worker.setLocalWork(SwKernels.GROUP_SIZE, 1, 1);
			gridScheduler = new GridScheduler("sw.scores", worker);
			plan = new TornadoExecutionPlan(taskGraph.snapshot());
		}
	}

	private static int roundUp(int n) {
		return (n + SwKernels.GROUP_SIZE - 1) / SwKernels.GROUP_SIZE * SwKernels.GROUP_SIZE;
	}

	private static int pow2(int n) {
		return n <= 1 ? 1 : Integer.highestOneBit(n - 1) << 1;
	}
}
