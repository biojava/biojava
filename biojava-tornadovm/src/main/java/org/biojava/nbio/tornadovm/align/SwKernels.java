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

import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;

/**
 * GPU kernel for batches of pairwise alignment scores with affine gaps, as BioJava's {@code NeedlemanWunsch}
 * (global) and {@code SmithWaterman} (local) aligners compute them in {@code AlignerHelper.setScorePoint} and
 * {@code setScoreVector}.
 * <p>
 * Score only (no traceback), in linear memory, one pair per thread, strip-mined. Integer arithmetic, so the scores are exactly
 * BioJava's on every device.
 */
public final class SwKernels {

	/** Threads per work-group. */
	public static final int GROUP_SIZE = 64;

	/** Indices into the {@code params} array. */
	public static final int P_ALPHABET = 0, P_GOP = 1, P_GEP = 2, P_LOCAL = 3, P_PAIRS = 4, P_STRIDE = 5;

	private SwKernels() {
	}

	/** Query rows computed per pass over the target by {@link #alignmentScores}. */
	public static final int STRIP = 32;

	/**
	 * Computes the alignment score of pair {@code p = globalIdx} of the batch: query {@code pairQuery[p]} against
	 * target {@code pairTarget[p]}.
	 * <p>
	 * The DP boundary row is kept per target position {@code y} in three interleaved buffers, at {@code y*stride+p}
	 * so that the threads of a warp access consecutive words: the substitution state {@code S0}, the deletion state
	 * {@code S1} and the maximum of the three states {@code M} (all the substitution recurrence needs). The DP is
	 * strip-mined: each pass over the target computes {@link #STRIP} query rows, keeping their per-row state (left
	 * cell and diagonal) in local memory, so the boundary buffers are read and written once per strip, not per row.
	 *
	 * @param ctx kernel context, one thread per pair
	 * @param residues all sequences, as alphabet codes, concatenated
	 * @param seqStart start of each sequence in residues, size nSeq+1
	 * @param pairQuery query sequence index of each pair
	 * @param pairTarget target sequence index of each pair
	 * @param subs substitution scores, {@code subs[a*K+b]} for query code a and target code b
	 * @param params [K, gop, gep, local (0/1), nPairs, stride]; gop and gep are negative, as in BioJava
	 * @param bufM scratch, max of the three states of the boundary row
	 * @param bufS0 scratch, substitution state of the boundary row
	 * @param bufS1 scratch, deletion state of the boundary row
	 * @param scores output, one score per pair
	 */
	public static void alignmentScores(KernelContext ctx, IntArray residues, IntArray seqStart,
			IntArray pairQuery, IntArray pairTarget, IntArray subs, IntArray params, IntArray bufM, IntArray bufS0,
			IntArray bufS1, IntArray scores) {
		int p = ctx.globalIdx;
		int lid = ctx.localIdx;
		int[] left0 = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int[] left2 = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int[] diag = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int[] qrow = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int nPairs = params.get(P_PAIRS);
		if (p < nPairs) {
			int alphabet = params.get(P_ALPHABET);
			int gop = params.get(P_GOP);
			int gep = params.get(P_GEP);
			int local = params.get(P_LOCAL);
			int stride = params.get(P_STRIDE);
			int min = Integer.MIN_VALUE - gop - gep;

			int qBegin = seqStart.get(pairQuery.get(p));
			int m = seqStart.get(pairQuery.get(p) + 1) - qBegin;
			int tBegin = seqStart.get(pairTarget.get(p));
			int n = seqStart.get(pairTarget.get(p) + 1) - tBegin;

			// row 0
			int s1col0 = gop;
			if (local == 1) {
				for (int y = 0; y <= n; y++) {
					bufM.set(y * stride + p, 0);
					bufS0.set(y * stride + p, 0);
					bufS1.set(y * stride + p, 0);
				}
			} else {
				bufS0.set(p, 0);
				bufS1.set(p, gop);
				bufM.set(p, gop > 0 ? gop : 0);
				int s2 = gop;
				for (int y = 1; y <= n; y++) {
					s2 = s2 + gep;
					bufS0.set(y * stride + p, min);
					bufS1.set(y * stride + p, min);
					bufM.set(y * stride + p, s2 > min ? s2 : min);
				}
			}

			int best = 0;
			for (int x0 = 0; x0 < m; x0 += STRIP) {
				int rows = m - x0 < STRIP ? m - x0 : STRIP;
				// column 0 of the strip rows, and the diagonal of their first cell
				diag[lid] = bufM.get(p);
				for (int r = 0; r < rows; r++) {
					int slot = r * GROUP_SIZE + lid;
					qrow[slot] = residues.get(qBegin + x0 + r) * alphabet;
					if (local == 1) {
						left0[slot] = 0;
						left2[slot] = 0;
						if (r > 0) {
							diag[slot] = 0;
						}
					} else {
						// column 0 of row x0+r+1: S1 = S1 of the row above + gep, S0 = S2 = min
						int above = s1col0;
						s1col0 = s1col0 + gep;
						left0[slot] = min;
						left2[slot] = min;
						if (r > 0) {
							diag[slot] = above > min ? above : min;
						}
					}
				}
				// the bottom row's column 0 becomes the next strip's top boundary
				if (local == 1) {
					bufM.set(p, 0);
					bufS0.set(p, 0);
					bufS1.set(p, 0);
				} else {
					bufS0.set(p, min);
					bufS1.set(p, s1col0);
					bufM.set(p, s1col0 > min ? s1col0 : min);
				}

				for (int y = 1; y <= n; y++) {
					int idx = y * stride + p;
					int up0 = bufS0.get(idx);
					int up1 = bufS1.get(idx);
					int upM = bufM.get(idx);
					int t = residues.get(tBegin + y - 1);
					for (int r = 0; r < rows; r++) {
						int slot = r * GROUP_SIZE + lid;
						int s0 = diag[slot] + subs.get(qrow[slot] + t);
						int s1 = up1 >= up0 + gop ? up1 + gep : up0 + gop + gep;
						int l0 = left0[slot];
						int l2 = left2[slot];
						int s2 = l0 + gop >= l2 ? l0 + gop + gep : l2 + gep;
						if (local == 1) {
							s0 = s0 <= 0 ? 0 : s0;
							s1 = s1 <= 0 ? 0 : s1;
							s2 = s2 <= 0 ? 0 : s2;
							best = s0 > best ? s0 : best;
						}
						int mx = s0 >= s1 ? s0 : s1;
						mx = mx >= s2 ? mx : s2;
						diag[slot] = upM;
						left0[slot] = s0;
						left2[slot] = s2;
						up0 = s0;
						up1 = s1;
						upM = mx;
					}
					bufS0.set(idx, up0);
					bufS1.set(idx, up1);
					bufM.set(idx, upM);
				}
			}

			if (local == 1) {
				scores.set(p, best);
			} else {
				scores.set(p, bufM.get(n * stride + p));
			}
		}
	}

	/**
	 * Same as {@link #alignmentScores} for a linear gap penalty (open penalty 0), as BioJava's linear
	 * {@code AlignerHelper.setScorePoint}: a single state, {@code H = max(up + gep, left + gep, diag + sub)}. Only
	 * {@code bufM} is used as the boundary row.
	 */
	public static void alignmentScoresLinear(KernelContext ctx, IntArray residues, IntArray seqStart,
			IntArray pairQuery, IntArray pairTarget, IntArray subs, IntArray params, IntArray bufM, IntArray bufS0,
			IntArray bufS1, IntArray scores) {
		int p = ctx.globalIdx;
		int lid = ctx.localIdx;
		int[] left = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int[] diag = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int[] qrow = ctx.allocateIntLocalArray(STRIP * GROUP_SIZE);
		int nPairs = params.get(P_PAIRS);
		if (p < nPairs) {
			int alphabet = params.get(P_ALPHABET);
			int gep = params.get(P_GEP);
			int local = params.get(P_LOCAL);
			int stride = params.get(P_STRIDE);

			int qBegin = seqStart.get(pairQuery.get(p));
			int m = seqStart.get(pairQuery.get(p) + 1) - qBegin;
			int tBegin = seqStart.get(pairTarget.get(p));
			int n = seqStart.get(pairTarget.get(p) + 1) - tBegin;

			// row 0: zeros (local), or 0, gep, 2*gep, ... (global)
			int h = 0;
			for (int y = 0; y <= n; y++) {
				bufM.set(y * stride + p, h);
				if (local == 0) {
					h = h + gep;
				}
			}

			int col0 = 0;
			int best = 0;
			for (int x0 = 0; x0 < m; x0 += STRIP) {
				int rows = m - x0 < STRIP ? m - x0 : STRIP;
				for (int r = 0; r < rows; r++) {
					int slot = r * GROUP_SIZE + lid;
					qrow[slot] = residues.get(qBegin + x0 + r) * alphabet;
					// diagonal of the first cell: column 0 of the row above; left: column 0 of this row
					diag[slot] = col0;
					if (local == 0) {
						col0 = col0 + gep;
					}
					left[slot] = col0;
				}
				bufM.set(p, col0);

				for (int y = 1; y <= n; y++) {
					int idx = y * stride + p;
					int up = bufM.get(idx);
					int t = residues.get(tBegin + y - 1);
					for (int r = 0; r < rows; r++) {
						int slot = r * GROUP_SIZE + lid;
						int d = up + gep;
						int ins = left[slot] + gep;
						int sub = diag[slot] + subs.get(qrow[slot] + t);
						int v = d >= sub && d >= ins ? d : (sub >= ins ? sub : ins);
						if (local == 1) {
							v = v <= 0 ? 0 : v;
							best = v > best ? v : best;
						}
						diag[slot] = up;
						left[slot] = v;
						up = v;
					}
					bufM.set(idx, up);
				}
			}
			scores.set(p, local == 1 ? best : bufM.get(n * stride + p));
		}
	}
}
