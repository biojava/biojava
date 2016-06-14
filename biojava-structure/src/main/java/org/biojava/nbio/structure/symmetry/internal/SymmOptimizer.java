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
package org.biojava.nbio.structure.symmetry.internal;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Optimizes a symmetry alignment by a Monte Carlo score optimization of the
 * repeat multiple alignment. The superposition of the repeats is not free
 * (felxible), because it is constrained on the symmetry axes found in the
 * structure. This is the main difference to the {@link MultipleMC} algorithm in
 * biojava. Another major difference is that the free Pool is shared for all
 * repeats, so that no residue can appear to more than one repeat at a time.
 * <p>
 * This algorithm does not use a unfiform distribution for selecting moves,
 * farther residues have more probability to be shrinked or gapped. This
 * modification of the algorithm improves convergence and running time.
 * <p>
 * Use call method to parallelize optimizations, or use optimize method instead.
 * Because gaps are allowed in the repeats, a {@link MultipleAlignment} format
 * is returned.
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SymmOptimizer {

	private static final Logger logger = LoggerFactory
			.getLogger(SymmOptimizer.class);

	private Random rnd;

	// Optimization parameters
	private int Rmin = 2; // min aligned repeats per column
	private int Lmin; // min repeat length
	private int maxIter; // max iterations
	private double C = 20; // probability of accept bad moves constant

	// Score function parameters
	private static final double Gopen = 20.0; // Penalty for opening gap
	private static final double Gextend = 10.0; // Penalty for extending gaps
	private double dCutoff;

	// Alignment Information
	private MultipleAlignment msa;
	private SymmetryAxes axes;
	private Atom[] atoms;
	private int order;
	private int length; // total alignment columns (block size)
	private int repeatCore; // core length (without gaps)

	// Aligned Residues and Score
	private List<List<Integer>> block; // residues aligned
	private List<Integer> freePool; // residues not aligned
	private double mcScore; // alignment score to optimize

	// Variables that store the history of the optimization - slower if on
	private static final boolean history = false;
	private static final int saveStep = 100;
	private static final String pathToHistory = "results/symm-opt/";
	private List<Long> timeHistory;
	private List<Integer> lengthHistory;
	private List<Double> rmsdHistory;
	private List<Double> tmScoreHistory;
	private List<Double> mcScoreHistory;

	/**
	 * Constructor with a seed MultipleAligment storing a refined symmetry
	 * alignment of the repeats. To perform the optimization use the call or
	 * optimize methods after instantiation.
	 *
	 * @param symmResult
	 *            CeSymmResult with all the information
	 * @throws RefinerFailedException
	 * @throws StructureException
	 */
	public SymmOptimizer(CeSymmResult symmResult) {

		this.axes = symmResult.getAxes();
		this.rnd = new Random(symmResult.getParams().getRndSeed());
		this.Lmin = symmResult.getParams().getMinCoreLength();
		this.dCutoff = symmResult.getParams().getDistanceCutoff();

		MultipleAlignmentEnsemble e = symmResult.getMultipleAlignment()
				.getEnsemble().clone();
		this.msa = e.getMultipleAlignment(0);

		this.atoms = msa.getAtomArrays().get(0);
		this.order = msa.size();
		this.repeatCore = msa.getCoreLength();

		// 50% of the structures aligned (minimum) or all (no gaps)
		if (symmResult.getParams().isGaps())
			Rmin = Math.max(order / 2, 2);
		else
			Rmin = order;

		maxIter = symmResult.getParams().getOptimizationSteps();
		if (maxIter < 1)
			maxIter = 100 * atoms.length;
	}

	private void initialize() throws StructureException, RefinerFailedException {

		if (order == 1)
			throw new RefinerFailedException(
					"Non-symmetric seed alignment: order = 1");
		if (repeatCore < 1)
			throw new RefinerFailedException(
					"Seed alignment too short: repeat core length < 1");

		// Initialize the history variables
		timeHistory = new ArrayList<Long>();
		lengthHistory = new ArrayList<Integer>();
		rmsdHistory = new ArrayList<Double>();
		mcScoreHistory = new ArrayList<Double>();
		tmScoreHistory = new ArrayList<Double>();

		C = 20 * order;

		// Initialize alignment variables
		block = msa.getBlock(0).getAlignRes();
		freePool = new ArrayList<Integer>();
		length = block.get(0).size();

		// Store the residues aligned in the block
		List<Integer> aligned = new ArrayList<Integer>();
		for (int su = 0; su < order; su++)
			aligned.addAll(block.get(su));

		// Add any residue not aligned to the free pool
		for (int i = 0; i < atoms.length; i++) {
			if (!aligned.contains(i))
				freePool.add(i);
		}
		checkGaps();

		// Set the MC score of the initial state (seed alignment)
		updateMultipleAlignment();
		mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend,
				dCutoff);
	}

	/**
	 * Optimization method based in a Monte-Carlo approach. Starting from the
	 * refined alignment uses 4 types of moves:
	 * <p>
	 * <li>1- Shift Row: if there are enough freePool residues available.
	 * <li>2- Expand Block: add another alignment column.
	 * <li>3- Shrink Block: move a block column to the freePool.
	 * <li>4- Insert gap: insert a gap in a position of the alignment.
	 *
	 * @throws StructureException
	 * @throws RefinerFailedException
	 *             if the alignment is not symmetric or too short.
	 */
	public MultipleAlignment optimize() throws StructureException,
			RefinerFailedException {

		initialize();

		// Save the optimal alignment
		List<List<Integer>> optBlock = new ArrayList<List<Integer>>();
		List<Integer> optFreePool = new ArrayList<Integer>();
		optFreePool.addAll(freePool);
		for (int k = 0; k < order; k++) {
			List<Integer> b = new ArrayList<Integer>();
			b.addAll(block.get(k));
			optBlock.add(b);
		}
		double optScore = mcScore;

		int conv = 0; // Number of steps without an alignment improvement
		int i = 1;
		int stepsToConverge = Math.max(maxIter / 50, 1000);
		long initialTime = System.nanoTime()/1000000;

		while (i < maxIter && conv < stepsToConverge) {

			// Save the state of the system
			List<List<Integer>> lastBlock = new ArrayList<List<Integer>>();
			List<Integer> lastFreePool = new ArrayList<Integer>();
			lastFreePool.addAll(freePool);
			for (int k = 0; k < order; k++) {
				List<Integer> b = new ArrayList<Integer>();
				b.addAll(block.get(k));
				lastBlock.add(b);
			}
			double lastScore = mcScore;
			int lastRepeatCore = repeatCore;

			boolean moved = false;

			while (!moved) {
				// Randomly select one of the steps to modify the alignment.
				// Because of biased moves, the probabilities are not the same
				double move = rnd.nextDouble();
				if (move < 0.4) {
					moved = shiftRow();
					logger.debug("did shift");
				} else if (move < 0.7) {
					moved = expandBlock();
					logger.debug("did expand");
				} else if (move < 0.85) {
					moved = shrinkBlock();
					logger.debug("did shrink");
				} else {
					moved = insertGap();
					logger.debug("did insert gap");
				}
			}

			// Get the properties of the new alignment
			updateMultipleAlignment();
			mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend,
					dCutoff);

			// Calculate change in the optimization Score
			double AS = mcScore - lastScore;
			double prob = 1.0;

			if (AS < 0) {

				// Probability of accepting bad move
				prob = probabilityFunction(AS, i, maxIter);
				double p = rnd.nextDouble();

				// Reject the move
				if (p > prob) {
					block = lastBlock;
					freePool = lastFreePool;
					length = block.get(0).size();
					repeatCore = lastRepeatCore;
					mcScore = lastScore;
					conv++; // no change in score if rejected

				} else
					conv = 0; // if accepted

			} else
				conv = 0; // if positive change

			logger.debug(i + ": --prob: " + prob + ", --score: " + AS
					+ ", --conv: " + conv);
			
			// Store as the optimal alignment if better
			if (mcScore > optScore) {
				optBlock = new ArrayList<List<Integer>>();
				optFreePool = new ArrayList<Integer>();
				optFreePool.addAll(freePool);
				for (int k = 0; k < order; k++) {
					List<Integer> b = new ArrayList<Integer>();
					b.addAll(block.get(k));
					optBlock.add(b);
				}
				optScore = mcScore;
			}

			if (history) {
				if (i % saveStep == 1) {
					// Get the correct superposition again
					updateMultipleAlignment();

					timeHistory.add(System.nanoTime()/1000000 - initialTime);
					lengthHistory.add(length);
					rmsdHistory.add(msa.getScore(MultipleAlignmentScorer.RMSD));
					tmScoreHistory.add(msa
							.getScore(MultipleAlignmentScorer.AVGTM_SCORE));
					mcScoreHistory.add(mcScore);
				}
			}

			i++;
		}
		
		// Use the optimal alignment of the trajectory
		block = optBlock;
		freePool = optFreePool;
		mcScore = optScore;
		
		// Superimpose and calculate final scores
		updateMultipleAlignment();
		msa.putScore(MultipleAlignmentScorer.MC_SCORE, mcScore);

		// Save the history to the results folder of the symmetry project
		if (history) {
			try {
				saveHistory(pathToHistory);
			} catch (Exception e) {
				logger.warn("Could not save history file: " + e.getMessage());
			}
		}

		return msa;
	}

	/**
	 * This method translates the internal data structures to a
	 * MultipleAlignment of the repeats in order to use the methods to score
	 * MultipleAlignments.
	 *
	 * @throws StructureException
	 * @throws RefinerFailedException
	 */
	private void updateMultipleAlignment() throws StructureException,
			RefinerFailedException {

		msa.clear();

		// Override the alignment with the new information
		Block b = msa.getBlock(0);
		b.setAlignRes(block);
		repeatCore = b.getCoreLength();
		if (repeatCore < 1)
			throw new RefinerFailedException(
					"Optimization converged to length 0");

		SymmetryTools.updateSymmetryTransformation(axes, msa);
	}

	/**
	 * Method that loops through all the alignment columns and checks that there
	 * are no more gaps than the maximum allowed: Rmin.
	 * <p>
	 * There must be at least Rmin residues different than null in every
	 * alignment column. In case there is a column with more gaps than allowed
	 * it will be shrinked (moved to freePool).
	 *
	 * @return true if any columns has been shrinked and false otherwise
	 */
	private boolean checkGaps() {

		List<Integer> shrinkColumns = new ArrayList<Integer>();
		// Loop for each column
		for (int res = 0; res < length; res++) {
			int gapCount = 0;
			// Loop for each repeat and count the gaps
			for (int su = 0; su < order; su++) {
				if (block.get(su).get(res) == null)
					gapCount++;
			}
			if ((order - gapCount) < Rmin) {
				// Add the column to the shrink list
				shrinkColumns.add(res);
			}
		}

		// Shrink the columns that have more gaps than allowed
		for (int col = shrinkColumns.size() - 1; col >= 0; col--) {
			for (int su = 0; su < order; su++) {
				Integer residue = block.get(su).get(shrinkColumns.get(col));
				block.get(su).remove((int) shrinkColumns.get(col));
				if (residue != null)
					freePool.add(residue);
				Collections.sort(freePool);
			}
			length--;
		}

		if (shrinkColumns.size() != 0)
			return true;
		else
			return false;
	}

	/**
	 * Insert a gap in one of the repeats into selected position (by higher
	 * distances) in the alignment. Calculates the average residue distance to
	 * make the choice. A gap is a null in the block.
	 *
	 * @throws StructureException
	 * @throws RefinerFailedException
	 */
	private boolean insertGap() throws StructureException,
			RefinerFailedException {

		// Let gaps only if the repeat is larger than the minimum length
		if (repeatCore <= Lmin)
			return false;

		// Select residue by maximum distance
		updateMultipleAlignment();
		Matrix residueDistances = MultipleAlignmentTools
				.getAverageResidueDistances(msa);

		double maxDist = Double.MIN_VALUE;
		int su = 0;
		int res = 0;
		for (int col = 0; col < length; col++) {
			for (int s = 0; s < order; s++) {
				if (residueDistances.get(s, col) != -1) {
					if (residueDistances.get(s, col) > maxDist) {
						// geometric distribution
						if (rnd.nextDouble() > 0.5) {
							su = s;
							res = col;
							maxDist = residueDistances.get(s, col);
						}
					}
				}
			}
		}

		// Insert the gap at the position
		Integer residueL = block.get(su).get(res);
		if (residueL != null) {
			freePool.add(residueL);
			Collections.sort(freePool);
		} else
			return false; // If there was a gap already in the position.

		block.get(su).set(res, null);
		checkGaps();
		return true;
	}

	/**
	 * Move all the block residues of one repeat one position to the left or
	 * right and move the corresponding boundary residues from the freePool to
	 * the block, and viceversa.
	 * <p>
	 * The boundaries are determined by any irregularity (either a gap or a
	 * discontinuity in the alignment.
	 */
	private boolean shiftRow() {

		int su = rnd.nextInt(order); // Select the repeat
		int rl = rnd.nextInt(2); // Select between moving right (0) or left (1)
		int res = rnd.nextInt(length); // Residue as a pivot

		// When the pivot residue is null try to add a residue from the freePool
		if (block.get(su).get(res) == null) {

			int right = res;
			int left = res;
			// Find the boundary to the right abd left
			while (block.get(su).get(right) == null && right < length - 1) {
				right++;
			}
			while (block.get(su).get(left) == null && left > 0) {
				left--;
			}

			// If they both are null the whole block is null
			if (block.get(su).get(left) == null
					&& block.get(su).get(right) == null) {
				return false;
			} else if (block.get(su).get(left) == null) {
				// Choose the sequentially previous residue of the known one
				Integer residue = block.get(su).get(right) - 1;
				if (freePool.contains(residue)) {
					block.get(su).set(res, residue);
					freePool.remove(residue);
				} else
					return false;
			} else if (block.get(su).get(right) == null) {
				// Choose the sequentially next residue of the known one
				Integer residue = block.get(su).get(left) + 1;
				if (freePool.contains(residue)) {
					block.get(su).set(res, residue);
					freePool.remove(residue);
				} else
					return false;
			} else {
				// If boundaries are consecutive swap null and position (R or L)
				if (block.get(su).get(right) == block.get(su).get(left) + 1) {
					switch (rl) {
					case 0: // to the right
						block.get(su).set(right - 1, block.get(su).get(right));
						block.get(su).set(right, null);
						break;
					case 1: // to the left
						block.get(su).set(left + 1, block.get(su).get(left));
						block.get(su).set(left, null);
						break;
					}
				} else {
					// Choose randomly a residue in between left and right to
					// add
					Integer residue = rnd.nextInt(block.get(su).get(right)
							- block.get(su).get(left) - 1)
							+ block.get(su).get(left) + 1;

					if (freePool.contains(residue)) {
						block.get(su).set(res, residue);
						freePool.remove(residue);
					}
				}
			}
			return true;
		}

		// When the residue is different than null
		switch (rl) {
		case 0: // Move to the right

			int leftBoundary = res - 1;
			int leftPrevRes = res;
			while (true) {
				if (leftBoundary < 0)
					break;
				else {
					if (block.get(su).get(leftBoundary) == null) {
						break; // gap
					} else if (block.get(su).get(leftPrevRes) > block.get(su)
							.get(leftBoundary) + 1) {
						break; // discontinuity
					}
				}
				leftPrevRes = leftBoundary;
				leftBoundary--;
			}
			leftBoundary++;

			int rightBoundary = res + 1;
			int rightPrevRes = res;
			while (true) {
				if (rightBoundary == length)
					break;
				else {
					if (block.get(su).get(rightBoundary) == null) {
						break; // gap
					} else if (block.get(su).get(rightPrevRes) + 1 < block.get(
							su).get(rightBoundary)) {
						break; // discontinuity
					}
				}
				rightPrevRes = rightBoundary;
				rightBoundary++;
			}
			rightBoundary--;

			// Residues at the boundary
			Integer residueR0 = block.get(su).get(rightBoundary);
			Integer residueL0 = block.get(su).get(leftBoundary);

			// Remove residue at the right of the block and add to the freePool
			block.get(su).remove(rightBoundary);
			if (residueR0 != null) {
				freePool.add(residueR0);
				Collections.sort(freePool);
			}

			// Add the residue at the left of the block
			residueL0 -= 1; // cannot be null, throw exception if it is
			if (freePool.contains(residueL0)) {
				block.get(su).add(leftBoundary, residueL0);
				freePool.remove(residueL0);
			} else {
				block.get(su).add(leftBoundary, null);
			}
			break;

		case 1: // Move to the left

			int leftBoundary1 = res - 1;
			int leftPrevRes1 = res;
			while (true) {
				if (leftBoundary1 < 0)
					break;
				else {
					if (block.get(su).get(leftBoundary1) == null) {
						break; // gap
					} else if (block.get(su).get(leftPrevRes1) > block.get(su)
							.get(leftBoundary1) + 1) {
						break; // discontinuity
					}
				}
				leftPrevRes1 = leftBoundary1;
				leftBoundary1--;
			}
			leftBoundary1++;

			int rightBoundary1 = res + 1;
			int rightPrevRes1 = res;
			while (true) {
				if (rightBoundary1 == length)
					break;
				else {
					if (block.get(su).get(rightBoundary1) == null) {
						break; // gap
					} else if (block.get(su).get(rightPrevRes1) + 1 < block
							.get(su).get(rightBoundary1)) {
						break; // discontinuity
					}
				}
				rightPrevRes1 = rightBoundary1;
				rightBoundary1++;
			}
			rightBoundary1--;

			// Residues at the boundary
			Integer residueR1 = block.get(su).get(rightBoundary1);
			Integer residueL1 = block.get(su).get(leftBoundary1);

			// Add the residue at the right of the block
			residueR1 += 1; // cannot be null
			if (freePool.contains(residueR1)) {
				if (rightBoundary1 == length - 1)
					block.get(su).add(residueR1);
				else
					block.get(su).add(rightBoundary1 + 1, residueR1);
				freePool.remove(residueR1);
			} else {
				block.get(su).add(rightBoundary1 + 1, null);
			}

			// Remove the residue at the left of the block
			block.get(su).remove(leftBoundary1);
			freePool.add(residueL1);
			Collections.sort(freePool);
			break;
		}
		checkGaps();
		return true;
	}

	/**
	 * It extends the alignment one position to the right or to the left of a
	 * randomly selected position by moving the consecutive residues of each
	 * repeat (if present) from the freePool to the block.
	 * <p>
	 * If there are not enough residues in the freePool it introduces gaps.
	 */
	private boolean expandBlock() {

		boolean moved = false;

		int rl = rnd.nextInt(2); // Select between right (0) or left (1)
		int res = rnd.nextInt(length); // Residue as a pivot

		switch (rl) {
		case 0:

			int rightBoundary = res;
			int[] previousPos = new int[order];
			for (int su = 0; su < order; su++)
				previousPos[su] = -1;

			// Search a position to the right that has at minimum Rmin
			while (length - 1 > rightBoundary) {
				int noncontinuous = 0;
				for (int su = 0; su < order; su++) {
					if (block.get(su).get(rightBoundary) == null) {
						continue;
					} else if (previousPos[su] == -1) {
						previousPos[su] = block.get(su).get(rightBoundary);
					} else if (block.get(su).get(rightBoundary) > previousPos[su] + 1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin)
					rightBoundary++;
				else
					break;
			}
			if (rightBoundary > 0)
				rightBoundary--;

			// Expand the block with the residues at the repeat boundaries
			for (int su = 0; su < order; su++) {
				Integer residueR = block.get(su).get(rightBoundary);
				if (residueR == null) {
					if (rightBoundary == length - 1)
						block.get(su).add(null);
					else
						block.get(su).add(rightBoundary + 1, null);
				} else if (freePool.contains(residueR + 1)) {
					Integer residueAdd = residueR + 1;
					if (rightBoundary == length - 1) {
						block.get(su).add(residueAdd);
					} else
						block.get(su).add(rightBoundary + 1, residueAdd);
					freePool.remove(residueAdd);
				} else {
					if (rightBoundary == length - 1)
						block.get(su).add(null);
					else
						block.get(su).add(rightBoundary + 1, null);
				}
			}
			length++;
			moved = true;
			break;

		case 1:

			int leftBoundary = res;
			int[] nextPos = new int[order];
			for (int su = 0; su < order; su++)
				nextPos[su] = -1;

			// Search a position to the right that has at minimum Rmin
			while (leftBoundary > 0) {
				int noncontinuous = 0;
				for (int su = 0; su < order; su++) {
					if (block.get(su).get(leftBoundary) == null) {
						continue;
					} else if (nextPos[su] == -1) {
						nextPos[su] = block.get(su).get(leftBoundary);
					} else if (block.get(su).get(leftBoundary) < nextPos[su] - 1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin)
					leftBoundary--;
				else
					break;
			}

			// Expand the block with the residues at the repeat boundaries
			for (int su = 0; su < order; su++) {
				Integer residueL = block.get(su).get(leftBoundary);
				if (residueL == null) {
					block.get(su).add(leftBoundary, null);
				} else if (freePool.contains(residueL - 1)) {
					Integer residueAdd = residueL - 1;
					block.get(su).add(leftBoundary, residueAdd);
					freePool.remove(residueAdd);
				} else {
					block.get(su).add(leftBoundary, null);
				}
			}
			length++;
			moved = true;
			break;
		}
		if (moved)
			return !checkGaps();
		return moved;
	}

	/**
	 * Deletes an alignment column at a randomly selected position.
	 *
	 * @throws StructureException
	 * @throws RefinerFailedException
	 */
	private boolean shrinkBlock() throws StructureException,
			RefinerFailedException {

		// Let shrink moves only if the repeat is larger enough
		if (repeatCore <= Lmin)
			return false;

		// Select column by maximum distance
		updateMultipleAlignment();
		Matrix residueDistances = MultipleAlignmentTools
				.getAverageResidueDistances(msa);

		double maxDist = Double.MIN_VALUE;
		double[] colDistances = new double[length];
		int res = 0;
		for (int col = 0; col < length; col++) {
			int normalize = 0;
			for (int s = 0; s < order; s++) {
				if (residueDistances.get(s, col) != -1) {
					colDistances[col] += residueDistances.get(s, col);
					normalize++;
				}
			}
			colDistances[col] /= normalize;
			if (colDistances[col] > maxDist) {
				// geometric distribution
				if (rnd.nextDouble() > 0.5) {
					maxDist = colDistances[col];
					res = col;
				}
			}
		}

		for (int su = 0; su < order; su++) {
			Integer residue = block.get(su).get(res);
			block.get(su).remove(res);
			if (residue != null)
				freePool.add(residue);
			Collections.sort(freePool);
		}
		length--;
		checkGaps();
		return true;
	}

	/**
	 * Calculates the probability of accepting a bad move given the iteration
	 * step and the score change.
	 * <p>
	 * Function: p=(C-AS)/(C*sqrt(step)) Added a normalization factor so that
	 * the probability approaches 0 when the maxIter is reached.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {

		double prob = (C + AS) / (C * Math.sqrt(m));
		double norm = (1 - (m * 1.0) / maxIter); // Normalization factor
		return Math.min(Math.max(prob * norm, 0.0), 1.0);
	}

	/**
	 * Save the evolution of the optimization process as a csv file.
	 */
	private void saveHistory(String folder) throws IOException {

		String name = msa.getStructureIdentifier(0).getIdentifier();
		FileWriter writer = new FileWriter(folder + name
				+ "-symm_opt.csv");
		writer.append("Step,Time,RepeatLength,RMSD,TMscore,MCscore\n");

		for (int i = 0; i < lengthHistory.size(); i++) {
			writer.append(i * saveStep + ",");
			writer.append(timeHistory.get(i) + ",");
			writer.append(lengthHistory.get(i) + ",");
			writer.append(rmsdHistory.get(i) + ",");
			writer.append(tmScoreHistory.get(i) + ",");
			writer.append(mcScoreHistory.get(i) + "\n");
		}

		writer.flush();
		writer.close();
	}

}
