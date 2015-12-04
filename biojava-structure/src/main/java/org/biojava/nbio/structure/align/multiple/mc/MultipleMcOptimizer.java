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
package org.biojava.nbio.structure.align.multiple.mc;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.Callable;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.util.CoreSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.align.multiple.util.MultipleSuperimposer;
import org.biojava.nbio.structure.jama.Matrix;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class takes a MultipleAlignment seed previously generated and runs a
 * Monte Carlo optimization in order to improve the overall score and highlight
 * common structural motifs.
 * <p>
 * The seed alignment can be flexible, non-topological or include CP, but this
 * optimization will not change the number of flexible parts {@link BlockSet}s
 * or non-topological regions {@link Block}. Thus, the definition of those parts
 * depend exclusively on the pairwise alignment (or user alignment) used to
 * generate the seed multiple alignment.
 * <p>
 * This class implements Callable, because multiple instances of the
 * optimization can be run in parallel.
 * 
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleMcOptimizer implements Callable<MultipleAlignment> {

	private static final Logger logger = LoggerFactory
			.getLogger(MultipleMcOptimizer.class);

	private Random rnd;
	private MultipleSuperimposer imposer;

	// Optimization parameters
	private int Rmin; // number of aligned structures without a gap
	private int Lmin; // Minimum alignment length of a Block
	private int convergenceSteps; // Steps without score change before stopping
	private double C; // Probability function constant

	// Score function parameters - they are defined by the user
	private double Gopen; // Penalty for opening gap
	private double Gextend; // Penalty for extending gaps
	private double dCutoff; // max allowed residue distance

	// Alignment Information
	private MultipleAlignment msa; // Alignment to optimize
	private List<SortedSet<Integer>> freePool; // unaligned residues
	private List<Atom[]> atomArrays;

	// Alignment Properties
	private int size; // number of structures in the alignment
	private int blockNr; // the number of Blocks in the alignment
	private double mcScore; // Optimization score, objective function

	// Variables that store the history of the optimization - slower if on
	private static final boolean history = false;
	private static final String pathToHistory = "McOptHistory.csv";
	private List<Integer> lengthHistory;
	private List<Double> rmsdHistory;
	private List<Double> scoreHistory;

	/**
	 * Constructor. Sets the internal variables from the parameters. To run the
	 * optimization use the call (in a different thread) or optimize methods.
	 * 
	 * @param seedAln
	 *            MultipleAlignment to be optimized.
	 * @param params
	 *            the parameter beam
	 * @param reference
	 *            the index of the most similar structure to all others
	 * @throws StructureException
	 */
	public MultipleMcOptimizer(MultipleAlignment seedAln,
			MultipleMcParameters params, int reference) {

		MultipleAlignmentEnsemble e = seedAln.getEnsemble().clone();
		msa = e.getMultipleAlignment(0);
		atomArrays = msa.getAtomArrays();
		size = seedAln.size();

		rnd = new Random(params.getRandomSeed());
		Gopen = params.getGapOpen();
		Gextend = params.getGapExtension();
		dCutoff = params.getDistanceCutoff();
		imposer = new CoreSuperimposer(reference);

		if (params.getConvergenceSteps() == 0) {
			List<Integer> lens = new ArrayList<Integer>();
			for (int i = 0; i < size; i++)
				lens.add(atomArrays.get(i).length);
			convergenceSteps = Collections.min(lens) * size;
		} else
			convergenceSteps = params.getConvergenceSteps();

		if (params.getMinAlignedStructures() == 0) {
			Rmin = Math.max(size / 3, 2); // 33% of the structures aligned
		} else {
			Rmin = Math
					.min(Math.max(params.getMinAlignedStructures(), 2), size);
		}
		C = 20 * size;
		Lmin = params.getMinBlockLen();

		// Delete all shorter than Lmin blocks, and empty blocksets
		List<Block> toDelete = new ArrayList<Block>();
		List<BlockSet> emptyBs = new ArrayList<BlockSet>();

		for (Block b : msa.getBlocks()) {
			if (b.getCoreLength() < Lmin) {
				toDelete.add(b);
				logger.warn("Deleting a Block: coreLength < Lmin.");
			}
		}
		for (Block b : toDelete) {
			for (BlockSet bs : msa.getBlockSets()) {
				bs.getBlocks().remove(b);
				if (bs.getBlocks().size() == 0)
					emptyBs.add(bs);
			}
		}
		for (BlockSet bs : emptyBs) {
			msa.getBlockSets().remove(bs);
		}

		blockNr = msa.getBlocks().size();
		if (blockNr < 1) {
			throw new IllegalArgumentException(
					"Optimization: empty seed alignment, no Blocks found.");
		}
	}

	@Override
	public MultipleAlignment call() throws Exception {
		return optimize();
	}

	/**
	 * Initialize the freePool and all the variables needed for the
	 * optimization.
	 * 
	 * @throws StructureException
	 */
	private void initialize() throws StructureException {

		// Initialize alignment variables
		freePool = new ArrayList<SortedSet<Integer>>();
		List<List<Integer>> aligned = new ArrayList<List<Integer>>();

		// Generate freePool residues from the ones not aligned
		for (int i = 0; i < size; i++) {
			List<Integer> residues = new ArrayList<Integer>();
			for (BlockSet bs : msa.getBlockSets()) {
				for (Block b : bs.getBlocks()) {
					for (int l = 0; l < b.length(); l++) {
						Integer residue = b.getAlignRes().get(i).get(l);
						if (residue != null)
							residues.add(residue);
					}
				}
			}
			aligned.add(residues);
			freePool.add(new TreeSet<Integer>());
		}

		// Add any residue not aligned to the free pool for every structure
		for (int i = 0; i < size; i++) {
			for (int k = 0; k < atomArrays.get(i).length; k++) {
				if (!aligned.get(i).contains(k))
					freePool.get(i).add(k);
			}
		}

		// Set the superposition and score for the seed aligment
		checkGaps();
		msa.clear();
		imposer.superimpose(msa);
		mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend,
				dCutoff);

		// Initialize the history variables
		if (history) {
			lengthHistory = new ArrayList<Integer>();
			rmsdHistory = new ArrayList<Double>();
			scoreHistory = new ArrayList<Double>();
		}
	}

	/**
	 * Optimization method based in a Monte-Carlo approach. Starting from the
	 * refined alignment uses 4 types of moves:
	 * <p>
	 * <ul>
	 * <li>Shift Row: if there are enough freePool residues available.
	 * <li>Expand Block: add another alignment column.
	 * <li>Shrink Block: move a block column to the freePool.
	 * <li>Insert gap: insert a gap in a random position of the alignment.
	 * </ul>
	 * </li>
	 */
	public MultipleAlignment optimize() throws StructureException {

		initialize();

		int conv = 0; // Number of steps without an alignment improvement
		int i = 1;
		int maxIter = convergenceSteps * 100;

		while (i < maxIter && conv < convergenceSteps) {

			// Save the state of the system
			MultipleAlignment lastMSA = msa.clone();
			List<SortedSet<Integer>> lastFreePool = new ArrayList<SortedSet<Integer>>();
			for (int k = 0; k < size; k++) {
				SortedSet<Integer> p = new TreeSet<Integer>();
				for (Integer l : freePool.get(k))
					p.add(l);
				lastFreePool.add(p);
			}
			double lastScore = mcScore;

			boolean moved = false;

			while (!moved) {
				// Randomly select one of the steps to modify the alignment
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

			// Get the score of the new alignment
			msa.clear();
			imposer.superimpose(msa);
			mcScore = MultipleAlignmentScorer.getMCScore(msa, Gopen, Gextend,
					dCutoff);

			double AS = mcScore - lastScore;
			double prob = 1.0;

			if (AS < 0) {

				// Probability of accepting the move
				prob = probabilityFunction(AS, i, maxIter);
				double p = rnd.nextDouble();
				// Reject the move
				if (p > prob) {
					msa = lastMSA;
					freePool = lastFreePool;
					mcScore = lastScore;
					conv++;

				} else
					conv = 0;

			} else
				conv = 0;

			logger.debug("Step: " + i + ": --prob: " + prob
					+ ", --score change: " + AS + ", --conv: " + conv);

			if (history) {
				if (i % 100 == 1) {
					lengthHistory.add(msa.length());
					rmsdHistory.add(MultipleAlignmentScorer.getRMSD(msa));
					scoreHistory.add(mcScore);
				}
			}

			i++;
		}

		// Return Multiple Alignment
		imposer.superimpose(msa);
		MultipleAlignmentScorer.calculateScores(msa);
		msa.putScore(MultipleAlignmentScorer.MC_SCORE, mcScore);

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
	 * Method that loops through all the alignment columns and checks that there
	 * are no more gaps than the maximum allowed, Rmin.
	 * <p>
	 * There must be at least Rmin residues different than null in every
	 * alignment column. In case there is a column with more gaps it will be
	 * shrinked (moved to freePool).
	 * 
	 * @return true if any columns has been shrinked and false otherwise
	 */
	private boolean checkGaps() {

		boolean shrinkedAny = false;

		List<List<Integer>> shrinkColumns = new ArrayList<List<Integer>>();
		// Loop for each Block
		for (Block b : msa.getBlocks()) {
			List<Integer> shrinkCol = new ArrayList<Integer>();
			// Loop for each column in the Block
			for (int res = 0; res < b.length(); res++) {
				int gapCount = 0;
				// count the gaps in the column
				for (int su = 0; su < size; su++) {
					if (b.getAlignRes().get(su).get(res) == null)
						gapCount++;
				}
				if ((size - gapCount) < Rmin) {
					// Add the column to the shrink list
					shrinkCol.add(res);
				}
			}
			shrinkColumns.add(shrinkCol);
		}
		// Shrink columns that have more gaps than allowed
		for (int b = 0; b < blockNr; b++) {
			for (int col = shrinkColumns.get(b).size() - 1; col >= 0; col--) {
				for (int str = 0; str < size; str++) {
					Block bk = msa.getBlock(b);
					Integer residue = bk.getAlignRes().get(str)
							.get(shrinkColumns.get(b).get(col));
					bk.getAlignRes().get(str)
							.remove((int) shrinkColumns.get(b).get(col));
					if (residue != null) {
						freePool.get(str).add(residue);
					}
				}
				shrinkedAny = true;
			}
		}
		return shrinkedAny;
	}

	/**
	 * Insert a gap in one of the structures in a random position of the
	 * alignment.
	 * <p>
	 * The distribution is not uniform, because positions with higher average
	 * distance to their aligned neighbors are more likely to be gapped.
	 * <p>
	 * A gap is a null in the Block position.
	 * 
	 * @return true if the alignment has been changed, false otherwise.
	 */
	private boolean insertGap() {

		// Select residue by maximum distance
		Matrix residueDistances = MultipleAlignmentTools
				.getAverageResidueDistances(msa);
		double maxDist = Double.MIN_VALUE;
		int structure = 0;
		int block = 0;
		int position = 0;
		int column = 0;
		for (int b = 0; b < blockNr; b++) {
			for (int col = 0; col < msa.getBlock(b).length(); col++) {
				for (int str = 0; str < size; str++) {
					if (residueDistances.get(str, column) != -1) {
						if (residueDistances.get(str, column) > maxDist) {
							// Geometric distribution
							if (rnd.nextDouble() > 0.5) {
								structure = str;
								block = b;
								position = col;
								maxDist = residueDistances.get(str, column);
							}
						}
					}
				}
				column++;
			}
		}
		Block bk = msa.getBlock(block);
		if (bk.getCoreLength() <= Lmin)
			return false;

		// Insert the gap at the position
		Integer residueL = bk.getAlignRes().get(structure).get(position);
		if (residueL != null) {
			freePool.get(structure).add(residueL);
		} else
			return false; // If there was a gap already in the position.

		bk.getAlignRes().get(structure).set(position, null);
		checkGaps();
		return true;
	}

	/**
	 * Move all the block residues of one subunit one position to the left or to
	 * the right and move the corresponding boundary residues from the freePool
	 * to the block.
	 * <p>
	 * The boundaries are determined by any irregularity (either a null or a
	 * discontinuity in the alignment).
	 * 
	 * @return true if the alignment has been changed, false otherwise.
	 */
	private boolean shiftRow() {

		int str = rnd.nextInt(size); // Select randomly the subunit
		int rl = rnd.nextInt(2); // Select between moving right (0) or left (1)
		int bk = rnd.nextInt(blockNr); // Select randomly the Block
		int res = rnd.nextInt(msa.getBlock(bk).length());

		Block block = msa.getBlock(bk);
		if (block.getCoreLength() <= Lmin)
			return false;

		// When the pivot residue is null try to add a residue from the freePool
		if (block.getAlignRes().get(str).get(res) == null) {
			// Residues not null at the right and left of the pivot null residue
			int rightRes = res;
			int leftRes = res;
			// Find the boundary to the right abd left
			while (block.getAlignRes().get(str).get(rightRes) == null
					&& rightRes < block.length() - 1) {
				rightRes++;
			}
			while (block.getAlignRes().get(str).get(leftRes) == null
					&& leftRes > 0) {
				leftRes--;
			}

			// If both are null return because the block is empty
			if (block.getAlignRes().get(str).get(leftRes) == null
					&& block.getAlignRes().get(str).get(rightRes) == null) {
				return false;
			} else if (block.getAlignRes().get(str).get(leftRes) == null) {
				// Choose the sequentially previous residue of the known one
				Integer residue = block.getAlignRes().get(str).get(rightRes) - 1;
				if (freePool.get(str).contains(residue)) {
					block.getAlignRes().get(str).set(res, residue);
					freePool.get(str).remove(residue);
				} else
					return false;
			} else if (block.getAlignRes().get(str).get(rightRes) == null) {
				// Choose the sequentially next residue of the known one
				Integer residue = block.getAlignRes().get(str).get(leftRes) + 1;
				if (freePool.contains(residue)) {
					block.getAlignRes().get(str).set(res, residue);
					freePool.get(str).remove(residue);
				} else
					return false;
			} else {
				// If boundaries are consecutive no residue can be added
				if (block.getAlignRes().get(str).get(rightRes) == block
						.getAlignRes().get(str).get(leftRes) + 1) {
					return false;
				} else {
					// Choose randomly a residue in between left and right
					Integer residue = rnd.nextInt(block.getAlignRes().get(str)
							.get(rightRes)
							- block.getAlignRes().get(str).get(leftRes) - 1)
							+ block.getAlignRes().get(str).get(leftRes) + 1;

					if (freePool.get(str).contains(residue)) {
						block.getAlignRes().get(str).set(res, residue);
						freePool.get(str).remove(residue);
					}
				}
			}
			return true;
		}

		// When residue different than null shift the whole block
		switch (rl) {
		case 0: // Move to the right

			// Find the nearest boundary to the left of the pivot
			int leftBoundary = res - 1;
			int leftPrevRes = res;
			while (true) {
				if (leftBoundary < 0)
					break;
				else {
					if (block.getAlignRes().get(str).get(leftBoundary) == null)
						break; // Break if there is a gap (this is the boundary)
					else if (block.getAlignRes().get(str).get(leftPrevRes) > block
							.getAlignRes().get(str).get(leftBoundary) + 1)
						break; // Break if there is a discontinuity
				}
				leftPrevRes = leftBoundary;
				leftBoundary--;
			}
			leftBoundary++;

			// Find the nearest boundary to the right of the pivot
			int rightBoundary = res + 1;
			int rightPrevRes = res;
			while (true) {
				if (rightBoundary == block.length())
					break;
				else {
					if (block.getAlignRes().get(str).get(rightBoundary) == null)
						break; // Break if there is a gap
					else if (block.getAlignRes().get(str).get(rightPrevRes) + 1 < block
							.getAlignRes().get(str).get(rightBoundary))
						break; // Discontinuity
				}
				rightPrevRes = rightBoundary;
				rightBoundary++;
			}
			rightBoundary--;

			// Residues at the boundary
			Integer resR0 = block.getAlignRes().get(str).get(rightBoundary);
			Integer resL0 = block.getAlignRes().get(str).get(leftBoundary);

			// Remove the residue at the right of the block
			block.getAlignRes().get(str).remove(rightBoundary);
			if (resR0 != null)
				freePool.get(str).add(resR0);

			// Add the residue at the left of the block
			if (resL0 != null)
				resL0 -= 1;

			if (freePool.get(str).contains(resL0)) {
				block.getAlignRes().get(str).add(leftBoundary, resL0);
				freePool.get(str).remove(resL0);
			} else
				block.getAlignRes().get(str).add(leftBoundary, null);

			break;

		case 1: // Move to the left

			// Find the nearest boundary to the left of the pivot
			int leftBoundary1 = res - 1;
			int leftPrevRes1 = res;
			while (true) {
				if (leftBoundary1 < 0)
					break;
				else {
					if (block.getAlignRes().get(str).get(leftBoundary1) == null)
						break; // Break if there is a gap (this is the boundary)
					else if (block.getAlignRes().get(str).get(leftPrevRes1) > block
							.getAlignRes().get(str).get(leftBoundary1) + 1)
						break; // Break if there is a discontinuity
				}
				leftPrevRes1 = leftBoundary1;
				leftBoundary1--;
			}
			leftBoundary1++;

			// Find the nearest boundary to the right of the pivot
			int rightBoundary1 = res + 1;
			int rightPrevRes1 = res;
			while (true) {
				if (rightBoundary1 == block.length())
					break;
				else {
					if (block.getAlignRes().get(str).get(rightBoundary1) == null)
						break; // Break if there is a gap
					else if (block.getAlignRes().get(str).get(rightPrevRes1) + 1 < block
							.getAlignRes().get(str).get(rightBoundary1))
						break; // Discontinuity
				}
				rightPrevRes1 = rightBoundary1;
				rightBoundary1++;
			}
			rightBoundary1--;

			// Residues at the boundary
			Integer resR1 = block.getAlignRes().get(str).get(rightBoundary1);
			Integer resL1 = block.getAlignRes().get(str).get(leftBoundary1);

			// Add the residue at the right of the block
			if (resR1 != null)
				resR1 += 1;

			if (freePool.contains(resR1)) {
				if (rightBoundary1 == block.length() - 1) {
					block.getAlignRes().get(str).add(resR1);
				} else
					block.getAlignRes().get(str).add(rightBoundary1 + 1, resR1);

				freePool.get(str).remove(resR1);
			} else
				block.getAlignRes().get(str).add(rightBoundary1 + 1, null);

			// Remove the residue at the left of the block
			block.getAlignRes().get(str).remove(leftBoundary1);
			if (resL1 != null)
				freePool.get(str).add(resL1);

			break;
		}
		checkGaps();
		return true;
	}

	/**
	 * It extends the alignment one position to the right or to the left of a
	 * randomly selected position by moving the consecutive residues of each
	 * subunit (if enough) from the freePool to the block.
	 * <p>
	 * If there are not enough residues in the freePool it introduces gaps.
	 */
	private boolean expandBlock() {

		int rl = rnd.nextInt(2); // Select expanding right (0) or left (1)
		int bk = rnd.nextInt(blockNr); // Select randomly the Block
		int res = rnd.nextInt(msa.getBlock(bk).length());

		Block block = msa.getBlock(bk);
		int gaps = 0; // store the number of gaps in the expansion

		switch (rl) {
		case 0:

			int rightBound = res;
			int[] previousPos = new int[size];
			for (int str = 0; str < size; str++)
				previousPos[str] = -1;

			// Search t the right for >Rmin non consecutive residues
			while (block.length() - 1 > rightBound) {
				int noncontinuous = 0;
				for (int str = 0; str < size; str++) {
					if (block.getAlignRes().get(str).get(rightBound) == null) {
						continue;
					} else if (previousPos[str] == -1) {
						previousPos[str] = block.getAlignRes().get(str)
								.get(rightBound);
					} else if (block.getAlignRes().get(str).get(rightBound) > previousPos[str] + 1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin)
					rightBound++;
				else
					break;
			}
			if (rightBound > 0)
				rightBound--;

			// Expand the block with the residues at the subunit boundaries
			for (int str = 0; str < size; str++) {
				Integer residueR = block.getAlignRes().get(str).get(rightBound);
				if (residueR == null) {
					if (rightBound == block.length() - 1) {
						block.getAlignRes().get(str).add(null);
					} else
						block.getAlignRes().get(str).add(rightBound + 1, null);
					gaps++;
				} else if (freePool.get(str).contains(residueR + 1)) {
					Integer residueAdd = residueR + 1;
					if (rightBound == block.length() - 1) {
						block.getAlignRes().get(str).add(residueAdd);
					} else {
						block.getAlignRes().get(str)
								.add(rightBound + 1, residueAdd);
					}
					freePool.get(str).remove(residueAdd);
				} else {
					if (rightBound == block.length() - 1)
						block.getAlignRes().get(str).add(null);
					else
						block.getAlignRes().get(str).add(rightBound + 1, null);
					gaps++;
				}
			}
			break;

		case 1:

			int leftBoundary = res;
			int[] nextPos = new int[size];
			for (int str = 0; str < size; str++)
				nextPos[str] = -1;

			// Search position to the right with >Rmin non consecutive residues
			while (leftBoundary > 0) {
				int noncontinuous = 0;
				for (int str = 0; str < size; str++) {
					if (block.getAlignRes().get(str).get(leftBoundary) == null)
						continue;
					else if (nextPos[str] == -1) {
						nextPos[str] = block.getAlignRes().get(str)
								.get(leftBoundary);
					} else if (block.getAlignRes().get(str).get(leftBoundary) < nextPos[str] - 1) {
						noncontinuous++;
					}
				}
				if (noncontinuous < Rmin)
					leftBoundary--;
				else
					break;
			}

			// Expand the block with the residues at the subunit boundaries
			for (int str = 0; str < size; str++) {
				Integer residueL = block.getAlignRes().get(str)
						.get(leftBoundary);
				if (residueL == null) {
					block.getAlignRes().get(str).add(leftBoundary, null);
					gaps++;
				} else if (freePool.get(str).contains(residueL - 1)) {
					Integer residueAdd = residueL - 1;
					block.getAlignRes().get(str).add(leftBoundary, residueAdd);
					freePool.get(str).remove(residueAdd);
				} else {
					block.getAlignRes().get(str).add(leftBoundary, null);
					gaps++;
				}
			}
			break;
		}
		if (size - gaps >= Rmin)
			return true;
		else
			checkGaps();
		return false;
	}

	/**
	 * Deletes an alignment column at a randomly selected position.
	 */
	private boolean shrinkBlock() {

		// Select column by maximum distance
		Matrix residueDistances = MultipleAlignmentTools
				.getAverageResidueDistances(msa);
		double[] colDistances = new double[msa.length()];
		double maxDist = Double.MIN_VALUE;
		int position = 0;
		int block = 0;
		int column = 0;
		for (int b = 0; b < msa.getBlocks().size(); b++) {
			for (int col = 0; col < msa.getBlock(b).length(); col++) {
				int normalize = 0;
				for (int s = 0; s < size; s++) {
					if (residueDistances.get(s, column) != -1) {
						colDistances[column] += residueDistances.get(s, column);
						normalize++;
					}
				}
				colDistances[column] /= normalize;
				if (colDistances[column] > maxDist) {
					if (rnd.nextDouble() > 0.5) {
						maxDist = colDistances[column];
						position = col;
						block = b;
					}
				}
				column++;
			}
		}
		Block currentBlock = msa.getBlock(block);
		if (currentBlock.getCoreLength() <= Lmin)
			return false;

		for (int str = 0; str < size; str++) {
			Integer residue = currentBlock.getAlignRes().get(str).get(position);

			currentBlock.getAlignRes().get(str).remove(position);
			if (residue != null)
				freePool.get(str).add(residue);
		}
		return true;
	}

	/**
	 * Calculates the probability of accepting a bad move given the iteration
	 * step and the score change.
	 * <p>
	 * Function: p=(C-AS)/(m*C) , slightly different from the CEMC algorithm.
	 * <p>
	 * Added a normalization factor so that the probability approaches 0 as the
	 * final of the optimization gets closer.
	 */
	private double probabilityFunction(double AS, int m, int maxIter) {

		double prob = (C + AS) / (m * C);
		double norm = (1 - (m * 1.0) / maxIter); // Normalization factor
		return Math.min(Math.max(prob * norm, 0.0), 1.0);
	}

	/**
	 * Save the evolution of the optimization process as a csv file.
	 */
	private void saveHistory(String filePath) throws IOException {

		FileWriter writer = new FileWriter(filePath);
		writer.append("Step,Length,RMSD,Score\n");

		for (int i = 0; i < lengthHistory.size(); i++) {
			writer.append("" + (i * 100));
			writer.append("," + lengthHistory.get(i));
			writer.append("," + rmsdHistory.get(i));
			writer.append("," + scoreHistory.get(i) + "\n");
		}
		writer.flush();
		writer.close();
	}
}
