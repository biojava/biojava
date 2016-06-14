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

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.MatrixListener;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.secstruc.SecStrucCalc;
import org.biojava.nbio.structure.secstruc.SecStrucTools;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Identify the symmetries in a structure by running an alignment of the
 * structure against itself disabling the diagonal of the identity alignment.
 * <p>
 * Iterating over previous results and disabling the diagonal of the previous
 * alignments can also be done with this implementation, which will generate a
 * set of self-alignments (disabled, because none improvements were shown, but
 * can be turn on manually).
 * <p>
 * Multiple levels of symmetry can be analyzed by finding symmetries in repeats
 * of previous results. This feature allows to find multiple symmetry axes.
 * <p>
 * The alignment is then refined to obtain a consistent alignment among all
 * residues of the structure and organized into different parts, called
 * symmetric repeats.
 * <p>
 * After refinement of the initial alignment, an optimization step can be used
 * to improve the overall score of the repeat multiple alignment.
 *
 * @author Andreas Prlic
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class CeSymm {

	/**
	 * Version History:
	 * <p>
	 * <ul>
	 * <li>1.0 - initial implementation of CE-Symm.
	 * <li>1.1 - enable multiple CE-Symm runs to calculate all self-alignments.
	 * <li>2.0 - refine the alignment for consistency of repeat definition.
	 * <li>2.1 - optimize the alignment to improve the score.
	 * <li>2.2 - run multiple symmetry levels recursively to find PG and
	 * hierarchical symmetries.
	 * </ul>
	 * </li>
	 */
	public static final String version = "2.2";
	public static final String algorithmName = "jCE-symm";
	private static final Logger logger = LoggerFactory.getLogger(CeSymm.class);
	private final static boolean multiPass = false; // multiple self-alignments

	/**
	 * Prevent instantiation. Static class.
	 */
	private CeSymm() {
	}

	private static Matrix align(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
			CESymmParameters params, Matrix origM, CECalculator calculator,
			int counter) throws StructureException {

		int fragmentLength = params.getWinSize();
		Atom[] ca2clone = StructureTools.cloneAtomArray(ca2);

		int rows = ca1.length;
		int cols = ca2.length;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		int blankWindowSize = fragmentLength;
		if (origM == null) {

			// Build alignment ca1 to ca2-ca2
			afpChain = calculator.extractFragments(afpChain, ca1, ca2clone);

			origM = SymmetryTools.blankOutPreviousAlignment(afpChain, ca2,
					rows, cols, calculator, null, blankWindowSize);

		} else {
			// we are doing an iteration on a previous alignment
			// mask the previous alignment
			origM = SymmetryTools.blankOutPreviousAlignment(afpChain, ca2,
					rows, cols, calculator, origM, blankWindowSize);
		}

		Matrix clone = (Matrix) origM.clone();

		// that's the matrix to run the alignment on..
		calculator.setMatMatrix(clone.getArray());

		calculator.traceFragmentMatrix(afpChain, ca1, ca2clone);

		final Matrix origMfinal = (Matrix) origM.clone();
		// Add a matrix listener to keep the blacked zones in max.
		calculator.addMatrixListener(new MatrixListener() {

			@Override
			public double[][] matrixInOptimizer(double[][] max) {

				// Check every entry of origM for blacked out regions
				for (int i = 0; i < max.length; i++) {
					for (int j = 0; j < max[i].length; j++) {
						if (origMfinal.getArray()[i][j] > 1e9) {
							max[i][j] = -origMfinal.getArray()[i][j];
						}
					}
				}
				return max;
			}

			@Override
			public boolean[][] initializeBreakFlag(boolean[][] brkFlag) {

				return brkFlag;
			}
		});

		calculator.nextStep(afpChain, ca1, ca2clone);

		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(version);

		afpChain.setDistanceMatrix(origM);

		return origMfinal;

	}

	@SuppressWarnings("unused")
	protected static CeSymmResult align(Atom[] atoms, CESymmParameters params)
			throws StructureException {

		CeSymmResult result = new CeSymmResult();
		result.setParams(params);
		result.setAtoms(atoms);

		// STEP 1: prepare all the information for the symmetry alignment
		Atom[] ca2 = StructureTools.duplicateCA2(atoms);
		int rows = atoms.length;
		int cols = ca2.length;

		if (rows == 0 || cols == 0) {
			throw new StructureException("Aligning empty structure");
		}

		Matrix origM = null;
		AFPChain myAFP = new AFPChain(algorithmName);
		CECalculator calculator = new CECalculator(params);
		Matrix lastMatrix = null;

		List<AFPChain> selfAlignments = new ArrayList<AFPChain>();
		AFPChain optimalAFP = null;

		// STEP 2: perform the self-alignments of the structure
		int i = 0;
		do {
			if (origM != null)
				myAFP.setDistanceMatrix((Matrix) origM.clone());

			origM = align(myAFP, atoms, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, atoms, ca2);
			myAFP.setTMScore(tmScore2);

			AFPChain newAFP = (AFPChain) myAFP.clone();
			newAFP = CeCPMain.postProcessAlignment(newAFP, atoms, ca2,
					calculator);

			// Calculate and set the TM score for the newAFP alignment
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, atoms, ca2);
			newAFP.setTMScore(tmScore3);

			// Determine if the alignment is significant, stop if false
			if (tmScore3 < params.getUnrefinedScoreThreshold()) {
				// If it is the first alignment save it anyway
				if (i == 0)
					selfAlignments.add(newAFP);
				// store final matrix
				lastMatrix = newAFP.getDistanceMatrix().copy();
				break;
			}

			// If it is a symmetric alignment add it to the allAlignments list
			selfAlignments.add(newAFP);

			i++;

		} while (i < params.getMaxSymmOrder() && multiPass);

		// We reached the maximum order, so blank out the final alignment
		if (lastMatrix == null && selfAlignments.size() > 1 && multiPass) {
			AFPChain last = selfAlignments.get(selfAlignments.size() - 1);
			lastMatrix = SymmetryTools.blankOutPreviousAlignment(last, ca2,
					last.getCa1Length(), last.getCa2Length(), calculator,
					origM, params.getWinSize());
			lastMatrix = lastMatrix.getMatrix(0, last.getCa1Length() - 1, 0,
					last.getCa2Length() - 1);
		}

		// Extract the structure identifier
		optimalAFP = selfAlignments.get(0);
		StructureIdentifier id = atoms[0].getGroup().getChain().getStructure()
				.getStructureIdentifier();
		optimalAFP.setName1(id.getIdentifier());
		optimalAFP.setName2(id.getIdentifier());

		// Store the optimal self-alignment
		result.setSelfAlignment(optimalAFP);
		result.setStructureId(id);

		// Determine the symmetry Type or get the one in params
		SymmetryType type = params.getSymmType();
		if (type == SymmetryType.AUTO) {
			if (result.getSelfAlignment().getBlockNum() == 1) {
				type = SymmetryType.OPEN;
				logger.info("Open Symmetry detected");
			} else {
				type = SymmetryType.CLOSED;
				logger.info("Close Symmetry detected");
			}
		}

		// Do not try the refinement if the self-alignment is not significant
		if (optimalAFP.getTMScore() < params.getUnrefinedScoreThreshold()){
			result.setNumRepeats(1);
			return result;
		}

		// STEP 3: order detection & symmetry refinement, apply consistency
		try {
			// ORDER DETECTION
			OrderDetector orderDetector = null;
			int order = 1;
			switch (params.getOrderDetectorMethod()) {
			case USER_INPUT:
				order = params.getUserOrder();
				break;
			case SEQUENCE_FUNCTION:
				// Does not work for OPEN alignments
				if (type == SymmetryType.CLOSED) {
					orderDetector = new SequenceFunctionOrderDetector(
							params.getMaxSymmOrder(), 0.4f);
					order = orderDetector.calculateOrder(
							result.getSelfAlignment(), atoms);
					break;
				}
			case ANGLE:
				// Does not work for OPEN alignments
				if (type == SymmetryType.CLOSED) {
					orderDetector = new AngleOrderDetectorPlus(
							params.getMaxSymmOrder());
					order = orderDetector.calculateOrder(
							result.getSelfAlignment(), atoms);
					break;
				}
			case GRAPH_COMPONENT:
				orderDetector = new GraphComponentOrderDetector();
				order = orderDetector.calculateOrder(result.getSelfAlignment(),
						atoms);
				break;
			}
			result.setNumRepeats(order);
			
			// REFINEMENT
			SymmetryRefiner refiner = null;
			switch (params.getRefineMethod()) {
			case NOT_REFINED:
				return result;
			case SEQUENCE_FUNCTION:
				// Does not work for OPEN alignments
				if (type == SymmetryType.CLOSED) {
					refiner = new SequenceFunctionRefiner();
					break;
				}
			case GRAPH_COMPONENT:
				refiner = new GraphComponentRefiner();
				break;
			}

			MultipleAlignment msa = refiner.refine(result.getSelfAlignment(),
					atoms, order);

			// Refinement succeeded, store results
			result.setMultipleAlignment(msa);
			result.setNumRepeats(msa.size());
			result.setRefined(true);

		} catch (RefinerFailedException e) {
			logger.info("Refinement failed: " + e.getMessage());
			return result;
		}

		// STEP 4: symmetry axes
		SymmetryAxes axes = new SymmetryAxes();
		int order = result.getMultipleAlignment().size();
		Matrix4d axis = result.getMultipleAlignment().getBlockSet(0)
				.getTransformations().get(1);
		axes.addAxis(axis, order, type);
		
		result.setAxes(axes);
		return result;
	}

	/**
	 * Analyze the symmetries of the input Atom array using the DEFAULT
	 * parameters.
	 *
	 * @param atoms
	 *            representative Atom array of the Structure
	 * @return CeSymmResult
	 * @throws StructureException
	 */
	public static CeSymmResult analyze(Atom[] atoms) throws StructureException {
		CESymmParameters params = new CESymmParameters();
		return analyze(atoms, params);
	}

	/**
	 * Analyze the symmetries of the input Atom array using the provided
	 * parameters.
	 *
	 * @param atoms
	 *            representative Atom array of the Structure
	 * @param param
	 *            CeSymmParameters bean
	 * @return CeSymmResult
	 * @throws StructureException
	 */
	public static CeSymmResult analyze(Atom[] atoms, CESymmParameters params)
			throws StructureException {

		if (atoms.length < 1)
			throw new IllegalArgumentException("Empty Atom array given.");

		// If the SSE information is needed, we calculate it if the user did not
		if (params.getSSEThreshold() > 0) {
			Structure s = atoms[0].getGroup().getChain().getStructure();
			if (SecStrucTools.getSecStrucInfo(s).isEmpty()) {
				logger.info("Calculating Secondary Structure...");
				SecStrucCalc ssp = new SecStrucCalc();
				ssp.calculate(s, true);
			}
		}

		CeSymmIterative iter = new CeSymmIterative(params);
		CeSymmResult result = iter.execute(atoms);

		if (result.isRefined()) {
			// Optimize the global alignment freely once more (final step)
			if (params.getOptimization() && result.getSymmLevels() > 1) {
				try {
					SymmOptimizer optimizer = new SymmOptimizer(result);
					MultipleAlignment optimized = optimizer.optimize();
					// Set the optimized MultipleAlignment and the axes
					result.setMultipleAlignment(optimized);
				} catch (RefinerFailedException e) {
					logger.info("Final optimization failed:" + e.getMessage());
				}
			}
			result.getMultipleAlignment().getEnsemble()
					.setStructureIdentifiers(result.getRepeatsID());
		}
		return result;
	}

	/**
	 * Analyze a single level of symmetry.
	 *
	 * @param atoms
	 *            Atom array of the current level
	 * @return CeSymmResult
	 * @throws StructureException
	 */
	public static CeSymmResult analyzeLevel(Atom[] atoms,
			CESymmParameters params) throws StructureException {

		if (atoms.length < 1)
			throw new IllegalArgumentException("Empty Atom array given.");

		CeSymmResult result = align(atoms, params);

		if (result.isRefined()) {
			// STEP 5: symmetry alignment optimization
			if (result.getParams().getOptimization()) {
				try {
					MultipleAlignment msa = result.getMultipleAlignment();
					SymmOptimizer optimizer = new SymmOptimizer(result);
					msa = optimizer.optimize();
					result.setMultipleAlignment(msa);
				} catch (RefinerFailedException e) {
					logger.debug("Optimization failed:" + e.getMessage());
				}
			}
		}
		return result;
	}

}
