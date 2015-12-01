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
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.MatrixListener;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.util.CoreSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.secstruc.SecStrucPred;
import org.biojava.nbio.structure.secstruc.SecStrucTools;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Identify the symmetries in a structure by running an alignment of the
 * structure against itself disabling the diagonal of the identity alignment.
 * <p>
 * Iterating recursively over all results and disabling the diagonal of each
 * previous result can also be done with this implementation, which will
 * generate a set of self-alignments.
 * <p>
 * The alignment is then refined to obtain a consistent alignment among all
 * residues of the structure and organized into different parts, called
 * symmetric subunits. After refinement of the initial alignment, an
 * optimization step can be used to improve the overall score of the subunit
 * multiple alignment.
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
	 * <li>2.0 - refine the alignment for consistency of subunit definition.
	 * <li>2.1 - optimize the alignment to improve the score.
	 * <li>2.2 - run multiple symmetry levels recursively to find PG and
	 * hierarchical symmetries.
	 * </ul>
	 * </li>
	 */
	public static final String version = "2.2";
	public static final String algorithmName = "jCE-symm";
	private static final Logger logger = LoggerFactory.getLogger(CeSymm.class);

	private MultipleAlignment msa;
	private List<AFPChain> selfAlignments;
	private SymmetryAxes axes;
	private boolean refined;

	private CESymmParameters params = new CESymmParameters();

	public CeSymm() {
		reset();
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

	protected AFPChain align(Atom[] ca1) throws StructureException {

		// STEP 1: prepare all the information for the symmetry alignment
		Atom[] ca2 = StructureTools.duplicateCA2(ca1);
		int rows = ca1.length;
		int cols = ca2.length;

		if (rows == 0 || cols == 0) {
			throw new StructureException("Aligning empty structure");
		}

		Matrix origM = null;
		AFPChain myAFP = new AFPChain();
		CECalculator calculator = new CECalculator(params);

		// Set multiple to true if multiple alignments are needed for refinement
		boolean multiple = (params.getRefineMethod() == RefineMethod.MULTIPLE);
		Matrix lastMatrix = null;

		// STEP 2: perform the self-alignments of the structure with CECP
		int i = 0;
		do {

			if (origM != null) {
				myAFP.setDistanceMatrix((Matrix) origM.clone());
			}
			origM = align(myAFP, ca1, ca2, params, origM, calculator, i);

			double tmScore2 = AFPChainScorer.getTMScore(myAFP, ca1, ca2);
			myAFP.setTMScore(tmScore2);

			// Clone the AFPChain
			AFPChain newAFP = (AFPChain) myAFP.clone();

			// Post process the alignment
			newAFP = CeCPMain
					.postProcessAlignment(newAFP, ca1, ca2, calculator);

			// Calculate and set the TM score for the newAFP alignment
			double tmScore3 = AFPChainScorer.getTMScore(newAFP, ca1, ca2);
			newAFP.setTMScore(tmScore3);
			logger.debug("Alignment " + (i + 1) + " score: "
					+ newAFP.getTMScore());
			// Determine if the alignment is significant, stop if true
			if (tmScore3 < params.getScoreThreshold()) {
				logger.debug("Not symmetric alignment with TM score: "
						+ newAFP.getTMScore());
				// If it is the first alignment save it anyway
				if (i == 0)
					selfAlignments.add(newAFP);
				// store final matrix &
				lastMatrix = newAFP.getDistanceMatrix().copy();
				break;
			}
			// If it is a symmetric alignment add it to the allAlignments list
			selfAlignments.add(newAFP);

			i++;

		} while (i < params.getMaxSymmOrder() && multiple);

		if (lastMatrix == null && selfAlignments.size() > 1) {
			// We reached the maximum order, so blank out the final alignment
			AFPChain last = selfAlignments.get(selfAlignments.size() - 1);
			lastMatrix = SymmetryTools.blankOutPreviousAlignment(last, ca2,
					last.getCa1Length(), last.getCa2Length(), calculator,
					origM, params.getWinSize());
			lastMatrix = lastMatrix.getMatrix(0, last.getCa1Length() - 1, 0,
					last.getCa2Length() - 1);
		}

		// Extract the optimal alignment
		AFPChain optimalAFP = selfAlignments.get(0);
		if (ca1.length != 0 && ca1[0].getGroup().getChain() != null
				&& ca1[0].getGroup().getChain().getStructure() != null) {
			String name = ca1[0].getGroup().getChain().getStructure().getName();
			optimalAFP.setName1(name);
			optimalAFP.setName2(name);
		}

		if (params.getRefineMethod() == RefineMethod.NOT_REFINED) {
			return optimalAFP;
		}

		// Determine the symmetry Type or get the one in params
		SymmetryType type = params.getSymmType();
		if (type == SymmetryType.AUTO) {
			if (optimalAFP.getBlockNum() == 1) {
				type = SymmetryType.OPEN;
				logger.info("Open Symmetry detected");
			} else {
				type = SymmetryType.CLOSE;
				logger.info("Close Symmetry detected");
			}
		}

		// STEP 3: symmetry refinement, apply consistency in the subunit
		// residues
		Refiner refiner = null;
		int order = 1;
		try {
			switch (type) {
			case CLOSE:
				OrderDetector orderDetector = null;
				switch (params.getOrderDetectorMethod()) {
				case SEQUENCE_FUNCTION:
					orderDetector = new SequenceFunctionOrderDetector(
							params.getMaxSymmOrder(), 0.4f);
					order = orderDetector.calculateOrder(optimalAFP, ca1);
					break;
				case USER_INPUT:
					order = params.getUserOrder();
					break;
				}
				refiner = new SingleRefiner();
				break;
			default: // case OPEN
				refiner = new OpenRefiner();
				order = params.getUserOrder();
				break;
			}

			optimalAFP = refiner.refine(selfAlignments, ca1, order);
			refined = true;

		} catch (RefinerFailedException e) {
			logger.info("Refinement failed: " + e.getMessage());
			return optimalAFP;
		}

		// STEP4: determine the symmetry axis and its subunit dependencies
		order = optimalAFP.getBlockNum();
		Matrix rot = optimalAFP.getBlockRotationMatrix()[0];
		Atom shift = optimalAFP.getBlockShiftVector()[0];
		Matrix4d axis = Calc.getTransformation(rot, shift);

		List<List<Integer>> superposition = new ArrayList<List<Integer>>();
		List<Integer> chain1 = new ArrayList<Integer>();
		List<Integer> chain2 = new ArrayList<Integer>();
		superposition.add(chain1);
		superposition.add(chain2);
		List<Integer> subunitTrans = new ArrayList<Integer>();

		switch (type) {
		case CLOSE:

			for (int bk = 0; bk < order; bk++) {
				chain1.add(bk);
				chain2.add((bk + 1) % order);
				subunitTrans.add(bk);
			}
			axes.addAxis(axis, superposition, subunitTrans, order);
			break;

		default: // case OPEN:

			subunitTrans.add(0);
			for (int bk = 0; bk < order - 1; bk++) {
				chain1.add(bk);
				chain2.add(bk + 1);
				subunitTrans.add(bk + 1);
			}
			axes.addAxis(axis, superposition, subunitTrans, order);
			break;
		}

		return optimalAFP;
	}

	public CESymmParameters getParameters() {
		return params;
	}

	public void setParameters(CESymmParameters parameters) {
		params = parameters;
	}

	public String getAlgorithmName() {
		return algorithmName;
	}

	public String getVersion() {
		return version;
	}

	public boolean isSignificant() throws StructureException {
		double symmetryThreshold = params.getScoreThreshold();
		return SymmetryTools.isSignificant(msa, symmetryThreshold);
	}

	/**
	 * Get the list of all self-alignments.
	 * 
	 * @return List of AFPChain self-alignments
	 */
	public List<AFPChain> getSelfAlignments() {
		return selfAlignments;
	}

	/**
	 * Analyze the symmetries of the input Atom array using the DEFAULT or
	 * previously set parameters.
	 * 
	 * @param atoms
	 *            representative Atom array of the Structure
	 * @return
	 * @throws StructureException
	 */
	public MultipleAlignment analyze(Atom[] atoms) throws StructureException {

		if (params == null)
			params = new CESymmParameters();
		return analyze(atoms, params);
	}

	/**
	 * Analyze the symmetries of the input Atom array using the provided
	 * parameters.
	 * 
	 * @param atoms
	 *            representative Atom array of the Structure
	 * @param param
	 *            CeSymm Parameter bean
	 * @return
	 * @throws StructureException
	 */
	public MultipleAlignment analyze(Atom[] atoms, CESymmParameters param)
			throws StructureException {

		reset();
		if (atoms.length < 1) {
			throw new IllegalArgumentException("Empty Atom array given.");
		}
		this.params = param;
		
		// If the SSE information is needed, we calculate it if the user did not
		if (params.getSSEThreshold() > 0) {
			Structure s = atoms[0].getGroup().getChain().getStructure();
			if (SecStrucTools.getSecStrucInfo(s).isEmpty()){
				SecStrucPred ssp = new SecStrucPred();
				ssp.predict(s, true);
			}
		}
		
		CeSymmIterative iter = new CeSymmIterative(param);
		msa = iter.execute(atoms);
		axes = iter.getSymmetryAxes();
		
		if (SymmetryTools.isRefined(msa)) {
			CoreSuperimposer imposer = new CoreSuperimposer();
			imposer.superimpose(msa);
			MultipleAlignmentScorer.calculateScores(msa);

			// Optimize the global alignment once more (as final step)
			if (this.params.getOptimization()) {
				try {
					SymmOptimizer optimizer = new SymmOptimizer(msa, axes,
							params, params.getRndSeed());
					msa = optimizer.optimize();
				} catch (RefinerFailedException e) {
					logger.info("Optimization failed:" + e.getMessage());
				}
			}
			msa.putScore("isRefined", 1.0);
		}
		return msa;
	}

	/**
	 * Analyze a single level of symmetry.
	 * 
	 * @param atoms
	 *            Atom array of the current level
	 * @return
	 * @throws StructureException
	 */
	public MultipleAlignment analyzeLevel(Atom[] atoms) throws StructureException {

		// Reset all the variables from previous calls
		reset();

		if (atoms.length < 1) {
			throw new IllegalArgumentException("Empty Atom array given.");
		}

		AFPChain selfAFP = align(atoms);

		if (refined) {
			msa = SymmetryTools.fromAFP(selfAFP, atoms);
			CoreSuperimposer imposer = new CoreSuperimposer();
			imposer.superimpose(msa);
			MultipleAlignmentScorer.calculateScores(msa);

			// STEP 5: symmetry alignment optimization
			if (this.params.getOptimization()) {
				try {
					SymmOptimizer optimizer = new SymmOptimizer(msa, axes,
							params, params.getRndSeed());
					msa = optimizer.optimize();
				} catch (RefinerFailedException e) {
					logger.debug("Optimization failed:" + e.getMessage());
				}
			}
			msa.putScore("isRefined", 1.0);
		} else {
			// Convert the optimal pairwise alignment to MSA
			MultipleAlignmentEnsemble e = new MultipleAlignmentEnsembleImpl(
					selfAFP, atoms, atoms, false);
			msa = e.getMultipleAlignment(0);
			logger.debug("Returning optimal self-alignment");
			msa.putScore("isRefined", 0.0);
		}
		return msa;
	}

	public SymmetryAxes getSymmetryAxes() {
		return axes;
	}

	/**
	 * Set the object to its construction state. This method resets all the
	 * member variables and sets the parameters to the default ones.
	 */
	private void reset() {
		refined = false;
		msa = null;
		selfAlignments = new ArrayList<AFPChain>();
		axes = new SymmetryAxes();
	}
}
