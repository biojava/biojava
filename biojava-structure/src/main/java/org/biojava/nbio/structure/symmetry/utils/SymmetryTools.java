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
package org.biojava.nbio.structure.symmetry.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CECalculator;
import org.biojava.nbio.structure.align.helper.AlignTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsemble;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.CoreSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.align.multiple.util.MultipleSuperimposer;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.Subunits;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Utility methods for the internal symmetry identification and manipulation.
 * <p>
 * Methods include: blank out regions of DP Matrix, build symmetry graphs, get
 * rotation symmetry angles, split repeats in quaternary structure chains,
 * convert between symmetry formats (full, repeats, rotations), determine if two
 * symmetry axes are equivalent, get groups from representative Atoms.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 *
 */
public class SymmetryTools {

	private static final Logger logger = LoggerFactory
			.getLogger(SymmetryTools.class);

	/** Prevent instantiation. */
	private SymmetryTools() {
	}

	/**
	 * Returns the "reset value" for graying out the main diagonal. If we're
	 * blanking out the main diagonal, this value is always Integer.MIN_VALUE.
	 * <p>
	 * This is possible if {@code gradientPolyCoeff = Integer.MIN_VALUE} and
	 * {@code gradientExpCoeff = 0}.
	 *
	 * @param unpenalizedScore
	 * @param nResFromMainDiag
	 * @param gradientPolyCoeff
	 * @param gradientExpCoeff
	 * @return
	 */
	private static double getResetVal(double unpenalizedScore,
			double nResFromMainDiag, double[] gradientPolyCoeff,
			double gradientExpCoeff) {

		if (unpenalizedScore == Double.NaN)
			return 0; // what else?

		// We can actually return a positive value if this is high enough
		double updateVal = unpenalizedScore;
		updateVal -= gradientExpCoeff * Math.pow(Math.E, -nResFromMainDiag);
		for (int p = 0; p < gradientPolyCoeff.length; p++) {
			updateVal -= gradientPolyCoeff[gradientPolyCoeff.length - 1 - p]
					* Math.pow(nResFromMainDiag, -p);
		}
		return updateVal;
	}

	/**
	 * Grays out the main diagonal of a duplicated distance matrix.
	 *
	 * @param ca2
	 * @param rows
	 *            Number of rows
	 * @param cols
	 *            Number of original columns
	 * @param calculator
	 *            Used to get the matrix if origM is null
	 * @param origM
	 *            starting matrix. If null, uses
	 *            {@link CECalculator#getMatMatrix()}
	 * @param blankWindowSize
	 *            Width of section to gray out
	 * @param gradientPolyCoeff
	 * @param gradientExpCoeff
	 * @return
	 */
	public static Matrix grayOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize,
			double[] gradientPolyCoeff, double gradientExpCoeff) {

		if (origM == null) {
			origM = new Matrix(calculator.getMatMatrix());
		}

		// symmetry hack, disable main diagonal

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				int diff = Math.abs(i - j);

				double resetVal = getResetVal(origM.get(i, j), diff,
						gradientPolyCoeff, gradientExpCoeff);

				if (diff < blankWindowSize) {
					origM.set(i, j, origM.get(i, j) + resetVal);

				}
				int diff2 = Math.abs(i - (j - ca2.length / 2)); // other side

				double resetVal2 = getResetVal(origM.get(i, j), diff2,
						gradientPolyCoeff, gradientExpCoeff);

				if (diff2 < blankWindowSize) {
					origM.set(i, j, origM.get(i, j) + resetVal2);

				}
			}
		}
		return origM;
	}

	public static Matrix grayOutPreviousAlignment(AFPChain afpChain,
			Atom[] ca2, int rows, int cols, CECalculator calculator,
			Matrix max, int blankWindowSize, double[] gradientPolyCoeff,
			double gradientExpCoeff) {

		max = grayOutCEOrig(ca2, rows, cols, calculator, max, blankWindowSize,
				gradientPolyCoeff, gradientExpCoeff);

		double[][] dist1 = calculator.getDist1();
		double[][] dist2 = calculator.getDist2();

		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated
		int breakPoint = ca2.length / 2;
		for (int bk = 0; bk < blockNum; bk++) {

			for (int i = 0; i < optLen[bk]; i++) {
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];

				int dist = blankWindowSize / 2;
				int start1 = Math.max(pos1 - dist, 0);
				int start2 = Math.max(pos2 - dist, 0);
				int end1 = Math.min(pos1 + dist, rows - 1);
				int end2 = Math.min(pos2 + dist, cols - 1);

				for (int i1 = start1; i1 < end1; i1++) {

					// blank diagonal of dist1
					for (int k = 0; k < blankWindowSize / 2; k++) {
						if (i1 - k >= 0) {
							double resetVal = getResetVal(
									max.get(i1 - k, i1 - k), 0,
									gradientPolyCoeff, gradientExpCoeff);
							dist1[i1 - k][i1 - k] = resetVal;
						} else if (i1 + k < rows) {
							double resetVal = getResetVal(
									max.get(i1 + k, i1 + k), 0,
									gradientPolyCoeff, gradientExpCoeff);
							dist1[i1 + k][i1 + k] = resetVal;
						}

					}

					for (int j2 = start2; j2 < end2; j2++) {
						double resetVal = getResetVal(max.get(i1, j2),
								Math.abs(i1 - j2), gradientPolyCoeff,
								gradientExpCoeff);
						max.set(i1, j2, resetVal);
						if (j2 < breakPoint) {
							double resetVal2 = getResetVal(
									max.get(i1, j2 + breakPoint),
									Math.abs(i1 - (j2 + breakPoint)),
									gradientPolyCoeff, gradientExpCoeff);
							max.set(i1, j2 + breakPoint, resetVal2);
						} else {
							double resetVal2 = getResetVal(
									max.get(i1, j2 - breakPoint),
									Math.abs(i1 - (j2 - breakPoint)),
									gradientPolyCoeff, gradientExpCoeff);
							max.set(i1, j2 - breakPoint, resetVal2);
						}
						for (int k = 0; k < blankWindowSize / 2; k++) {
							if (j2 - k >= 0) {
								if (j2 - k < breakPoint) {
									double resetVal2 = getResetVal(
											max.get(j2 - k, j2 - k), 0,
											gradientPolyCoeff, gradientExpCoeff);
									dist2[j2 - k][j2 - k] = resetVal2;
								} else {
									double resetVal2 = getResetVal(max.get(j2
											- k - breakPoint, j2 - k), 0,
											gradientPolyCoeff, gradientExpCoeff);
									dist2[j2 - k - breakPoint][j2 - k
											- breakPoint] = resetVal2;
								}
							} else if (j2 + k < cols) {
								if (j2 + k < breakPoint) {
									double resetVal2 = getResetVal(
											max.get(j2 + k, j2 + k), 0,
											gradientPolyCoeff, gradientExpCoeff);
									dist2[j2 + k][j2 + k] = resetVal2;
								} else {
									double resetVal2 = getResetVal(max.get(j2
											+ k - breakPoint, j2 + k), 0,
											gradientPolyCoeff, gradientExpCoeff);
									dist2[j2 + k - breakPoint][j2 + k
											- breakPoint] = resetVal2;
								}
							}
						}
					}
				}

			}
		}
		calculator.setDist1(dist1);
		calculator.setDist2(dist2);
		return max;

	}

	public Matrix getDkMatrix(Atom[] ca1, Atom[] ca2, int fragmentLength,
			double[] dist1, double[] dist2, int rows, int cols) {

		Matrix diffDistMax = Matrix.identity(ca1.length, ca2.length);

		for (int i = 0; i < rows; i++) {
			double score1 = 0;
			for (int x = 0; x < fragmentLength; x++) {
				score1 += dist1[i + x];
			}
			for (int j = 0; j < cols; j++) {
				double score2 = 0;
				for (int y = 0; y < fragmentLength; y++) {
					score2 += dist2[j + y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar,
				// i.e. the difference is close to 0
				diffDistMax.set(i, j, Math.abs(score1 - score2));
			}
		}

		// symmetry hack, disable main diagonal

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				int diff = Math.abs(i - j);

				if (diff < 15) {
					diffDistMax.set(i, j, 99);
				}
				int diff2 = Math.abs(i - (j - ca2.length / 2));
				if (diff2 < 15) {
					diffDistMax.set(i, j, 99);
				}
			}
		}
		return diffDistMax;

	}

	public static Matrix blankOutPreviousAlignment(AFPChain afpChain,
			Atom[] ca2, int rows, int cols, CECalculator calculator,
			Matrix max, int blankWindowSize) {
		return grayOutPreviousAlignment(afpChain, ca2, rows, cols, calculator,
				max, blankWindowSize, new double[] { Integer.MIN_VALUE }, 0.0);
	}

	public static Matrix blankOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM, int blankWindowSize) {
		return grayOutCEOrig(ca2, rows, cols, calculator, origM,
				blankWindowSize, new double[] { Integer.MIN_VALUE }, 0.0);
	}

	public static Matrix getDkMatrix(Atom[] ca1, Atom[] ca2, int k,
			int fragmentLength) {

		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);
		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length - fragmentLength - k + 1;
		int cols = ca2.length - fragmentLength - k + 1;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		Matrix m2 = new Matrix(rows, cols);

		for (int i = 0; i < rows; i++) {
			double score1 = 0;
			for (int x = 0; x < fragmentLength; x++) {
				score1 += dist1[i + x];
			}
			for (int j = 0; j < cols; j++) {
				double score2 = 0;
				for (int y = 0; y < fragmentLength; y++) {
					score2 += dist2[j + y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar,
				// i.e. the difference is close to 0
				m2.set(i, j, Math.abs(score1 - score2));
			}
		}
		return m2;
	}

	public static boolean[][] blankOutBreakFlag(AFPChain afpChain, Atom[] ca2,
			int rows, int cols, CECalculator calculator, boolean[][] breakFlag,
			int blankWindowSize) {

		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();

		// ca2 is circularly permutated at this point.
		int breakPoint = ca2.length / 2;

		for (int bk = 0; bk < blockNum; bk++) {

			// Matrix m= afpChain.getBlockRotationMatrix()[bk];
			// Atom shift = afpChain.getBlockShiftVector()[bk];
			for (int i = 0; i < optLen[bk]; i++) {
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];
				// blank out area around these positions...

				int dist = blankWindowSize;
				int start1 = Math.max(pos1 - dist, 0);
				int start2 = Math.max(pos2 - dist, 0);
				int end1 = Math.min(pos1 + dist, rows - 1);
				int end2 = Math.min(pos2 + dist, cols - 1);

				for (int i1 = start1; i1 < end1; i1++) {

					for (int j2 = start2; j2 < end2; j2++) {
						// System.out.println(i1 + " " + j2 + " (***)");
						breakFlag[i1][j2] = true;
						if (j2 < breakPoint) {
							breakFlag[i1][j2 + breakPoint] = true;
						}
					}
				}

			}
		}

		return breakFlag;
	}

	/**
	 * Returns the <em>magnitude</em> of the angle between the first and second
	 * blocks of {@code afpChain}, measured in degrees. This is always a
	 * positive value (unsigned).
	 *
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	public static double getAngle(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		Matrix rotation = afpChain.getBlockRotationMatrix()[0];
		return Math.acos(rotation.trace() - 1) * 180 / Math.PI;
	}

	/**
	 * Converts a set of AFP alignments into a Graph of aligned residues, where
	 * each vertex is a residue and each edge means the connection between the
	 * two residues in one of the alignments.
	 *
	 * @param afps
	 *            List of AFPChains
	 * @param atoms
	 *            Atom array of the symmetric structure
	 * @param undirected
	 *            if true, the graph is undirected
	 *
	 * @return adjacency List of aligned residues
	 */
	public static List<List<Integer>> buildSymmetryGraph(List<AFPChain> afps,
			Atom[] atoms, boolean undirected) {

		List<List<Integer>> graph = new ArrayList<List<Integer>>();

		for (int n = 0; n < atoms.length; n++) {
			graph.add(new ArrayList<Integer>());
		}

		for (int k = 0; k < afps.size(); k++) {
			for (int i = 0; i < afps.get(k).getOptAln().length; i++) {
				for (int j = 0; j < afps.get(k).getOptAln()[i][0].length; j++) {
					Integer res1 = afps.get(k).getOptAln()[i][0][j];
					Integer res2 = afps.get(k).getOptAln()[i][1][j];
					graph.get(res1).add(res2);
					if (undirected)
						graph.get(res2).add(res1);
				}
			}
		}
		return graph;
	}

	/**
	 * Converts a self alignment into a directed jGraphT of aligned residues,
	 * where each vertex is a residue and each edge means the equivalence
	 * between the two residues in the self-alignment.
	 *
	 * @param selfAlignment
	 *            AFPChain
	 *
	 * @return alignment Graph
	 */
	public static UndirectedGraph<Integer, DefaultEdge> buildSymmetryGraph(
			AFPChain selfAlignment) {

		UndirectedGraph<Integer, DefaultEdge> graph = new SimpleGraph<Integer, DefaultEdge>(
				DefaultEdge.class);

		for (int i = 0; i < selfAlignment.getOptAln().length; i++) {
			for (int j = 0; j < selfAlignment.getOptAln()[i][0].length; j++) {
				Integer res1 = selfAlignment.getOptAln()[i][0][j];
				Integer res2 = selfAlignment.getOptAln()[i][1][j];
				graph.addVertex(res1);
				graph.addVertex(res2);
				graph.addEdge(res1, res2);
			}
		}
		return graph;
	}

	/**
	 * Method that converts the symmetric units of a structure into different
	 * structures, so that they can be individually visualized.
	 *
	 * @param symmetry
	 *            CeSymmResult
	 * @throws StructureException
	 * @result List of structures, by repeat index sequentially
	 * 
	 */
	public static List<Structure> divideStructure(CeSymmResult symmetry)
			throws StructureException {

		if (!symmetry.isRefined())
			throw new IllegalArgumentException("The symmetry result "
					+ "is not refined, repeats cannot be defined");

		int order = symmetry.getMultipleAlignment().size();
		Atom[] atoms = StructureTools.cloneAtomArray(symmetry.getAtoms());
		List<StructureIdentifier> repeatsId = symmetry.getRepeatsID();
		List<Structure> repeats = new ArrayList<Structure>(order);

		// Create new structure containing the repeat atoms
		for (int i = 0; i < order; i++) {

			Structure s = new StructureImpl();
			s.setStructureIdentifier(repeatsId.get(i));

			Block align = symmetry.getMultipleAlignment().getBlock(0);

			// Get the start and end of the repeat
			int res1 = align.getStartResidue(i);
			int res2 = align.getFinalResidue(i);

			Atom[] repeat = Arrays.copyOfRange(atoms, res1, res2 + 1);

			Chain newCh = new ChainImpl();
			newCh.setChainID(repeat[0].getGroup().getChainId());

			for (int k = 0; k < repeat.length; k++) {
				Group g = (Group) repeat[k].getGroup().clone();
				newCh.addGroup(g);
			}
			s.addChain(newCh);
			repeats.add(s);
		}
		return repeats;
	}

	/**
	 * Method that converts a repeats symmetric alignment into an alignment of
	 * whole structures.
	 * <p>
	 * Example: if the structure has repeats A,B and C, the original alignment
	 * is A-B-C, and the returned alignment is ABC-BCA-CAB.
	 *
	 * @param symm
	 *            CeSymmResult
	 * @return MultipleAlignment of the full structure superpositions
	 */
	public static MultipleAlignment toFullAlignment(CeSymmResult symm) {

		if (!symm.isRefined())
			throw new IllegalArgumentException("The symmetry result "
					+ "is not refined, repeats cannot be defined");

		MultipleAlignment full = symm.getMultipleAlignment().clone();

		for (int str = 1; str < full.size(); str++) {
			// Create a new Block with swapped AlignRes (move first to last)
			Block b = full.getBlock(full.getBlocks().size() - 1).clone();
			b.getAlignRes().add(b.getAlignRes().get(0));
			b.getAlignRes().remove(0);
			full.getBlockSet(0).getBlocks().add(b);
		}
		return full;
	}

	/**
	 * Method that converts a symmetry alignment into an alignment of the
	 * repeats only, as new independent structures.
	 * <p>
	 * This method changes the structure identifiers, the Atom arrays and
	 * re-scles the aligned residues in the Blocks corresponding to those
	 * changes.
	 * <p>
	 * Application: display superimposed repeats in Jmol.
	 *
	 * @param result
	 *            CeSymmResult of symmetry
	 * @return MultipleAlignment of the repeats
	 * @throws StructureException
	 */
	public static MultipleAlignment toRepeatsAlignment(CeSymmResult result)
			throws StructureException {

		if (!result.isRefined())
			throw new IllegalArgumentException("The symmetry result "
					+ "is not refined, repeats cannot be defined");

		MultipleAlignment msa = result.getMultipleAlignment();
		MultipleAlignmentEnsemble newEnsemble = msa.getEnsemble().clone();

		List<Structure> repSt = SymmetryTools.divideStructure(result);

		MultipleAlignment repeats = newEnsemble.getMultipleAlignment(0);
		Block block = repeats.getBlock(0);
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();

		for (Structure s : repSt)
			atomArrays.add(StructureTools.getRepresentativeAtomArray(s));

		newEnsemble.setAtomArrays(atomArrays);

		for (int su = 0; su < block.size(); su++) {
			Integer start = block.getStartResidue(su);

			// Normalize aligned residues
			for (int res = 0; res < block.length(); res++) {
				Integer residue = block.getAlignRes().get(su).get(res);
				if (residue != null)
					residue -= start;
				block.getAlignRes().get(su).set(res, residue);
			}
		}

		return repeats;
	}

	/**
	 * Converts a refined symmetry AFPChain alignment into the standard
	 * representation of symmetry in a MultipleAlignment, that contains the
	 * entire Atom array of the strcuture and the symmetric repeats are orgaized
	 * in different rows in a single Block.
	 *
	 * @param symm
	 *            AFPChain created with a symmetry algorithm and refined
	 * @param atoms
	 *            Atom array of the entire structure
	 * @return MultipleAlignment format of the symmetry
	 * @throws StructureException
	 */
	public static MultipleAlignment fromAFP(AFPChain symm, Atom[] atoms)
			throws StructureException {

		if (!symm.getAlgorithmName().contains("symm")) {
			throw new IllegalArgumentException(
					"The input alignment is not a symmetry alignment.");
		}

		MultipleAlignmentEnsemble e = new MultipleAlignmentEnsembleImpl(symm,
				atoms, atoms, false);
		e.setAtomArrays(new ArrayList<Atom[]>());
		StructureIdentifier name = null;
		if (e.getStructureIdentifiers() != null) {
			if (!e.getStructureIdentifiers().isEmpty())
				name = e.getStructureIdentifiers().get(0);
		} else
			name = atoms[0].getGroup().getChain().getStructure()
					.getStructureIdentifier();

		e.setStructureIdentifiers(new ArrayList<StructureIdentifier>());

		MultipleAlignment result = new MultipleAlignmentImpl();
		BlockSet bs = new BlockSetImpl(result);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		int order = symm.getBlockNum();
		for (int su = 0; su < order; su++) {
			List<Integer> residues = e.getMultipleAlignment(0).getBlock(su)
					.getAlignRes().get(0);
			b.getAlignRes().add(residues);
			e.getStructureIdentifiers().add(name);
			e.getAtomArrays().add(atoms);
		}
		e.getMultipleAlignments().set(0, result);
		result.setEnsemble(e);

		CoreSuperimposer imposer = new CoreSuperimposer();
		imposer.superimpose(result);
		updateSymmetryScores(result);

		return result;
	}

	/**
	 * Determines if two symmetry axis are equivalent inside the error
	 * threshold. It only takes into account the direction of the vector where
	 * the rotation is made: the angle and translation are not taken into
	 * account.
	 *
	 * @param axis1
	 * @param axis2
	 * @param epsilon
	 *            error allowed in the axis comparison
	 * @return true if equivalent, false otherwise
	 */
	@Deprecated
	public static boolean equivalentAxes(Matrix4d axis1, Matrix4d axis2,
			double epsilon) {

		AxisAngle4d rot1 = new AxisAngle4d();
		rot1.set(axis1);
		AxisAngle4d rot2 = new AxisAngle4d();
		rot2.set(axis2);

		// rot1.epsilonEquals(rot2, error); //that also compares angle
		// L-infinite distance without comparing the angle (epsilonEquals)
		List<Double> sameDir = new ArrayList<Double>();
		sameDir.add(Math.abs(rot1.x - rot2.x));
		sameDir.add(Math.abs(rot1.y - rot2.y));
		sameDir.add(Math.abs(rot1.z - rot2.z));

		List<Double> otherDir = new ArrayList<Double>();
		otherDir.add(Math.abs(rot1.x + rot2.x));
		otherDir.add(Math.abs(rot1.y + rot2.y));
		otherDir.add(Math.abs(rot1.z + rot2.z));

		Double error = Math.min(Collections.max(sameDir),
				Collections.max(otherDir));

		return error < epsilon;
	}

	/**
	 * Given a symmetry result, it calculates the overall global symmetry,
	 * factoring out the alignment and detection steps of
	 * {@link QuatSymmetryDetector} algorithm.
	 *
	 * @param result
	 *            symmetry result
	 * @return global symmetry results
	 * @throws StructureException
	 */
	public static QuatSymmetryResults getQuaternarySymmetry(CeSymmResult result)
			throws StructureException {

		// Obtain the clusters of aligned Atoms and repeat variables
		List<Atom[]> alignedCA = result.getMultipleAlignment().getAtomArrays();
		MultipleAlignment msa = result.getMultipleAlignment();
		List<Integer> corePos = MultipleAlignmentTools.getCorePositions(result
				.getMultipleAlignment().getBlock(0));

		List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
		List<Integer> folds = new ArrayList<Integer>();
		List<Boolean> pseudo = new ArrayList<Boolean>();
		List<String> chainIds = new ArrayList<String>();
		List<Integer> models = new ArrayList<Integer>();
		List<Double> seqIDmin = new ArrayList<Double>();
		List<Double> seqIDmax = new ArrayList<Double>();
		List<Integer> clusterIDs = new ArrayList<Integer>();
		int fold = 1;

		for (int str = 0; str < alignedCA.size(); str++) {
			Atom[] array = alignedCA.get(str);
			List<Point3d> points = new ArrayList<Point3d>();
			List<Integer> alignedRes = msa.getBlock(0).getAlignRes().get(str);
			for (int pos = 0; pos < alignedRes.size(); pos++) {
				Integer residue = alignedRes.get(pos);
				if (residue == null)
					continue;
				else if (!corePos.contains(pos))
					continue;
				Atom a = array[residue];
				points.add(new Point3d(a.getCoords()));
			}
			caCoords.add(points.toArray(new Point3d[points.size()]));
			if (alignedCA.size() % fold == 0) {
				folds.add(fold); // the folds are the common denominators
			}
			fold++;
			pseudo.add(false);
			chainIds.add(alignedCA.get(str)[0].getGroup().getChainId());
			models.add(0);
			seqIDmax.add(1.0);
			seqIDmin.add(1.0);
			clusterIDs.add(0);
		}

		// Create directly the repeats, because we know the aligned CA
		Subunits globalSubunits = new Subunits(caCoords, clusterIDs, pseudo,
				seqIDmin, seqIDmax, folds, chainIds, models);

		// Quaternary Symmetry Detection
		QuatSymmetryParameters param = new QuatSymmetryParameters();

		QuatSymmetryResults gSymmetry = QuatSymmetryDetector.calcQuatSymmetry(
				globalSubunits, param);

		return gSymmetry;
	}

	/**
	 * Returns true a symmetry multiple alignment has been refined, false
	 * otherwise.
	 * <p>
	 * For a refined alignment only one Block with no repeated residues is
	 * necessary. Sufficient condition is not tested (only known from the
	 * algorithm or CeSymmResult).
	 *
	 * @param symm
	 *            the symmetry alignment
	 * @return true if the alignment is refined
	 */
	@Deprecated
	public static boolean isRefined(MultipleAlignment symm) {

		if (symm.getBlocks().size() > 1) {
			return false;
		} else if (symm.size() < 2)
			return false;
		else {
			List<Integer> alreadySeen = new ArrayList<Integer>();
			List<List<Integer>> align = symm.getBlock(0).getAlignRes();
			for (int str = 0; str < symm.size(); str++) {
				for (int res = 0; res < align.get(str).size(); res++) {
					Integer residue = align.get(str).get(res);
					if (residue == null)
						continue;
					if (alreadySeen.contains(residue)) {
						return false;
					} else {
						alreadySeen.add(residue);
					}
				}
			} // end of all repeats
			return true;
		}
	}

	/**
	 * Returns true if the symmetry alignment is significant, false otherwise.
	 * <p>
	 * For a symmetry alignment to be significant, the alignment has to be
	 * refined and the TM-score has to be higher than the threshold.
	 * <p>
	 * It is recommended to use the {@link CeSymmResult#isSignificant()} method
	 * instead.
	 *
	 * @param msa
	 * @param symmetryThreshold
	 * @return
	 * @throws StructureException
	 */
	@Deprecated
	public static boolean isSignificant(MultipleAlignment msa,
			double symmetryThreshold) throws StructureException {

		// Order/refinement check
		if (!SymmetryTools.isRefined(msa))
			return false;

		// TM-score cutoff
		double tm = 0.0;
		if (msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE) == null)
			tm = MultipleAlignmentScorer.getAvgTMScore(msa);
		else
			tm = msa.getScore(MultipleAlignmentScorer.AVGTM_SCORE);

		if (tm < symmetryThreshold)
			return false;

		return true;
	}

	/**
	 * Returns the List of Groups of the corresponding representative Atom
	 * array. The representative Atom array needs to fulfill: no two Atoms are
	 * from the same Group and Groups are sequential (connected in the original
	 * Structure), except if they are from different Chains.
	 *
	 * @param rAtoms
	 *            array of representative Atoms (CA, P, etc).
	 * @return List of Groups
	 */
	public static List<Group> getGroups(Atom[] rAtoms) {

		List<Group> groups = new ArrayList<Group>(rAtoms.length);

		for (Atom a : rAtoms) {
			Group g = a.getGroup();
			if (g != null)
				groups.add(g);
			else
				logger.info("Group not found for representative Atom {}", a);
		}
		return groups;
	}

	/**
	 * Calculates the set of symmetry operation Matrices (transformations) of
	 * the new alignment, based on the symmetry relations in the SymmetryAxes
	 * object. It sets the transformations to the input MultipleAlignment and
	 * SymmetryAxes objects. If the SymmetryAxes object is null, the
	 * superposition of the repeats is done without symmetry constraints.
	 * <p>
	 * This method also sets the scores (RMSD and TM-score) after the new
	 * superposition has been updated.
	 *
	 * @param axes
	 *            SymmetryAxes object. It will be modified.
	 * @param msa
	 *            MultipleAlignment. It will be modified.
	 */
	public static void updateSymmetryTransformation(SymmetryAxes axes,
			MultipleAlignment msa) throws StructureException {

		List<List<Integer>> block = msa.getBlocks().get(0).getAlignRes();
		int length = block.get(0).size();
		Atom[] atoms = StructureTools.cloneAtomArray(msa.getEnsemble()
				.getAtomArrays().get(0));

		if (axes != null) {
			for (int level = 0; level < axes.getNumLevels(); level++) {

				// Calculate the aligned atom arrays to superimpose
				List<Atom> list1 = new ArrayList<Atom>();
				List<Atom> list2 = new ArrayList<Atom>();

				for (int firstRepeat : axes.getFirstRepeats(level)) {

					Matrix4d transform = axes.getRepeatTransform(firstRepeat);

					List<List<Integer>> relation = axes.getRepeatRelation(
							level, firstRepeat);

					for (int index = 0; index < relation.get(0).size(); index++) {
						int p1 = relation.get(0).get(index);
						int p2 = relation.get(1).get(index);

						for (int k = 0; k < length; k++) {
							Integer pos1 = block.get(p1).get(k);
							Integer pos2 = block.get(p2).get(k);
							if (pos1 != null && pos2 != null) {
								Calc.transform(atoms[pos1], transform);
								Calc.transform(atoms[pos2], transform);
								list1.add(atoms[pos1]);
								list2.add(atoms[pos2]);
							}
						}
					}
				}

				Atom[] arr1 = list1.toArray(new Atom[list1.size()]);
				Atom[] arr2 = list2.toArray(new Atom[list2.size()]);

				// Calculate the new transformation information
				if (arr1.length > 0 && arr2.length > 0) {
					SVDSuperimposer svd = new SVDSuperimposer(arr1, arr2);
					Matrix4d axis = svd.getTransformation();
					axes.updateAxis(level, axis);
				}

				// Get the transformations from the SymmetryAxes
				List<Matrix4d> transformations = new ArrayList<Matrix4d>();
				for (int su = 0; su < msa.size(); su++) {
					transformations.add(axes.getRepeatTransform(su));
				}
				msa.getBlockSet(0).setTransformations(transformations);
			}
		} else {
			MultipleSuperimposer imposer = new CoreSuperimposer();
			imposer.superimpose(msa);
		}
		updateSymmetryScores(msa);
	}

	/**
	 * Update the scores (TM-score and RMSD) of a symmetry multiple alignment.
	 * This method does not redo the superposition of the alignment.
	 *
	 * @param symm
	 *            Symmetry Multiple Alignment of Repeats
	 * @throws StructureException
	 */
	public static void updateSymmetryScores(MultipleAlignment symm)
			throws StructureException {

		// Multiply by the order of symmetry to normalize score
		double tmScore = MultipleAlignmentScorer.getAvgTMScore(symm)
				* symm.size();
		double rmsd = MultipleAlignmentScorer.getRMSD(symm);

		symm.putScore(MultipleAlignmentScorer.AVGTM_SCORE, tmScore);
		symm.putScore(MultipleAlignmentScorer.RMSD, rmsd);
	}

	/**
	 * Returns the representative Atom Array of the first model, if the
	 * structure is NMR, or the Array for each model, if it is a biological
	 * assembly with multiple models.
	 * 
	 * @param structure
	 * @return representative Atom[]
	 */
	public static Atom[] getRepresentativeAtoms(Structure structure) {

		if (structure.isNmr())
			return StructureTools.getRepresentativeAtomArray(structure);

		else {

			// Get Atoms of all models
			List<Atom> atomList = new ArrayList<Atom>();
			for (int m = 0; m < structure.nrModels(); m++) {
				for (Chain c : structure.getModel(m))
					atomList.addAll(Arrays.asList(StructureTools
							.getRepresentativeAtomArray(c)));
			}
			return atomList.toArray(new Atom[0]);
		}

	}

}
