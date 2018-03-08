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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.secstruc.SecStrucElement;
import org.biojava.nbio.structure.secstruc.SecStrucTools;
import org.biojava.nbio.structure.secstruc.SecStrucType;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Iterative version of CeSymm that aims at identifying all symmetry axis of a
 * structure.
 * <p>
 * Works in the following way:
 * <ul>
 * <li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on one of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Map back all residues in a multiple alignment of the repeats.
 * </ul>
 * </li>
 *
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class CeSymmIterative {

	private static final Logger logger = LoggerFactory
			.getLogger(CeSymmIterative.class);

	private CESymmParameters params;
	private Graph<Integer, DefaultEdge> alignGraph; // cumulative
	private List<CeSymmResult> levels; // symmetry at each level

	/**
	 * For the iterative algorithm to work properly the refinement and
	 * optimization options should be turned on, because the alignment has to be
	 * consistent at every recursive step.
	 *
	 * @param param
	 *            CeSymm parameters, make sure they are cloned
	 */
	public CeSymmIterative(CESymmParameters param) {
		params = param;
		alignGraph = new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);
		levels = new ArrayList<CeSymmResult>();
	}

	/**
	 * This method uses iteratively CeSymm to calculate all symmetries in the
	 * input array of atoms and organize them in a multiple alignment of the
	 * repeats.
	 *
	 * @param atoms
	 *            atoms
	 * @return CeSymmResult
	 *
	 * @throws StructureException
	 */
	public CeSymmResult execute(Atom[] atoms) throws StructureException {

		// First iterate through all levels and then reconstruct all repeats
		iterate(atoms);
		return reconstructSymmResult(atoms);

	}

	/**
	 * This method runs iteratively the analysis of one level of symmetry with
	 * CeSymm on the input Atom array until no more symmetries exist.
	 *
	 * @param atoms
	 *            representative Atom array of the Structure
	 * @return true if any symmetry was found, false if asymmetric
	 * @throws StructureException
	 */
	private void iterate(Atom[] atoms) throws StructureException {

		logger.debug("Starting new iteration...");

		// Return if the Atom array is too short
		if ((atoms.length <= params.getWinSize()
				|| atoms.length <= params.getMinCoreLength())
				&& !levels.isEmpty()) {
			logger.debug("Aborting iteration due to insufficient Atom "
					+ "array length: %d", atoms.length);
			return;
		}

		// Return if the maximum levels of symmetry have been reached
		if (params.getSymmLevels() > 0) {
			if (levels.size() == params.getSymmLevels())
				return;
		}

		// Perform one level CeSymm alignment
		CeSymmResult result = CeSymm.analyzeLevel(atoms, params);

		if (params.getRefineMethod() == RefineMethod.NOT_REFINED
				|| !result.isSignificant()) {
			if (levels.isEmpty())
				levels.add(result);
			return;
		}

		// Generate the Atoms of one of the symmetric repeat
		Integer start = null;
		int it = 0;
		while (start == null) {
			start = result.getMultipleAlignment().getBlocks().get(0)
					.getAlignRes().get(0).get(it);
			it++;
		}
		Integer end = null;
		it = result.getMultipleAlignment().getBlocks().get(0).getAlignRes()
				.get(0).size() - 1;
		while (end == null) {
			end = result.getMultipleAlignment().getBlocks().get(0)
					.getAlignRes().get(0).get(it);
			it--;
		}
		Atom[] atomsR = Arrays.copyOfRange(atoms, start, end + 1);

		// Check the SSE requirement
		if (countHelixStrandSSE(atomsR) < params.getSSEThreshold()) {
			if (levels.isEmpty())
				levels.add(result);
			return;
		}

		// If symmetric store the residue dependencies in alignment graph
		Block b = result.getMultipleAlignment().getBlock(0);
		for (int pos = 0; pos < b.length(); pos++) {
			for (int su = 0; su < b.size() - 1; su++) {
				Integer pos1 = b.getAlignRes().get(su).get(pos);
				Integer pos2 = b.getAlignRes().get(su + 1).get(pos);
				// Add edge from lower to higher positions
				if (pos1 != null && pos2 != null) {
					alignGraph.addVertex(pos1);
					alignGraph.addVertex(pos2);
					alignGraph.addEdge(pos1, pos2);
				}
			}
		}

		// Iterate further on those Atoms (of the first repeat only)
		levels.add(result);
		iterate(atomsR);
	}

	/**
	 * After all the analysis iterations have finished, the final Result object
	 * is reconstructed using the cumulative alignment graph.
	 *
	 * @param atoms
	 *            the original structure atoms
	 * @return CeSymmResult reconstructed symmetry result
	 * @throws StructureException
	 */
	private CeSymmResult reconstructSymmResult(Atom[] atoms)
			throws StructureException {
		
		// If one level, nothing to build or calculate
		if (levels.size() == 1)
			return levels.get(0);
		
		CeSymmResult result = new CeSymmResult();
		result.setSelfAlignment(levels.get(0).getSelfAlignment());
		result.setStructureId(levels.get(0).getStructureId());
		result.setAtoms(levels.get(0).getAtoms());
		result.setParams(levels.get(0).getParams());

		// Initialize a new multiple alignment
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.getEnsemble().setStructureIdentifiers(
				new ArrayList<StructureIdentifier>());
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		// Calculate the connected groups of the alignment graph
		ConnectivityInspector<Integer, DefaultEdge> inspector = new ConnectivityInspector<Integer, DefaultEdge>(
				alignGraph);
		List<Set<Integer>> comps = inspector.connectedSets();
		List<ResidueGroup> groups = new ArrayList<ResidueGroup>(comps.size());
		for (Set<Integer> comp : comps)
			groups.add(new ResidueGroup(comp));

		// Calculate the total number of repeats
		int order = 1;
		for (CeSymmResult sr : levels)
			order *= sr.getMultipleAlignment().size();
		for (int su = 0; su < order; su++)
			b.getAlignRes().add(new ArrayList<Integer>());

		// Construct the resulting MultipleAlignment from ResidueGroups
		for (ResidueGroup group : groups) {
			if (group.order() != order)
				continue;
			group.combineWith(b.getAlignRes());
		}

		// The reconstruction failed, so the top level is returned
		if (b.length() == 0)
			return levels.get(0);

		for (int su = 0; su < order; su++) {
			Collections.sort(b.getAlignRes().get(su));
			msa.getEnsemble().getAtomArrays().add(atoms);
			msa.getEnsemble().getStructureIdentifiers()
					.add(result.getStructureId());
		}

		result.setMultipleAlignment(msa);
		result.setRefined(true);
		result.setNumRepeats(order);

		SymmetryAxes axes = recoverAxes(result);
		result.setAxes(axes);

		// Set the transformations and scores of the final alignment
		SymmetryTools
				.updateSymmetryTransformation(result.getAxes(), msa);
		double tmScore = MultipleAlignmentScorer.getAvgTMScore(msa)
				* msa.size();
		double rmsd = MultipleAlignmentScorer.getRMSD(msa);
		msa.putScore(MultipleAlignmentScorer.AVGTM_SCORE, tmScore);
		msa.putScore(MultipleAlignmentScorer.RMSD, rmsd);

		return result;
	}

	/**
	 * The symmetry axes of each level are recovered after the symmetry analysis
	 * iterations have finished, using the stored MultipleAlignment at each
	 * symmetry level.
	 * @return SymmetryAxes
	 */
	private SymmetryAxes recoverAxes(CeSymmResult result) {

		SymmetryAxes axes = new SymmetryAxes();

		for (int m = 0; m < levels.size(); m++) {

			MultipleAlignment align = levels.get(m).getMultipleAlignment();
			Matrix4d axis = align.getBlockSet(0).getTransformations().get(1);
			SymmetryType type = levels.get(m).getAxes().getElementaryAxis(0).getSymmType();
			int order = align.size();

			axes.addAxis(axis, order, type);
		}
		return axes;
	}

	/**
	 * Calculate the number of helix and strand SSE of a repeat.
	 *
	 * @param atoms
	 *            Atom array of the repeat found
	 * @return int number of helix or strand SSE
	 */
	private static int countHelixStrandSSE(Atom[] atoms) {

		List<SecStrucElement> sses = SecStrucTools
				.getSecStrucElements(SymmetryTools.getGroups(atoms));
		int count = 0;

		//keep track of different helix types
		boolean helix = false;
		int hEnd = 0;

		for (SecStrucElement sse : sses) {
			SecStrucType t = sse.getType();
			if (t.isBetaStrand()) {
				helix = false;
				count++;
			} else if (t.isHelixType()){
				if (helix){
					// If this helix is contiguous to the previous
					if (sse.getRange().getStart().getSeqNum() + 1 == hEnd)
						hEnd = sse.getRange().getEnd().getSeqNum();
					else
						count++;
				} else
					count++;
			} else
				helix = false;
		}
		return count;
	}

}
