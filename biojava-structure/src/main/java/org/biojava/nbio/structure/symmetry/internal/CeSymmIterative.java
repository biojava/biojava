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
import org.biojava.nbio.structure.PassthroughIdentifier;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.CoreSuperimposer;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.secstruc.SecStrucElement;
import org.biojava.nbio.structure.secstruc.SecStrucTools;
import org.biojava.nbio.structure.secstruc.SecStrucType;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.jgrapht.DirectedGraph;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.experimental.dag.DirectedAcyclicGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Iterative version of CeSymm that aims at identifying all symmetry axis
 * (internal or quaternary) of a particular structure.
 * <p>
 * Works in the following way:
 * <ul>
 * <li>Run CeSymm on the original structure.
 * <li>Calculate the symmetric unit boundaries.
 * <li>Run CeSymm on one of the symmetric units to find further symmetries.
 * <li>Repeat the last two steps until no more significant results are found.
 * <li>Map back all residues in a multiple alignment of the subunits.
 * <li>Run a final optimization of all symmetric units correctly superimposed.
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
	private Atom[] allAtoms;
	private UndirectedGraph<Integer, DefaultEdge> alignGraph; // alignment graph
	private List<MultipleAlignment> levels; // msa at each level
	private CeSymmResult result;

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
		levels = new ArrayList<MultipleAlignment>();
	}

	/**
	 * This method uses iteratively CeSymm to calculate all symmetries in the
	 * input array of atoms and organize them in a multiple alignment of the
	 * subunits.
	 * 
	 * @param atoms
	 *            atoms
	 * @return CeSymmResult
	 * 
	 * @throws StructureException
	 */
	public CeSymmResult execute(Atom[] atoms) throws StructureException {

		allAtoms = atoms;

		// True if symmetry found
		boolean symm = iterate(atoms);

		if (symm) {
			buildAlignment();
			recoverAxes();
		}

		return result;
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
	private boolean iterate(Atom[] atoms) throws StructureException {

		logger.debug("Starting new iteration...");

		// Return if the Atom array is too short
		if (atoms.length <= params.getWinSize()
				|| atoms.length <= params.getMinCoreLength()) {
			logger.debug("Aborting iteration due to insufficient Atom "
					+ "array length: %d", atoms.length);
			return !levels.isEmpty();
		}

		// Return if the maximum levels of symmetry have been reached
		if (params.getSymmLevels() > 0) {
			if (levels.size() == params.getSymmLevels())
				return true;
		}

		// Perform one level CeSymm alignment
		CeSymmResult r = CeSymm.analyzeLevel(atoms, params);
		if (result == null)
			result = r;

		if (params.getRefineMethod() == RefineMethod.NOT_REFINED)
			return false;
		else if (!r.isSignificant())
			return !levels.isEmpty();

		// Generate the Atoms of one of the symmetric subunit
		Integer start = null;
		int it = 0;
		while (start == null) {
			start = r.getMultipleAlignment().getBlocks().get(0)
					.getAlignRes().get(0).get(it);
			it++;
		}
		Integer end = null;
		it = r.getMultipleAlignment().getBlocks().get(0).getAlignRes()
				.get(0).size() - 1;
		while (end == null) {
			end = r.getMultipleAlignment().getBlocks().get(0)
					.getAlignRes().get(0).get(it);
			it--;
		}
		Atom[] atomsR = Arrays.copyOfRange(atoms, start, end + 1);

		// Check the SSE requirement
		if (countHelixStrandSSE(atomsR) < params.getSSEThreshold())
			return !levels.isEmpty();

		// If symmetric store the residue dependencies in alignment graph
		Block b = r.getMultipleAlignment().getBlock(0);
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

		// Iterate further on those Atoms (of the first subunit only)
		levels.add(r.getMultipleAlignment());
		return iterate(atomsR);
	}

	/**
	 * After all the analysis iteratives have finished, the final
	 * MultipleAlignment object is constructed using the alignment graph.
	 * 
	 * @throws StructureException
	 */
	private void buildAlignment() throws StructureException {

		// Initialize a new multiple alignment
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.getEnsemble().setAtomArrays(new ArrayList<Atom[]>());
		msa.getEnsemble().setAlgorithmName(CeSymm.algorithmName);
		msa.getEnsemble().setVersion(CeSymm.version);
		msa.getEnsemble().setStructureIdentifiers(
				new ArrayList<StructureIdentifier>());

		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(new ArrayList<List<Integer>>());

		// Calculate the connected groups of the alignment graph
		ConnectivityInspector<Integer, DefaultEdge> inspector = 
				new ConnectivityInspector<Integer, DefaultEdge>(alignGraph);
		List<Set<Integer>> comps = inspector.connectedSets();
		List<ResidueGroup> groups = new ArrayList<ResidueGroup>(comps.size());
		for (Set<Integer> comp : comps)
			groups.add(new ResidueGroup(comp));

		// Calculate thr order of symmetry from levels
		int order = 1;
		for (MultipleAlignment m : levels)
			order *= m.size();
		for (int su = 0; su < order; su++)
			b.getAlignRes().add(new ArrayList<Integer>());

		// Construct the resulting MultipleAlignment from components
		for (ResidueGroup group : groups) {
			if (group.order() != order)
				continue;
			group.combineWith(b.getAlignRes());
		}
		
		for (int su = 0; su < order; su++) {
			Collections.sort(b.getAlignRes().get(su));
			// TODO Set the identifier to the true range of the repeat
			msa.getEnsemble().getStructureIdentifiers()
					.add(new PassthroughIdentifier("S" + (su + 1)));
			msa.getEnsemble().getAtomArrays().add(allAtoms);
		}
		
		CoreSuperimposer imposer = new CoreSuperimposer();
		imposer.superimpose(msa);
		MultipleAlignmentScorer.calculateScores(msa);
		
		result.setMultipleAlignment(msa);
		result.setRefined(true);
		result.setSymmOrder(order);
		result.setSymmLevels(levels.size());
	}

	/**
	 * The symmetry axes of each level are recovered after the symmetry analysis
	 * iterations have finished, using the stored MultipleAlignment at each
	 * symmetry level.
	 */
	private void recoverAxes() {

		SymmetryAxes axes = new SymmetryAxes();

		int size = result.getSymmOrder();
		int parents = 1;

		for (int m = 0; m < levels.size(); m++) {

			MultipleAlignment align = levels.get(m);
			Matrix4d axis = align.getBlockSet(0).getTransformations().get(1);

			int subsize = align.size();
			parents *= subsize;
			size /= subsize;

			List<Integer> subunitTransform = new ArrayList<Integer>();
			for (int i = 0; i < size * parents; i++) {
				subunitTransform.add(0);
			}

			List<List<Integer>> superpose = new ArrayList<List<Integer>>();
			superpose.add(new ArrayList<Integer>());
			superpose.add(new ArrayList<Integer>());

			for (int su = 0; su < subsize - 1; su++) {
				for (int s = 0; s < size; s++) {
					Integer subIndex1 = su * size + s;
					Integer subIndex2 = (su + 1) * size + s;
					superpose.get(0).add(subIndex1);
					superpose.get(1).add(subIndex2);
				}
			}

			for (int p = 0; p < parents; p++) {
				for (int s = 0; s < size; s++) {
					subunitTransform.set(p * size + s, p % subsize);
				}
			}
			axes.addAxis(axis, superpose, subunitTransform, subsize);
		}
		result.setAxes(axes);
	}

	/**
	 * Calculate the number of helix and strand SSE of a subunit.
	 * 
	 * @param atoms
	 *            Atom array of the subunit found
	 * @return int number of helix or strand SSE
	 */
	private static int countHelixStrandSSE(Atom[] atoms) {

		List<SecStrucElement> sses = SecStrucTools
				.getSecStrucElements(SymmetryTools.getGroups(atoms));
		int count = 0;
		for (SecStrucElement sse : sses) {
			SecStrucType t = sse.getType();
			if (t.isBetaStrand() || t.isHelixType()) {
				count++;
			}
		}
		return count;
	}

}
