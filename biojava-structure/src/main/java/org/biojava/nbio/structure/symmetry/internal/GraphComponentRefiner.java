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
import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.vecmath.GMatrix;
import javax.vecmath.GVector;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;

/**
 * The GraphRefiner transforms the self-alignment into a Graph and extracts its
 * maximally connected Components. It then refines the alignment by combining
 * the compatible Components with the following heuristic:
 *
 * <pre>
 * Given a set of components and their pairwise compatibilities, iteratively
 * add the most compatible component, which is compatible to all the components
 * already added, to the final alignment.
 * </pre>
 *
 * @author Aleix Lafita
 *
 */
public class GraphComponentRefiner implements SymmetryRefiner {

	@Override
	public MultipleAlignment refine(AFPChain selfAlignment, Atom[] atoms, int order)
			throws StructureException, RefinerFailedException {

		// Construct the alignment graph with jgrapht
		Graph<Integer, DefaultEdge> graph = SymmetryTools
				.buildSymmetryGraph(selfAlignment);

		// Find the maximally connected components of the graph
		ConnectivityInspector<Integer, DefaultEdge> inspector = new ConnectivityInspector<Integer, DefaultEdge>(
				graph);
		List<Set<Integer>> components = inspector.connectedSets();

		// Filter components with size != order, and transform to ResidueGroups
		List<ResidueGroup> groups = new ArrayList<ResidueGroup>();
		for (Set<Integer> comp : components) {
			if (comp.size() == order) {
				ResidueGroup group = new ResidueGroup(comp);
				groups.add(group);
			}
		}
		int size = groups.size();
		if (size == 0)
			throw new RefinerFailedException("Could not find any components "
					+ "in the alignment Graph of size != order.");

		// Create a square matrix of component compatibility
		GMatrix matrix = new GMatrix(size, size);
		for (int i = 0; i < size; i++) {
			for (int j = i; j < size; j++) {
				// The diagonal is always 0
				if (i == j){
					matrix.setElement(i, j, 0);
					continue;
				}
				// If compatible put 1, otherwise 0
				ResidueGroup g1 = groups.get(i);
				ResidueGroup g2 = groups.get(j);
				if (g1.isCompatible(g2)) {
					matrix.setElement(i, j, 1);
					matrix.setElement(j, i, 1);
				} else {
					matrix.setElement(i, j, 0);
					matrix.setElement(j, i, 0);
				}
			}
		}

		// The compatibility score is the sum of rows of the matrix
		List<Integer> rowScores = new ArrayList<Integer>(size);
		for (int i = 0; i < size; i++) {
			GVector row = new GVector(size);
			matrix.getRow(i, row);
			// because element={0,1}, this is the sum
			int rowScore = (int) row.normSquared();
			rowScores.add(rowScore);
		}

		// Refined multiple alignment Block as a result
		List<List<Integer>> alignRes = new ArrayList<List<Integer>>(order);
		for (int i = 0; i < order; i++)
			alignRes.add(new ArrayList<Integer>());

		// Iterate until no more groups left to add (all groups score 0)
		while (true) {
			// Take the most compatible ResidueGroup
			Integer max = Collections.max(rowScores);
			int index = rowScores.indexOf(max);

			// Add the group to the alignment Block
			groups.get(index).combineWith(alignRes);

			// Zero all the scores of incompatible groups
			boolean allZero = true;
			for (int i=0; i<size; i++){
				if (matrix.getElement(index, i) < 1.0)
					rowScores.set(i, 0);
				else if (rowScores.get(i) != 0)
					allZero = false;
			}
			if (allZero)
				break;
		}

		for (int i = 0; i < order; i++)
			Collections.sort(alignRes.get(i));

		int length = alignRes.get(0).size();
		if (length == 0)
			throw new RefinerFailedException("Empty alignment");

		int[][][] optAln = new int[order][2][length];
		for (int bk = 0; bk < order; bk++) {
			optAln[bk] = new int[2][];
			optAln[bk][0] = new int[length];
			optAln[bk][1] = new int[length];
			for (int pos = 0; pos < length; pos++) {
				optAln[bk][0][pos] = alignRes.get(bk).get(pos);
				optAln[bk][1][pos] = alignRes.get((bk + 1) % order).get(pos);
			}
		}
		AFPChain afp = AlignmentTools.replaceOptAln(optAln, selfAlignment, atoms, atoms);
		return SymmetryTools.fromAFP(afp, atoms);
	}

}
