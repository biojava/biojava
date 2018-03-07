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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;

/**
 * The GraphOrderDetector transforms the self-alignment into a Graph and
 * extracts its maximally connected Components.
 * <p>
 * The order reported is the one that maximizes the number of residues aligned,
 * i.e. the highest order (component size) times the frequency of the Component.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class GraphComponentOrderDetector implements OrderDetector {

	@Override
	public int calculateOrder(AFPChain selfAlignment, Atom[] ca)
			throws RefinerFailedException {

		// Construct the alignment graph with jgrapht
		Graph<Integer, DefaultEdge> graph = SymmetryTools
				.buildSymmetryGraph(selfAlignment);

		// Find the maximally connected components of the graph
		ConnectivityInspector<Integer, DefaultEdge> inspector =
				new ConnectivityInspector<Integer, DefaultEdge>(graph);
		List<Set<Integer>> components = inspector.connectedSets();

		// The order maximizes the residues aligned
		Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
		for (Set<Integer> c : components) {
			if (counts.containsKey(c.size()))
				counts.put(c.size(), counts.get(c.size()) + c.size());
			else
				counts.put(c.size(), c.size());
		}
		int maxCounts = 0;
		int order = 1;
		for (Integer ord : counts.keySet()) {
			if (counts.get(ord) > maxCounts) {
				order = ord;
				maxCounts = counts.get(ord);
			}
		}
		return order;
	}

}
