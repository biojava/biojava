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
package org.biojava.nbio.structure.symmetry.core;

import javax.vecmath.Point3d;

import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class SubunitGraph {

	private static final Logger logger = LoggerFactory
			.getLogger(SubunitGraph.class);

	private static final double DISTANCE_CUTOFF = 8;
	private static final int MIN_CONTACTS = 10;
	private List<Point3d[]> caCoords = null;

	public SubunitGraph(List<Point3d[]> caCoords) {
		this.caCoords = caCoords;
	}

	public UndirectedGraph<Integer, DefaultEdge> getProteinGraph() {
		int n = caCoords.size();

		// add vertex for each chain center
		UndirectedGraph<Integer, DefaultEdge> graph = new SimpleGraph<Integer, DefaultEdge>(DefaultEdge.class);
		for (int i = 0; i < n; i++) {
			graph.addVertex(i);
		}
		
		// add edges if there are 10 or more contact of Calpha atoms
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				int numContacts = calcContactNumber(caCoords.get(i), caCoords.get(j));
				logger.debug("Calpha contacts between subunits {},{}: {}", i, j, numContacts);
				if (numContacts >= MIN_CONTACTS) {
					graph.addEdge(i, j);
				}
			}
		}

		return graph;
	}

	private int calcContactNumber(Point3d[] a, Point3d[] b) {
		double distCutoffSq = DISTANCE_CUTOFF * DISTANCE_CUTOFF;
		int contacts = 0;
		for (Point3d pa : a) {
			for (Point3d pb : b) {
				if (pa.distanceSquared(pb) < distCutoffSq) {
					contacts++;
				}
			}
		}
		return contacts;
	}
}
