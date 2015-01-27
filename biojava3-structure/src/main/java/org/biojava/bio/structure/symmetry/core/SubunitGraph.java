package org.biojava.bio.structure.symmetry.core;

import java.util.List;

import javax.vecmath.Point3d;

import org.biojava.bio.structure.symmetry.utils.Graph;
import org.biojava.bio.structure.symmetry.utils.SimpleGraph;

public class SubunitGraph {
	private static double DISTANCE_CUTOFF = 8;
	private static int MIN_CONTACTS = 10;
	private List<Point3d[]> caCoords = null;

	public SubunitGraph(List<Point3d[]> caCoords) {
		this.caCoords = caCoords;
	}

	public Graph<Integer> getProteinGraph() {
		int n = caCoords.size();

		// add vertex for each chain center
		Graph<Integer> graph = new SimpleGraph<Integer>(); 
		for (int i = 0; i < n; i++) {
			graph.addVertex(i);
		}

		// add edges if there are 10 or more contact of Calpha atoms
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
//				System.out.println("contacts: " + calcContactNumber(caCoords.get(i), caCoords.get(j)));
				if (calcContactNumber(caCoords.get(i), caCoords.get(j)) >= MIN_CONTACTS) {
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
