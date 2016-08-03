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

import org.biojava.nbio.structure.geometry.CalcPoint;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import java.util.ArrayList;
import java.util.List;

public class HelixExtender {
	private QuatSymmetrySubunits subunits = null;
	private Helix helix = null;

	public HelixExtender(QuatSymmetrySubunits subunits, Helix helix) {
		this.subunits = subunits;
		this.helix = helix;
	}

	public Point3d[] extendHelix(int steps) {
		List<List<Integer>> layerLines = helix.getLayerLines();

		// get list of subunit indices to be used for helix extension
		List<Integer> indices = new ArrayList<Integer>();
		for (List<Integer> line: layerLines) {
			if (steps < 0) {
				indices.add(line.get(0));
			} else if (steps > 0) {
				indices.add(line.get(line.size()-1));
			}
		}
		System.out.println("Extending subunits: " + indices);

		List<Point3d> points = new ArrayList<Point3d>();
		Matrix4d transformation = helix.getTransformation();
		for (int index: indices) {
	    	Point3d[] trace = subunits.getTraces().get(index);
	    	Point3d[] copy = CalcPoint.clonePoint3dArray(trace);
		    for (int i = 0; i < Math.abs(steps); i++) {
		    	CalcPoint.transform(transformation, copy);
		    }
		    for (Point3d p: copy) {
		    	points.add(p);
		    }
		}
		return points.toArray(new Point3d[0]);
	}

}
