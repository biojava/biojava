package org.biojava.bio.structure.symmetry.core;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.bio.structure.symmetry.geometry.SuperPosition;

public class HelixExtender {
	private Subunits subunits = null;
	private Helix helix = null;
	
	public HelixExtender(Subunits subunits, Helix helix) {
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
	    	Point3d[] copy = SuperPosition.clonePoint3dArray(trace);
		    for (int i = 0; i < Math.abs(steps); i++) {
		    	SuperPosition.transform(transformation, copy);
		    }
		    for (Point3d p: copy) {
		    	points.add(p);
		    }
		}
		return points.toArray(new Point3d[0]);
	}
	
}
