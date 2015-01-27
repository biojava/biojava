/**
 * 
 */
package org.biojava.bio.structure.symmetry.geometry;

import java.util.List;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;

/**
 * @author Peter
 *
 */
public interface Polyhedron {

	public Point3d[] getVertices();
	public List<int[]> getLineLoops();
	public double getCirumscribedRadius();
	public int getViewCount();
	public String getViewName(int index);
	public Matrix3d getViewMatrix(int index);
}
