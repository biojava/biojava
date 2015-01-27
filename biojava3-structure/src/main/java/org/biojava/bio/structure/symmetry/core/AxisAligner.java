package org.biojava.bio.structure.symmetry.core;

import java.util.List;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public abstract class AxisAligner {
	
	/**
	 * Returns an instance of AxisAligner for differnt type of QuatSymmetryResults (factory method)
	 * @param results symmetry results
	 * @return instance of AxisAligner
	 */
	public static AxisAligner getInstance(QuatSymmetryResults results) {
		String symmetry = results.getSymmetry();
		
		if (symmetry.equals("H")) {
			return new HelixAxisAligner(results);
		} else {
			return new RotationAxisAligner(results);
		}
	}

	public abstract String getSymmetry();
	
	public abstract Matrix4d getTransformation();

	public abstract Matrix3d getRotationMatrix();

	public abstract Matrix4d getReverseTransformation();

	public abstract Vector3d getPrincipalRotationAxis();

	public abstract Vector3d getRotationReferenceAxis();

	public abstract Vector3d[] getPrincipalAxesOfInertia();

	public abstract Vector3d getDimension();

	/**
	 * Returns the radius for drawing polyhedra
	 * @return double radius
	 */
	public abstract double getRadius();

	/**
	 * Returns a transformation matrix transform polyhedra for Cn structures.
	 * The center in this matrix is the geometric center, rather then the centroid.
	 * In Cn structures those are usually not the same.
	 * @return
	 */
	public abstract Matrix4d getGeometicCenterTransformation();

	/**
	 * Returns the geometric center of polyhedron. In the case of the Cn 
	 * point group, the centroid and geometric center are usually not
	 * identical.
	 * @return
	 */
	public abstract Point3d getGeometricCenter();

	public abstract Point3d getCentroid();

	public abstract Subunits getSubunits();

	public abstract List<List<Integer>> getOrbits();

}