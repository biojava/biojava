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
package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.List;


public class Octahedron implements Polyhedron {
	private static double TETRAHEDRAL_ANGLE = Math.acos(-1.0/3.0);
	private static int[] lineLoop1 = {2,4,3,5,2,1,3,0,5,1,4,0,2};
	private double cirumscribedRadius = 1.0;

	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	@Override
	public double getCirumscribedRadius() {
		return cirumscribedRadius;
	}

	/**
	 * Set the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @param cirumscribedRadius the cirumscribedRadius to set
	 */
	public void setCirumscribedRadius(double cirumscribedRadius) {
		this.cirumscribedRadius = cirumscribedRadius;
	}
	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each
	 * of the octahedron's faces
	 * @return the inscribedRadius
	 */
	public double getInscribedRadius() {
		double side = getSideLengthFromCircumscribedRadius(cirumscribedRadius);
		return getInscribedRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of an inscribed sphere, that is tangent to each
	 * of the octahedron's faces
	 * @param inscribedRadius the inscribedRadius to set
	 */
	public void setInscribedRadius(double radius) {
		double side = getSideLengthFromInscribedRadius(radius);
		this.cirumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the radius of a sphere, that is tangent to each
	 * of the octahedron's edges
	 *
	 * @return the midRadius
	 */
	public double getMidRadius() {
		double side = getSideLengthFromCircumscribedRadius(cirumscribedRadius);
		return getMiddleRadiusFromSideLength(side);
	}

	/**
	 * Sets the radius of radius of a sphere, that is tangent to each
	 * of the octahedron's edges
	 * @param midRadius the midRadius to set
	 */
	public void setMidRadius(double radius) {
		double side = getSideLengthFromMiddleRadius(radius);
		this.cirumscribedRadius = getCircumscribedRadiusFromSideLength(side);
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center
	 * @param n
	 * @param radius
	 * @param center
	 * @return
	 */
	@Override
	public Point3d[] getVertices() {
		Point3d[] octahedron = new Point3d[6];
	    octahedron[0] = new Point3d(-cirumscribedRadius, 0, 0);
	    octahedron[1] = new Point3d( cirumscribedRadius, 0, 0);
	    octahedron[2] = new Point3d(0, -cirumscribedRadius, 0);
	    octahedron[3] = new Point3d(0,  cirumscribedRadius, 0);
	    octahedron[4] = new Point3d(0, 0, -cirumscribedRadius);
	    octahedron[5] = new Point3d(0, 0,  cirumscribedRadius);

		return octahedron;
	};

	@Override
	public List<int[]> getLineLoops() {
		return Arrays.asList(lineLoop1);
	}

	public Point3d getC4Axis(double scale) {
		return new Point3d(0, 0, cirumscribedRadius*scale);
	}
	public Point3d getC3Axis(double scale) {
		double s = 1/Math.sqrt(1 + Math.sqrt(2));
		return new Point3d(cirumscribedRadius*scale*s, cirumscribedRadius*scale*s, cirumscribedRadius*scale*s);
	}
	public Point3d getC2Axis(double scale) {
		double s = 1/Math.sqrt(2);
		return new Point3d(cirumscribedRadius*scale*s, cirumscribedRadius*scale*s, 0);
	}

	@Override
	public int getViewCount() {
		return 3;
	}

	@Override
	public String getViewName(int index) {
		String name;
		switch (index) {
		case 0:  name = "4-fold axis vertex-centered";
		break;
		case 1:  name = "3-fold axis face-centered";
		break;
		case 2:  name = "2-fold axis edge-centered";
		break;
		default: throw new IllegalArgumentException("getViewMatrix: index out of range:" + index);
		}
		return name;
	}

	@Override
	public Matrix3d getViewMatrix(int index) {
		Matrix3d m = new Matrix3d();
		switch (index) {
		case 0:
			m.setIdentity(); // C4 vertex-centered
			break;
		case 1:
			m.rotX(-0.5 * TETRAHEDRAL_ANGLE); // C3 face-centered  2.0*Math.PI/3
			Matrix3d m1 = new Matrix3d();
			m1.rotZ(Math.PI/4);
			m.mul(m1);
			break;
		case 2:
			m.rotY(Math.PI/4); // side face-centered
			break;
		default:
			throw new IllegalArgumentException("getViewMatrix: index out of range:" + index);
		}
		return m;
	}

	private static double getSideLengthFromInscribedRadius(double radius) {
		return radius * 6 / Math.sqrt(6);
	}

	private static double getInscribedRadiusFromSideLength(double sideLength) {
		return sideLength / 6 * Math.sqrt(6);
	}

	private static double getSideLengthFromMiddleRadius(double radius) {
		return radius * 2;
	}

	private static double getMiddleRadiusFromSideLength(double sideLength) {
		return sideLength / 2;
	}

	private static double getSideLengthFromCircumscribedRadius(double radius) {
		return radius * 2 / Math.sqrt(2);
	}

	private static double getCircumscribedRadiusFromSideLength(double sideLength) {
		return sideLength / 2 * Math.sqrt(2);
	}
}
