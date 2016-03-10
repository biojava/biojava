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
/**
 *
 */
package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.List;


/**
 * @author Peter
 *
 */
public class Prism implements Polyhedron {
	private int n = 3;
	private double circumscribedRadius = 1.0;
	private double height = 1.0;

	public Prism(int n) {
		this.n = n;
	}

	/**
	 * @return the height
	 */
	public double getHeight() {
		return height;
	}

	/**
	 * @param height the height to set
	 */
	public void setHeight(double height) {
		this.height = height;
	}

	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	@Override
	public double getCirumscribedRadius() {
		return circumscribedRadius;
	}

	/**
	 * Set the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @param cirumscribedRadius the cirumscribedRadius to set
	 */
	public void setCirumscribedRadius(double cirumscribedRadius) {
		this.circumscribedRadius = cirumscribedRadius;
	}
	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each
	 * of the icosahedron's faces
	 * @return the inscribedRadius
	 */
	public double getInscribedRadius() {
		double side = getSideLengthFromCircumscribedRadius(circumscribedRadius, n);
		return getInscribedRadiusFromSideLength(side, n);
	}

	/**
	 * Sets the radius of an inscribed sphere, that is tangent to each
	 * of the icosahedron's faces
	 * @param inscribedRadius the inscribedRadius to set
	 */
	public void setInscribedRadius(double radius) {
		double side = getSideLengthFromInscribedRadius(radius, n);
		this.circumscribedRadius = getCircumscribedRadiusFromSideLength(side, n);
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center
	 * @return
	 */
	@Override
	public Point3d[] getVertices() {
		Point3d[] polygon = new Point3d[2*n];
		Matrix3d m = new Matrix3d();

		Point3d center = new Point3d(0, 0, height/2);

		for (int i = 0; i < n; i++) {
			polygon[i] = new Point3d(0, circumscribedRadius, 0);
			m.rotZ(i*2*Math.PI/n);
			m.transform(polygon[i]);
			polygon[n+i] = new  Point3d(polygon[i]);
			polygon[i].sub(center);
			polygon[n+i].add(center);
		}

		return polygon;
	};

	@Override
	public List<int[]> getLineLoops() {
		List<int[]> list = new ArrayList<int[]>();
		int[] l1 = new int[2*n+2];
		for (int i = 0; i < n; i++) {
			l1[i] = i;
		}
		l1[n] = 0;
		for (int i = 0; i < n; i++) {
			l1[n+i+1] = n+i;
		}
		l1[2*n+1] = l1[n+1];
		list.add(l1);

		for (int i = 1; i < n; i++) {
			int[] l2 = new int[2];
			l2[0] = i;
			l2[1] = n+i;
			list.add(l2);
		}
		return list;
	}

	/**
	 * Returns the vertices of an n-fold polygon of given radius and center
	 * @return
	 */
	public static Point3d[] getPolygonVertices(int n, double radius, Point3d center) {
		Point3d[] polygon = new Point3d[n];
		Matrix3d m = new Matrix3d();

		for (int i = 0; i < n; i++) {
			polygon[i] = new Point3d(0, radius, 0);
			m.rotZ(i*2*Math.PI/n);
			m.transform(polygon[i]);
			polygon[i].add(center);
		}
		return polygon;
	}

	@Override
	public int getViewCount() {
		return 4;
	}

	@Override
	public String getViewName(int index) {
		String name;
		switch (index) {
		case 0:  name = "Front C" + n + " axis";
		break;
		case 1:  name = "Side edge-centered";
		break;
		case 2:  name = "Side face-centered";
		break;
		case 3:  name = "Back C" + n + " axis";
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
			m.setIdentity(); // front
			break;
		case 1:
			m.rotX(Math.PI/2); // side edge-centered
			break;
		case 2:
			m.rotY(Math.PI/n); // side face-centered
			Matrix3d m1 = new Matrix3d();
			m1.rotX(Math.PI/2);
			m.mul(m1);
			break;
		case 3:
			m.set(flipX()); // back
			break;
		default:
			throw new IllegalArgumentException("getViewMatrix: index out of range:" + index);
		}
		return m;
	}

	// http://www.mathopenref.com/polygonincircle.html
	private static double getSideLengthFromInscribedRadius(double radius, int n) {
		if (n == 2) {
			return radius;
		}
		return radius * 2 * Math.tan(Math.PI/n);
	}

	private static double getInscribedRadiusFromSideLength(double length, int n) {
		if (n == 2) {
			return length;
		}
		return length / (2 * Math.tan(Math.PI/n));
	}

	// http://www.mathopenref.com/polygonradius.html
	private static double getSideLengthFromCircumscribedRadius(double radius, int n) {
		if (n == 2) {
			return radius;
		}
		return radius * (2 * Math.sin(Math.PI/n));
	}

	private static double getCircumscribedRadiusFromSideLength(double length, int n) {
		if (n == 2) {
			return length;
		}
		return length / (2 * Math.sin(Math.PI/n));
	}

	private static Matrix3d flipX() {
		Matrix3d rot = new Matrix3d();
		rot.m00 = 1;
		rot.m11 = -1;
		rot.m22 = -1;
		return rot;
	}

}
