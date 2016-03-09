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


public class RectangularPrism implements Polyhedron {
	private static int[] lineLoop1 = {0,1,2,3,0,4,5,6,7,4};
	private static int[] lineLoop2 = {1,5};
	private static int[] lineLoop3 = {2,6};
	private static int[] lineLoop4 = {3,7};
	private double length = 1.0;
	private double width = 1.0;
	private double height = 1.0;
	private static String[] viewNames = {"Front", "Left", "Back", "Right", "Top", "Bottom"};

	public RectangularPrism(double length, double width, double height) {
		this.length = length;
		this.width = width;
		this.height = height;
	}

	/**
	 * Returns the radius of a circumscribed sphere, that goes
	 * through all vertices
	 * @return the cirumscribedRadius
	 */
	public double getLength() {
		return length;
	}

	/**
	 * Returns the radius of an inscribed sphere, that is tangent to each
	 * of the octahedron's faces
	 * @return the inscribedRadius
	 */
	public double getWidth() {
		return width;
	}

	/**
	 * Returns the radius of a sphere, that is tangent to each
	 * of the octahedron's edges
	 *
	 * @return the midRadius
	 */
	public double getHeight() {
		return height;
	}

	/**
	 * Returns the radius of a circumscribed sphere (length of diagonal of
	 * rectangular prism/2, that goes through at least four vertices
	 * @return the cirumscribedRadius
	 */
	@Override
	public double getCirumscribedRadius() {
		return 0.5* Math.sqrt(width*width + height*height + length*length);
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
		double x = 0.5 * width;
		double y = 0.5 * height;
		double z = 0.5 * length;
		Point3d[] vertices = new Point3d[8];
		vertices[0] = new Point3d(-x, -y,  z);
		vertices[1] = new Point3d(-x,  y,  z);
		vertices[2] = new Point3d( x,  y,  z);
		vertices[3] = new Point3d( x, -y,  z);
		vertices[4] = new Point3d(-x, -y, -z);
		vertices[5] = new Point3d(-x,  y, -z);
		vertices[6] = new Point3d( x,  y, -z);
		vertices[7] = new Point3d( x, -y, -z);

		return vertices;
	};

	@Override
	public List<int[]> getLineLoops() {
		return Arrays.asList(lineLoop1, lineLoop2, lineLoop3, lineLoop4);
	}

	@Override
	public int getViewCount() {
		return viewNames.length;
	}

	@Override
	public String getViewName(int index) {
		return viewNames[index];
	}

	@Override
	public Matrix3d getViewMatrix(int index) {
		Matrix3d m = new Matrix3d();
		switch (index) {
		case 0:  m.setIdentity(); // front
		break;
		case 1:  m.rotY(Math.PI/2); // left
		break;
		case 2:  m.rotY(Math.PI); // back
		break;
		case 3:  m.rotY(-Math.PI/2); // right
		break;
		case 4:  m.rotX(Math.PI/2); // top
		break;
		case 5:  m.rotX(-Math.PI/2); // bottom
		break;
		default: throw new IllegalArgumentException("getViewMatrix: index out of range:" + index);
		}
		return m;
	}
}
