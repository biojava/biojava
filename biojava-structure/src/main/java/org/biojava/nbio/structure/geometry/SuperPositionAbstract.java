package org.biojava.nbio.structure.geometry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

/**
 * The SuperPositionAbstract contains common code shared by all SuperPosition
 * algorithm implementations.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public abstract class SuperPositionAbstract implements SuperPosition {

	protected boolean centered;

	public SuperPositionAbstract(boolean centered) {
		this.centered = centered;
	}

	@Override
	public Matrix4d superposeAndTransform(Point3d[] fixed, Point3d[] moved) {
		Matrix4d rotTrans = superpose(fixed, moved);
		CalcPoint.transform(rotTrans, moved);
		return rotTrans;
	}

	/**
	 * Check that the input to the superposition algorithms is valid.
	 * 
	 * @param fixed
	 * @param moved
	 */
	protected void checkInput(Point3d[] fixed, Point3d[] moved) {
		if (fixed.length != moved.length)
			throw new IllegalArgumentException(
					"Point arrays to superpose are of different lengths.");
	}

	/**
	 * @param centered
	 *            true if the point arrays are centered at the origin (faster),
	 *            false otherwise
	 */
	public void setCentered(boolean centered) {
		this.centered = centered;
	}

}
