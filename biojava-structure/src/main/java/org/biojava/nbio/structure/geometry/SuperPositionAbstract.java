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
