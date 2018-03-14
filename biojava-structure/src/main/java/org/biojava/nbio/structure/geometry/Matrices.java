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

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.jama.Matrix;

/**
 * Matrices contains static methods to operate and transform matrices used in 3D
 * geometry (transformation matrices and rotation matrices).
 * <p>
 * This class complements and extends the functionallity of vecmath and JAMA.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class Matrices {
	
	/** Prevent instantiation */
	private Matrices(){}

	/**
	 * Convert a transformation matrix into a JAMA rotation matrix. Because the
	 * JAMA matrix is a pre-multiplication matrix and the Vecmath matrix is a
	 * post-multiplication one, the rotation matrix is transposed to ensure that
	 * the transformation they produce is the same.
	 *
	 * @param transform
	 *            Matrix4d with transposed rotation matrix
	 * @return rotation matrix as JAMA object
	 */
	public static Matrix getRotationJAMA(Matrix4d transform) {

		Matrix rot = new Matrix(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				rot.set(j, i, transform.getElement(i, j)); // transposed
			}
		}
		return rot;
	}
	
	/**
	 * Convert a transformation matrix into a rotation matrix.
	 *
	 * @param transform
	 *            Matrix4d
	 * @return rotation matrix
	 */
	public static Matrix3d getRotationMatrix(Matrix4d transform) {

		Matrix3d rot = new Matrix3d();
		transform.setRotationScale(rot);
		return rot;
	}

	/**
	 * Extract the translational vector of a transformation matrix.
	 *
	 * @param transform
	 *            Matrix4d
	 * @return Vector3d translation vector
	 */
	public static Vector3d getTranslationVector(Matrix4d transform) {
		Vector3d transl = new Vector3d();
		transform.get(transl);
		return transl;
	}
	
	/**
	 * Convert JAMA rotation and translation to a Vecmath transformation matrix.
	 * Because the JAMA matrix is a pre-multiplication matrix and the Vecmath
	 * matrix is a post-multiplication one, the rotation matrix is transposed to
	 * ensure that the transformation they produce is the same.
	 *
	 * @param rot
	 *            3x3 Rotation matrix
	 * @param trans
	 *            3x1 Translation matrix
	 * @return 4x4 transformation matrix
	 */
	public static Matrix4d getTransformation(Matrix rot, Matrix trans) {
		return new Matrix4d(new Matrix3d(rot.getColumnPackedCopy()),
				new Vector3d(trans.getColumnPackedCopy()), 1.0);
	}

}
