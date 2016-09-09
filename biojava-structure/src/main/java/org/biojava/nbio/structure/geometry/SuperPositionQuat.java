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

import javax.vecmath.*;

/**
 * The SuperPositionQuat implements a quaternion based algorithm to superpose
 * arrays of Points in 3D.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public final class SuperPositionQuat extends SuperPositionAbstract {

	/**
	 * Constructor for the quaternion based superposition algorithm.
	 * 
	 * @param centered
	 *            true if the point arrays are centered at the origin (faster),
	 *            false otherwise
	 */
	public SuperPositionQuat(boolean centered) {
		super(centered);
	}

	@Override
	public Matrix4d superpose(Point3d[] fixed, Point3d[] moved) {

		checkInput(fixed, moved);
		
		if (centered) {
			Quat4d q = UnitQuaternions.relativeOrientation(fixed, moved);
			Matrix4d rotTrans = new Matrix4d();
			rotTrans.set(q);
			return rotTrans;
		}

		// translate to origin
		Point3d[] xref = CalcPoint.clonePoint3dArray(fixed);
		Point3d xtrans = CalcPoint.centroid(xref);
		xtrans.negate();
		CalcPoint.translate(xtrans, xref);

		Point3d[] yref = CalcPoint.clonePoint3dArray(moved);
		Point3d ytrans = CalcPoint.centroid(yref);
		ytrans.negate();
		CalcPoint.translate(ytrans, yref);

		// calculate rotational component (rotation around origin)
		Quat4d q = UnitQuaternions.relativeOrientation(xref, yref);
		Matrix4d rotTrans = new Matrix4d();
		rotTrans.set(q);

		// combine with moved -> origin translation
		Matrix4d trans = new Matrix4d();
		trans.setIdentity();
		trans.setTranslation(new Vector3d(ytrans));
		rotTrans.mul(rotTrans, trans);

		// combine with origin -> fixed translation
		xtrans.negate();
		Matrix4d transInverse = new Matrix4d();
		transInverse.setIdentity();
		transInverse.setTranslation(new Vector3d(xtrans));
		rotTrans.mul(transInverse, rotTrans);

		return rotTrans;

	}

	@Deprecated
	public static Matrix4d superposeAtOrigin(Point3d[] x, Point3d[] y,
			AxisAngle4d axisAngle) {
		Quat4d q = UnitQuaternions.relativeOrientation(x, y);
		Matrix4d rotTrans = new Matrix4d();
		rotTrans.setIdentity();
		rotTrans.set(q);
		axisAngle.set(q);
		Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		if (axis.lengthSquared() < 1.0E-6) {
			// System.err.println("Error: SuperPosition.superposeAtOrigin: axis vector undefined!");
			axisAngle.x = 0;
			axisAngle.y = 0;
			axisAngle.z = 1;
			axisAngle.angle = 0;
		} else {
			axis.normalize();
			axisAngle.x = axis.x;
			axisAngle.y = axis.y;
			axisAngle.z = axis.z;
		}
		CalcPoint.transform(rotTrans, x);

		return rotTrans;
	}

	@Override
	public double getRmsd(Point3d[] x, Point3d[] y) {
		Point3d[] yclone = CalcPoint.clonePoint3dArray(y);
		superposeAndTransform(x, yclone);
		return CalcPoint.rmsd(x, yclone);
	}

}
