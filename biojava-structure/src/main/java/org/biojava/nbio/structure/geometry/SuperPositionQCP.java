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
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implementation of the Quaternion-Based Characteristic Polynomial algorithm
 * for RMSD and Superposition calculations.
 * <p>
 * Usage:
 * <p>
 * The input consists of 2 Point3d arrays of equal length. The input coordinates
 * are not changed.
 * 
 * <pre>
 *    Point3d[] x = ...
 *    Point3d[] y = ...
 *    SuperPositionQCP qcp = new SuperPositionQCP();
 *    qcp.set(x, y);
 * </pre>
 * <p>
 * or with weighting factors [0 - 1]]
 * 
 * <pre>
 *    double[] weights = ...
 *    qcp.set(x, y, weights);
 * </pre>
 * <p>
 * For maximum efficiency, create a SuperPositionQCP object once and reuse it.
 * <p>
 * A. Calculate rmsd only
 * 
 * <pre>
 * double rmsd = qcp.getRmsd();
 * </pre>
 * <p>
 * B. Calculate a 4x4 transformation (rotation and translation) matrix
 * 
 * <pre>
 * Matrix4d rottrans = qcp.getTransformationMatrix();
 * </pre>
 * <p>
 * C. Get transformated points (y superposed onto the reference x)
 * 
 * <pre>
 * Point3d[] ySuperposed = qcp.getTransformedCoordinates();
 * </pre>
 * <p>
 * Citations:
 * <p>
 * Liu P, Agrafiotis DK, & Theobald DL (2011) Reply to comment on: "Fast
 * determination of the optimal rotation matrix for macromolecular
 * superpositions." Journal of Computational Chemistry 32(1):185-186.
 * [http://dx.doi.org/10.1002/jcc.21606]
 * <p>
 * Liu P, Agrafiotis DK, & Theobald DL (2010) "Fast determination of the optimal
 * rotation matrix for macromolecular superpositions." Journal of Computational
 * Chemistry 31(7):1561-1563. [http://dx.doi.org/10.1002/jcc.21439]
 * <p>
 * Douglas L Theobald (2005) "Rapid calculation of RMSDs using a
 * quaternion-based characteristic polynomial." Acta Crystallogr A
 * 61(4):478-480. [http://dx.doi.org/10.1107/S0108767305015266 ]
 * <p>
 * This is an adoption of the original C code QCProt 1.4 (2012, October 10) to
 * Java. The original C source code is available from
 * http://theobald.brandeis.edu/qcp/ and was developed by
 * <p>
 * Douglas L. Theobald Department of Biochemistry MS 009 Brandeis University 415
 * South St Waltham, MA 02453 USA
 * <p>
 * dtheobald@brandeis.edu
 * <p>
 * Pu Liu Johnson & Johnson Pharmaceutical Research and Development, L.L.C. 665
 * Stockton Drive Exton, PA 19341 USA
 * <p>
 * pliu24@its.jnj.com
 * <p>
 * 
 * @author Douglas L. Theobald (original C code)
 * @author Pu Liu (original C code)
 * @author Peter Rose (adopted to Java)
 * @author Aleix Lafita (adopted to Java)
 */
public final class SuperPositionQCP extends SuperPositionAbstract {

	private static final Logger logger = LoggerFactory.getLogger(SuperPositionQCP.class);

	private double evec_prec = 1E-6;
	private double eval_prec = 1E-11;

	private Point3d[] x;
	private Point3d[] y;

	private double[] weight;
	private double wsum;

	private Point3d[] xref;
	private Point3d[] yref;
	private Point3d xtrans;
	private Point3d ytrans;

	private double e0;
	private Matrix3d rotmat = new Matrix3d();
	private Matrix4d transformation = new Matrix4d();
	private double rmsd = 0;
	private double Sxy, Sxz, Syx, Syz, Szx, Szy;
	private double SxxpSyy, Szz, mxEigenV, SyzmSzy, SxzmSzx, SxymSyx;
	private double SxxmSyy, SxypSyx, SxzpSzx;
	private double Syy, Sxx, SyzpSzy;
	private boolean rmsdCalculated = false;
	private boolean transformationCalculated = false;
	private boolean centered = false;

	/**
	 * Default constructor for the quaternion based superposition algorithm.
	 * 
	 * @param centered
	 *            true if the point arrays are centered at the origin (faster),
	 *            false otherwise
	 */
	public SuperPositionQCP(boolean centered) {
		super(centered);
	}

	/**
	 * Constructor with option to set the precision values.
	 * 
	 * @param centered
	 *            true if the point arrays are centered at the origin (faster),
	 *            false otherwise
	 * @param evec_prec
	 *            required eigenvector precision
	 * @param eval_prec
	 *            required eigenvalue precision
	 */
	public SuperPositionQCP(boolean centered, double evec_prec, double eval_prec) {
		super(centered);
		this.evec_prec = evec_prec;
		this.eval_prec = eval_prec;
	}

	/**
	 * Sets the two input coordinate arrays. These input arrays must be of equal
	 * length. Input coordinates are not modified.
	 * 
	 * @param x
	 *            3d points of reference coordinate set
	 * @param y
	 *            3d points of coordinate set for superposition
	 */
	private void set(Point3d[] x, Point3d[] y) {
		this.x = x;
		this.y = y;
		rmsdCalculated = false;
		transformationCalculated = false;
	}

	/**
	 * Sets the two input coordinate arrays and weight array. All input arrays
	 * must be of equal length. Input coordinates are not modified.
	 * 
	 * @param x
	 *            3d points of reference coordinate set
	 * @param y
	 *            3d points of coordinate set for superposition
	 * @param weight
	 *            a weight in the inclusive range [0,1] for each point
	 */
	private void set(Point3d[] x, Point3d[] y, double[] weight) {
		this.x = x;
		this.y = y;
		this.weight = weight;
		rmsdCalculated = false;
		transformationCalculated = false;
	}

	/**
	 * Return the RMSD of the superposition of input coordinate set y onto x.
	 * Note, this is the fasted way to calculate an RMSD without actually
	 * superposing the two sets. The calculation is performed "lazy", meaning
	 * calculations are only performed if necessary.
	 * 
	 * @return root mean square deviation for superposition of y onto x
	 */
	private double getRmsd() {
		if (!rmsdCalculated) {
			calcRmsd(x, y);
			rmsdCalculated = true;
		}
		return rmsd;
	}

	/**
	 * Weighted superposition.
	 * 
	 * @param fixed
	 * @param moved
	 * @param weight
	 *            array of weigths for each equivalent point position
	 * @return
	 */
	public Matrix4d weightedSuperpose(Point3d[] fixed, Point3d[] moved, double[] weight) {
		set(moved, fixed, weight);
		getRotationMatrix();
		if (!centered) {
			calcTransformation();
		} else {
			transformation.set(rotmat);
		}
		return transformation;
	}

	private Matrix3d getRotationMatrix() {
		getRmsd();
		if (!transformationCalculated) {
			calcRotationMatrix();
			transformationCalculated = true;
		}
		return rotmat;
	}

	/**
	 * Calculates the RMSD value for superposition of y onto x. This requires
	 * the coordinates to be precentered.
	 * 
	 * @param x
	 *            3d points of reference coordinate set
	 * @param y
	 *            3d points of coordinate set for superposition
	 */
	private void calcRmsd(Point3d[] x, Point3d[] y) {
		if (centered) {
			innerProduct(y, x);
		} else {
			// translate to origin
			xref = CalcPoint.clonePoint3dArray(x);
			xtrans = CalcPoint.centroid(xref);
			logger.debug("x centroid: " + xtrans);
			xtrans.negate();
			CalcPoint.translate(new Vector3d(xtrans), xref);

			yref = CalcPoint.clonePoint3dArray(y);
			ytrans = CalcPoint.centroid(yref);
			logger.debug("y centroid: " + ytrans);
			ytrans.negate();
			CalcPoint.translate(new Vector3d(ytrans), yref);
			innerProduct(yref, xref);
		}
		calcRmsd(wsum);
	}

	/**
	 * Superposition coords2 onto coords1 -- in other words, coords2 is rotated,
	 * coords1 is held fixed
	 */
	private void calcTransformation() {

		// transformation.set(rotmat,new Vector3d(0,0,0), 1);
		transformation.set(rotmat);
		// long t2 = System.nanoTime();
		// System.out.println("create transformation: " + (t2-t1));
		// System.out.println("m3d -> m4d");
		// System.out.println(transformation);

		// combine with x -> origin translation
		Matrix4d trans = new Matrix4d();
		trans.setIdentity();
		trans.setTranslation(new Vector3d(xtrans));
		transformation.mul(transformation, trans);
		// System.out.println("setting xtrans");
		// System.out.println(transformation);

		// combine with origin -> y translation
		ytrans.negate();
		Matrix4d transInverse = new Matrix4d();
		transInverse.setIdentity();
		transInverse.setTranslation(new Vector3d(ytrans));
		transformation.mul(transInverse, transformation);
		// System.out.println("setting ytrans");
		// System.out.println(transformation);
	}

	/**
	 * Calculates the inner product between two coordinate sets x and y
	 * (optionally weighted, if weights set through
	 * {@link #set(Point3d[], Point3d[], double[])}). It also calculates an
	 * upper bound of the most positive root of the key matrix.
	 * http://theobald.brandeis.edu/qcp/qcprot.c
	 * 
	 * @param coords1
	 * @param coords2
	 * @return
	 */
	private void innerProduct(Point3d[] coords1, Point3d[] coords2) {
		double x1, x2, y1, y2, z1, z2;
		double g1 = 0.0, g2 = 0.0;

		Sxx = 0;
		Sxy = 0;
		Sxz = 0;
		Syx = 0;
		Syy = 0;
		Syz = 0;
		Szx = 0;
		Szy = 0;
		Szz = 0;

		if (weight != null) {
			wsum = 0;
			for (int i = 0; i < coords1.length; i++) {

				wsum += weight[i];

				x1 = weight[i] * coords1[i].x;
				y1 = weight[i] * coords1[i].y;
				z1 = weight[i] * coords1[i].z;

				g1 += x1 * coords1[i].x + y1 * coords1[i].y + z1 * coords1[i].z;

				x2 = coords2[i].x;
				y2 = coords2[i].y;
				z2 = coords2[i].z;

				g2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

				Sxx += (x1 * x2);
				Sxy += (x1 * y2);
				Sxz += (x1 * z2);

				Syx += (y1 * x2);
				Syy += (y1 * y2);
				Syz += (y1 * z2);

				Szx += (z1 * x2);
				Szy += (z1 * y2);
				Szz += (z1 * z2);
			}
		} else {
			for (int i = 0; i < coords1.length; i++) {
				g1 += coords1[i].x * coords1[i].x + coords1[i].y * coords1[i].y + coords1[i].z * coords1[i].z;
				g2 += coords2[i].x * coords2[i].x + coords2[i].y * coords2[i].y + coords2[i].z * coords2[i].z;

				Sxx += coords1[i].x * coords2[i].x;
				Sxy += coords1[i].x * coords2[i].y;
				Sxz += coords1[i].x * coords2[i].z;

				Syx += coords1[i].y * coords2[i].x;
				Syy += coords1[i].y * coords2[i].y;
				Syz += coords1[i].y * coords2[i].z;

				Szx += coords1[i].z * coords2[i].x;
				Szy += coords1[i].z * coords2[i].y;
				Szz += coords1[i].z * coords2[i].z;
			}
			wsum = coords1.length;
		}

		e0 = (g1 + g2) * 0.5;
	}

	private int calcRmsd(double len) {
		double Sxx2 = Sxx * Sxx;
		double Syy2 = Syy * Syy;
		double Szz2 = Szz * Szz;

		double Sxy2 = Sxy * Sxy;
		double Syz2 = Syz * Syz;
		double Sxz2 = Sxz * Sxz;

		double Syx2 = Syx * Syx;
		double Szy2 = Szy * Szy;
		double Szx2 = Szx * Szx;

		double SyzSzymSyySzz2 = 2.0 * (Syz * Szy - Syy * Szz);
		double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

		double c2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
		double c1 = 8.0 * (Sxx * Syz * Szy + Syy * Szx * Sxz + Szz * Sxy * Syx - Sxx * Syy * Szz - Syz * Szx * Sxy
				- Szy * Syx * Sxz);

		SxzpSzx = Sxz + Szx;
		SyzpSzy = Syz + Szy;
		SxypSyx = Sxy + Syx;
		SyzmSzy = Syz - Szy;
		SxzmSzx = Sxz - Szx;
		SxymSyx = Sxy - Syx;
		SxxpSyy = Sxx + Syy;
		SxxmSyy = Sxx - Syy;

		double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

		double c0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
				+ (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
				+ (-(SxzpSzx) * (SyzmSzy) + (SxymSyx) * (SxxmSyy - Szz))
						* (-(SxzmSzx) * (SyzpSzy) + (SxymSyx) * (SxxmSyy + Szz))
				+ (-(SxzpSzx) * (SyzpSzy) - (SxypSyx) * (SxxpSyy - Szz))
						* (-(SxzmSzx) * (SyzmSzy) - (SxypSyx) * (SxxpSyy + Szz))
				+ (+(SxypSyx) * (SyzpSzy) + (SxzpSzx) * (SxxmSyy + Szz))
						* (-(SxymSyx) * (SyzmSzy) + (SxzpSzx) * (SxxpSyy + Szz))
				+ (+(SxypSyx) * (SyzmSzy) + (SxzmSzx) * (SxxmSyy - Szz))
						* (-(SxymSyx) * (SyzpSzy) + (SxzmSzx) * (SxxpSyy - Szz));

		mxEigenV = e0;

		int i;
		for (i = 1; i < 51; ++i) {
			double oldg = mxEigenV;
			double x2 = mxEigenV * mxEigenV;
			double b = (x2 + c2) * mxEigenV;
			double a = b + c1;
			double delta = ((a * mxEigenV + c0) / (2.0 * x2 * mxEigenV + b + a));
			mxEigenV -= delta;

			if (Math.abs(mxEigenV - oldg) < Math.abs(eval_prec * mxEigenV))
				break;
		}

		if (i == 50) {
			logger.warn(String.format("More than %d iterations needed!", i));
		} else {
			logger.info(String.format("%d iterations needed!", i));
		}

		/*
		 * the fabs() is to guard against extremely small, but *negative*
		 * numbers due to floating point error
		 */
		rmsd = Math.sqrt(Math.abs(2.0 * (e0 - mxEigenV) / len));

		return 1;
	}

	private int calcRotationMatrix() {
		double a11 = SxxpSyy + Szz - mxEigenV;
		double a12 = SyzmSzy;
		double a13 = -SxzmSzx;
		double a14 = SxymSyx;
		double a21 = SyzmSzy;
		double a22 = SxxmSyy - Szz - mxEigenV;
		double a23 = SxypSyx;
		double a24 = SxzpSzx;
		double a31 = a13;
		double a32 = a23;
		double a33 = Syy - Sxx - Szz - mxEigenV;
		double a34 = SyzpSzy;
		double a41 = a14;
		double a42 = a24;
		double a43 = a34;
		double a44 = Szz - SxxpSyy - mxEigenV;
		double a3344_4334 = a33 * a44 - a43 * a34;
		double a3244_4234 = a32 * a44 - a42 * a34;
		double a3243_4233 = a32 * a43 - a42 * a33;
		double a3143_4133 = a31 * a43 - a41 * a33;
		double a3144_4134 = a31 * a44 - a41 * a34;
		double a3142_4132 = a31 * a42 - a41 * a32;
		double q1 = a22 * a3344_4334 - a23 * a3244_4234 + a24 * a3243_4233;
		double q2 = -a21 * a3344_4334 + a23 * a3144_4134 - a24 * a3143_4133;
		double q3 = a21 * a3244_4234 - a22 * a3144_4134 + a24 * a3142_4132;
		double q4 = -a21 * a3243_4233 + a22 * a3143_4133 - a23 * a3142_4132;

		double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

		/*
		 * The following code tries to calculate another column in the adjoint
		 * matrix when the norm of the current column is too small. Usually this
		 * commented block will never be activated. To be absolutely safe this
		 * should be uncommented, but it is most likely unnecessary.
		 */
		if (qsqr < evec_prec) {
			q1 = a12 * a3344_4334 - a13 * a3244_4234 + a14 * a3243_4233;
			q2 = -a11 * a3344_4334 + a13 * a3144_4134 - a14 * a3143_4133;
			q3 = a11 * a3244_4234 - a12 * a3144_4134 + a14 * a3142_4132;
			q4 = -a11 * a3243_4233 + a12 * a3143_4133 - a13 * a3142_4132;
			qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

			if (qsqr < evec_prec) {
				double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
				double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
				double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

				q1 = a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
				q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
				q3 = a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
				q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
				qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

				if (qsqr < evec_prec) {
					q1 = a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
					q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
					q3 = a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
					q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
					qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

					if (qsqr < evec_prec) {
						/*
						 * if qsqr is still too small, return the identity
						 * matrix.
						 */
						rotmat.setIdentity();

						return 0;
					}
				}
			}
		}

		double normq = Math.sqrt(qsqr);
		q1 /= normq;
		q2 /= normq;
		q3 /= normq;
		q4 /= normq;

		logger.debug("q: " + q1 + " " + q2 + " " + q3 + " " + q4);

		double a2 = q1 * q1;
		double x2 = q2 * q2;
		double y2 = q3 * q3;
		double z2 = q4 * q4;

		double xy = q2 * q3;
		double az = q1 * q4;
		double zx = q4 * q2;
		double ay = q1 * q3;
		double yz = q3 * q4;
		double ax = q1 * q2;

		rotmat.m00 = a2 + x2 - y2 - z2;
		rotmat.m01 = 2 * (xy + az);
		rotmat.m02 = 2 * (zx - ay);

		rotmat.m10 = 2 * (xy - az);
		rotmat.m11 = a2 - x2 + y2 - z2;
		rotmat.m12 = 2 * (yz + ax);

		rotmat.m20 = 2 * (zx + ay);
		rotmat.m21 = 2 * (yz - ax);
		rotmat.m22 = a2 - x2 - y2 + z2;

		return 1;
	}

	@Override
	public double getRmsd(Point3d[] fixed, Point3d[] moved) {
		set(moved, fixed);
		return getRmsd();
	}

	@Override
	public Matrix4d superpose(Point3d[] fixed, Point3d[] moved) {
		set(moved, fixed);
		getRotationMatrix();
		if (!centered) {
			calcTransformation();
		} else {
			transformation.set(rotmat);
		}
		return transformation;
	}

	/**
	 * @param fixed
	 * @param moved
	 * @param weight
	 *            array of weigths for each equivalent point position
	 * @return weighted RMSD.
	 */
	public double getWeightedRmsd(Point3d[] fixed, Point3d[] moved, double[] weight) {
		set(moved, fixed, weight);
		return getRmsd();
	}

	/**
	 * The QCP method can be used as a two-step calculation: first compute the
	 * RMSD (fast) and then compute the superposition.
	 * 
	 * This method assumes that the RMSD of two arrays of points has been
	 * already calculated using {@link #getRmsd(Point3d[], Point3d[])} method
	 * and calculates the transformation of the same two point arrays.
	 * 
	 * @param fixed
	 * @param moved
	 * @return transformation matrix as a Matrix4d to superpose moved onto fixed
	 *         point arrays
	 */
	public Matrix4d superposeAfterRmsd() {

		if (!rmsdCalculated) {
			throw new IllegalStateException("The RMSD was not yet calculated. Use the superpose() method instead.");
		}

		getRotationMatrix();
		if (!centered) {
			calcTransformation();
		} else {
			transformation.set(rotmat);
		}
		return transformation;
	}

}
