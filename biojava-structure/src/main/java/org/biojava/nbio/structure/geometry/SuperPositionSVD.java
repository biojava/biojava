/*
 *                  BioJava development code
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
 * Created on Dec 4, 2005
 *
 */
package org.biojava.nbio.structure.geometry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.jama.SingularValueDecomposition;

/**
 * A class that calculates the superposition between two sets of points using an
 * SVD Matrix Decomposition. It was introduced by Wolfgang Kabsch, hence the
 * alternative name Kabsh algorithm. Inspired by the biopython SVDSuperimposer
 * class.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 * @since 1.5
 * @version %I% %G%
 * 
 */
public class SuperPositionSVD extends SuperPositionAbstract {

	/**
	 * Constructor for the SVD superposition algorithm.
	 * 
	 * @param centered
	 *            true if the point arrays are centered at the origin (faster),
	 *            false otherwise
	 */
	public SuperPositionSVD(boolean centered) {
		super(centered);
	}

	@Override
	public Matrix4d superpose(Point3d[] fixed, Point3d[] moved) {

		checkInput(fixed, moved);

		Point3d cena = CalcPoint.centroid(fixed);
		Point3d cenb = CalcPoint.centroid(moved);

		double[][] centAcoords = new double[][] { { cena.x, cena.y, cena.z } };
		Matrix centroidA = new Matrix(centAcoords);

		double[][] centBcoords = new double[][] { { cenb.x, cenb.y, cenb.z } };
		Matrix centroidB = new Matrix(centBcoords);

		// center at centroid

		cena.negate();
		cenb.negate();

		Point3d[] ats1 = CalcPoint.clonePoint3dArray(fixed);
		CalcPoint.translate(new Vector3d(cena), ats1);

		Point3d[] ats2 = CalcPoint.clonePoint3dArray(moved);
		CalcPoint.translate(new Vector3d(cenb), ats2);

		double[][] coordSet1 = new double[ats1.length][3];
		double[][] coordSet2 = new double[ats2.length][3];

		// copy the atoms into the internal coords;
		for (int i = 0; i < ats1.length; i++) {
			coordSet1[i] = new double[3];
			ats1[i].get(coordSet1[i]);
			coordSet2[i] = new double[3];
			ats2[i].get(coordSet2[i]);
		}

		// now this is the bridge to the Jama package:
		Matrix a = new Matrix(coordSet1);
		Matrix b = new Matrix(coordSet2);

		// # correlation matrix

		Matrix b_trans = b.transpose();
		Matrix corr = b_trans.times(a);

		SingularValueDecomposition svd = corr.svd();

		Matrix u = svd.getU();
		// v is alreaady transposed ! difference to numermic python ...
		Matrix vt = svd.getV();

		Matrix vt_orig = (Matrix) vt.clone();
		Matrix u_transp = u.transpose();

		Matrix rot_nottrans = vt.times(u_transp);
		Matrix rot = rot_nottrans.transpose();

		// check if we have found a reflection

		double det = rot.det();

		if (det < 0) {
			vt = vt_orig.transpose();
			vt.set(2, 0, (0 - vt.get(2, 0)));
			vt.set(2, 1, (0 - vt.get(2, 1)));
			vt.set(2, 2, (0 - vt.get(2, 2)));

			Matrix nv_transp = vt.transpose();
			rot_nottrans = nv_transp.times(u_transp);
			rot = rot_nottrans.transpose();

		}

		Matrix cb_tmp = centroidB.times(rot);
		Matrix tran = centroidA.minus(cb_tmp);
		
		return Matrices.getTransformation(rot, tran);

	}

	@Override
	public double getRmsd(Point3d[] x, Point3d[] y) {
		Point3d[] yclone = CalcPoint.clonePoint3dArray(y);
		superposeAndTransform(x, yclone);
		return CalcPoint.rmsd(x, yclone);
	}

}
