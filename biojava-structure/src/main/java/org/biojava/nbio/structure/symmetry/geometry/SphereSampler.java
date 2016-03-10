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

import javax.vecmath.*;
import java.util.ArrayList;
import java.util.List;


// Generate the permutations and sign changes for a Triple.
class Permute {

	private List<Point3i> triples = new ArrayList<Point3i>();

	Permute(Point3i t) {
		Point3i tmp = new Point3i();
		tmp.x = t.x;
		tmp.y = t.y;
		tmp.z = t.z;
		triples.add(tmp);
		int n = 1;

		if (t.x != 0) {
			for (int i = 0; i < n; ++i) {
				Tuple3i m = triples.get(i);

				triples.add(new Point3i(-m.x, m.y, m.z));
			}
			n *= 2;
		}

		if (t.y != 0) {
			for (int i = 0; i < n; ++i) {
				Point3i m = triples.get(i);
				triples.add(new Point3i(m.x, -m.y, m.z));
			}
			n *= 2;
		}

		if (t.z != 0) {
			for (int i = 0; i < n; ++i) {
				Point3i m = triples.get(i);
				triples.add(new Point3i(m.x, m.y, -m.z));
			}
			n *= 2;
		}
		if (t.x == t.y && t.y == t.z) {
			return;
		}

		for (int i = 0; i < n; ++i) {
			Point3i m = triples.get(i);
			triples.add(new Point3i(m.y, m.z, m.x));
			triples.add(new Point3i(m.z, m.x, m.y));
		}
		n *= 3;

		if (t.x == t.y || t.y == t.z) {
			return;
		}

		for (int i = 0; i < n; ++i) {
			Point3i m = triples.get(i);
			triples.add(new Point3i(m.y, m.x, m.z));
		}
		n *= 2;
	}

	public int size() {
		return triples.size();
	}

	public Point3i get(int i) {
		return triples.get(i);
	}
};

/**
 *
 * @author Peter
 */
public final class SphereSampler {

	private static final List<Quat4d> orientations ;



	// The rotational symmetries of the cube. (Not normalized, since
	// PackSet.Add does this.)
	private static final double cubeSyms[][] = {
		{ 1, 0, 0, 0 },
		// 180 deg rotations about 3 axes
		{ 0, 1, 0, 0 },
		{ 0, 0, 1, 0 },
		{ 0, 0, 0, 1 },
		// +/- 120 degree rotations about 4 leading diagonals
		{ 1, 1, 1, 1 }, { 1, 1, 1, -1 }, { 1, 1, -1, 1 }, { 1, 1, -1, -1 },
		{ 1, -1, 1, 1 }, { 1, -1, 1, -1 },
		{ 1, -1, -1, 1 },
		{ 1, -1, -1, -1 },
		// +/- 90 degree rotations about 3 axes
		{ 1, 1, 0, 0 }, { 1, -1, 0, 0 }, { 1, 0, 1, 0 }, { 1, 0, -1, 0 },
		{ 1, 0, 0, 1 }, { 1, 0, 0, -1 },
		// 180 degree rotations about 6 face diagonals
		{ 0, 1, 1, 0 }, { 0, 1, -1, 0 }, { 0, 1, 0, 1 }, { 0, 1, 0, -1 },
		{ 0, 0, 1, 1 }, { 0, 0, 1, -1 }, };



	private static final double delta = 0.15846;
	private static final double sigma = 0.00;
	private static final int ntot = 7416;
	private static final int ncell = 309;
	private static final int nent = 18;

	private static final int[] k = { 0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5};
	private static final int[] l = { 0, 1, 0, 2, 2, 1, 3, 3, 0, 2, 2, 4, 4, 4, 1, 3, 3, 5};
	private static final int[] m = { 0, 1, 0, 0, 2, 1, 1, 3, 0, 0, 2, 0, 2, 4, 1, 1, 3, 1};

	private static final int[] mult = { 1, 8, 6, 12, 8, 24, 24, 8, 6, 24, 24, 12, 24, 8, 24, 48, 24, 24};

	static
	{

		List<Quat4d> myorientations = new ArrayList<Quat4d>();

		for (int i = 0; i < IcosahedralSampler.getSphereCount(); i++) {
			myorientations.add(IcosahedralSampler.getQuat4d(i));
		}
		List<Quat4d> grid = new ArrayList<Quat4d>();
		int ncell1 = 0;
		for (int n = 0; n < nent; ++n) {
			Permute p = new Permute(new Point3i(k[n], l[n], m[n]));
			assert (mult[n] == p.size());
			for (int i = 0; i < mult[n]; ++i) {
				Point3i t = p.get(i);
				grid.add(new Quat4d(1.0, pind(0.5 * t.x, delta, sigma), pind(
						0.5 * t.y, delta, sigma), pind(0.5 * t.z, delta, sigma)));
			}
			ncell1 += mult[n];
		}
		assert (ncell1 == ncell);
		int nc = grid.size();
		assert (nc == ncell);
		for (int n = 1; n < 24; ++n) {
			Quat4d q = new Quat4d(cubeSyms[n][0], cubeSyms[n][1],
					cubeSyms[n][2], cubeSyms[n][3]);
			for (int i = 0; i < nc; ++i) {
				Quat4d qs = new Quat4d();
				qs.mul(q, grid.get(i));
				grid.add(qs);
				//	s.add(times(q, s.getOrientation(i)), s.getWeight(i)); // this data set has no weights
			}
		}
		assert (grid.size() == ntot);
		myorientations.addAll(grid);

		orientations = myorientations;

	}

	// this class cannot be instantiated
	private SphereSampler() {
	};

	public static int getSphereCount() {

		return orientations.size();
	}

	public static Quat4d getQuat4d(int index) {

		return orientations.get(index);
	}

	public static void getAxisAngle(int index, AxisAngle4f axisAngle) {

		axisAngle.set(orientations.get(index));
	}

	public static void getAxisAngle(int index, AxisAngle4d axisAngle) {

		axisAngle.set(orientations.get(index));
	}

	// Convert from index to position. The sinh scaling tries to compensate
	// for the bunching up that occurs when [1 x y z] is projected onto the
	// unit sphere.
	private static double pind(double ind, double delta, double sigma) {
		return (sigma == 0) ? ind * delta : Math.sinh(sigma * ind * delta)
				/ sigma;
	}

}
