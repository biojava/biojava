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

package org.biojava.nbio.structure.symmetry.core;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * @see http://en.wikipedia.org/wiki/Rotation_group_SO(3)
 * @author Peter
 */
public class RotationGroup implements Iterable<Rotation> {
	private List<Rotation> rotations = new ArrayList<Rotation>();
	private int principalAxisIndex = 0;
	private int higherOrderRotationAxis = 0;
	private int twoFoldsPerpendicular = 0;
	private int highestOrder = 0;
	private String pointGroup = "C1";
	private double symmetryDeviation = 0;
	private boolean complete = true;
	private boolean modified = true;


	public int getOrder() {
		return rotations.size();
	}

	public Rotation getRotation(int index) {
		return rotations.get(index);
	}

	public void addRotation(Rotation rotation) {
		rotations.add(rotation);
		modified = true;
	}

	public void setC1(int n) {
		Rotation r = new Rotation();
		List<Integer> permutation = new ArrayList<Integer>(n);
		for (int i = 0; i < n; i++) {
			permutation.add(i);
		}
		r.setPermutation(permutation);
		Matrix4d m = new Matrix4d();
		m.setIdentity();
		r.setTransformation(m);
		r.setAxisAngle(new AxisAngle4d());
		r.setFold(1);
		r.setScores(new QuatSymmetryScores());
		rotations.add(r);
		pointGroup = "C1";
	}

	public void removeRotation(int index) {
		rotations.remove(index);
		modified = true;
	}

	public void complete() {
		if (modified) {
			if (rotations.size() > 0) {
				findHighestOrderAxis();
				setEAxis();
				calcAxesDirections();
				findHigherOrderAxes();
				findTwoFoldsPerpendicular();
				calcPointGroup();
				sortByFoldDecending();
			}
			modified = false;
		}
	}

	public String getPointGroup() {
		if (modified) {
			if (rotations.size() == 0) {
				return "C1";
			}
			complete();
		}
		return pointGroup;
	}

	/**
	 * Returns QuatSymmetryScores averaged over all rotations
	 * (except the first rotation, which is the unit operation E)
	 * @return mean scores average over rotations
	 */
	public QuatSymmetryScores getScores() {
		QuatSymmetryScores scores = new QuatSymmetryScores();

		int n = rotations.size()-1;

		if (n > 0) {
			double[] values = new double[n];

			// minRmsd
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getMinRmsd();
			}
			scores.setMinRmsd(minScores(values));

			// maxRmsd
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getMaxRmsd();
			}
			scores.setMaxRmsd(maxScores(values));

			// Rmsd
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getRmsd();
			}
			scores.setRmsd(averageScores(values));

			// minTm
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getMinTm();
			}
			scores.setMinTm(minScores(values));

			// maxTm
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getMaxTm();
			}
			scores.setMaxTm(maxScores(values));

			// Tm
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getTm();
			}
			scores.setTm(averageScores(values));

			// Rmsd subunit centers
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getRmsdCenters();
			}
			scores.setRmsdCenters(averageScores(values));
			// TmIntra
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getTmIntra();
			}
			scores.setTmIntra(averageScores(values));

			// RmsdIntra
			for (int i = 1; i < rotations.size(); i++) {
				values[i-1] = rotations.get(i).getScores().getRmsdIntra();
			}
			scores.setRmsdIntra(averageScores(values));

			// SymDeviation
			scores.setSymDeviation(symmetryDeviation);
		}
		return scores;
	}

	/**
	 * @param symmetryDeviation the symmetryDeviation to set
	 */
	public void setSymmetryDeviation(double symmetryDeviation) {
		this.symmetryDeviation = symmetryDeviation;
	}

	public boolean isComplete() {
		return complete;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Rotations: " + rotations.size() + "\n");
		for (Rotation s: rotations) {
			sb.append(s.toString()).append("\n");
		}
		return sb.toString();
	}

	private double averageScores(double[] scores) {
		double sum = 0;
		for (double s: scores) {
			sum += s;
		}
		return sum/scores.length;
	}

	private double minScores(double[] scores) {
		double score = Double.MAX_VALUE;
		for (double s: scores) {
			score = Math.min(score, s);
		}
		return score;
	}

	private double maxScores(double[] scores) {
		double score = Double.MIN_VALUE;
		for (double s: scores) {
			score = Math.max(score, s);
		}
		return score;
	}

	private void findHighestOrderAxis() {
		highestOrder = 1;
		principalAxisIndex = 0;
		double rmsd  = Double.MAX_VALUE;

		for (int i = 0; i < rotations.size(); i++) {
			Rotation s = rotations.get(i);
			if (s.getFold() > highestOrder) {
				highestOrder = s.getFold();
				principalAxisIndex = i;
				rmsd = s.getTraceRmsd();
			} else if (s.getFold() >= highestOrder && s.getTraceRmsd() < rmsd) {
				highestOrder = s.getFold();
				principalAxisIndex = i;
				rmsd = s.getTraceRmsd();
			}
		}
	}

	/**
	 * Add E operation to the highest order rotation axis. By definition
	 * E belongs to the highest order axis.
	 */
	private void setEAxis() {
		Rotation e = rotations.get(0);
		Rotation h = rotations.get(principalAxisIndex);
		e.setAxisAngle(new AxisAngle4d(h.getAxisAngle()));
		e.getAxisAngle().angle = 0.0;
		e.setFold(h.getFold());
	}

	private void findHigherOrderAxes() {
		higherOrderRotationAxis = 0;
		for (Rotation s: rotations) {
			if (s.getFold() > 2 && s.getDirection() == 1) {
				higherOrderRotationAxis++;
			}
		}
	}

	private void calcAxesDirections() {
		if (highestOrder == 1) {
			for (Rotation s: rotations) {
				s.setDirection(0);
			}
			return;
		}

		AxisAngle4d pa = rotations.get(principalAxisIndex).getAxisAngle();
		Vector3d pv = new Vector3d(pa.x, pa.y, pa.z);

		for (Rotation s: rotations) {
			AxisAngle4d axis = s.getAxisAngle();
			Vector3d av = new Vector3d(axis.x, axis.y, axis.z);
			if (Math.abs(pv.dot(av)) > 0.9f) {
				// co-linear with principal axis
				s.setDirection(0);
			} else {
				// not co-linear or perpendicular to principal axis
				s.setDirection(1);
			}
		}
		rotations.get(0).setDirection(0); // set the E axis to the principal axis (by definition)
	}

	private void findTwoFoldsPerpendicular() {
		twoFoldsPerpendicular = 0;
		for (Rotation s: rotations) {
			if (s.getFold() == 2 && s.getDirection() == 1) {
				twoFoldsPerpendicular++;
			}
		}
	}


	public int getHigherOrderRotationAxis(){
		return higherOrderRotationAxis;
	}

	public int getTwoFoldsPerpendicular(){
		return twoFoldsPerpendicular;
	}

	public int getPrincipalAxisIndex(){
		return principalAxisIndex;
	}

	private void calcPointGroup() {
		complete = false;
		if (higherOrderRotationAxis > 1) {
			// cubic groups
			if (highestOrder == 5) {
				// rotational icosahedral symmetry or chiral icosahedral symmetry
				pointGroup = "I";
				complete = rotations.size() == 60;
			} else if (highestOrder == 4) {
				// rotational octahedral symmetry or chiral octahedral symmetry
				pointGroup = "O";
				complete = rotations.size() == 24;
			} else if (highestOrder == 3) {
				// rotational tetrahedral symmetry or chiral tetrahedral symmetry
				pointGroup = "T";
				complete = rotations.size() == 12;
			}
		} else {
			// Cn and Dn groups
			// if E is not counted, subtract 1
			if (Math.abs(twoFoldsPerpendicular - highestOrder) <= 1 && highestOrder > 1) {
				pointGroup = "D" + highestOrder;
				complete = rotations.size() == 2 * highestOrder;
			} else {
				pointGroup = "C" + highestOrder;
				complete = rotations.size() == highestOrder;
			}
		}

		if (!complete) {
			// the rotation group is incomplete, remove partial results. This happens
			// when a structure is symmetric, some subunits are below the rmsd threshold,
			// and some are just above the rmsd threshold
			int n = 0;
			if (rotations.size() > 0) {
				n = rotations.get(0).getPermutation().size();
				rotations.clear();
			}
			setC1(n);
		}
	}

	public void sortByFoldDecending() {
		Collections.sort(rotations, new Comparator<Rotation>() {
			@Override
			public int compare(Rotation o1, Rotation o2) {
				int delta = o1.getDirection() - o2.getDirection();
				if (delta != 0) {
					return delta;
				}
				delta = Math.round(Math.signum(o2.getFold() - o1.getFold()));
				if (delta != 0) {
					return delta;
				}

				delta = (int)(Math.signum(o1.getAxisAngle().angle - o2.getAxisAngle().angle));
				return delta;
			}
		});
	}

	@Override
	public Iterator<Rotation> iterator() {
		return rotations.iterator();
	}
}
