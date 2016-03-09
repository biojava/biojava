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
import java.util.List;

/**
 *
 * @author Peter
 */
public class Rotation {
	private double subunitRmsd = Double.MAX_VALUE;
	private double traceRmsd = Double.MAX_VALUE;
	private double traceTmScoreMin = Double.MAX_VALUE;
	private QuatSymmetryScores scores = new QuatSymmetryScores();
	private List<Integer> permutation;
	private Matrix4d transformation;
	private AxisAngle4d axisAngle;
	private int direction;
	private int fold;

	/**
	 * @return the subunitRmsd
	 */
	public double getSubunitRmsd() {
		return subunitRmsd;
	}

	/**
	 * @param subunitRmsd the subunitRmsd to set
	 */
	public void setSubunitRmsd(double subunitRmsd) {
		this.subunitRmsd = subunitRmsd;
	}

	/**
	 * @return the traceRmsd
	 */
	public double getTraceRmsd() {
		return traceRmsd;
	}

	/**
	 * @param traceRmsd the traceRmsd to set
	 */
	public void setTraceRmsd(double traceRmsd) {
		this.traceRmsd = traceRmsd;
	}

	/**
	 * @param traceTmScoreMin the traceTmScore to set
	 */
	public void setTraceTmScoreMin(double traceTmScoreMin) {
		this.traceTmScoreMin = traceTmScoreMin;
	}

	/**
	 * @return the traceTmScoreMin
	 */
	public double getTraceTmScoreMin() {
		return traceTmScoreMin;
	}


	/**
	 * @return the permutation
	 */
	public List<Integer> getPermutation() {
		return permutation;
	}

	/**
	 * @param permutation the permutation to set
	 */
	public void setPermutation(List<Integer> permutation) {
		this.permutation = permutation;
	}

	/**
	 * @return the transformation
	 */
	public Matrix4d getTransformation() {
		return transformation;
	}

	/**
	 * @param transformation the transformation to set
	 */
	public void setTransformation(Matrix4d transformation) {
		this.transformation = transformation;
	}

	/**
	 * @return the fold
	 */
	public int getFold() {
		return fold;
	}

	/**
	 * @param fold the fold to set
	 */
	public void setFold(int fold) {
		this.fold = fold;
	}

	/**
	 * @return the scores
	 */
	public QuatSymmetryScores getScores() {
		return scores;
	}

	/**
	 * @param scores the scores to set
	 */
	public void setScores(QuatSymmetryScores scores) {
		this.scores = scores;
	}

	/**
	 * @return the direction
	 */
	public int getDirection() {
		return direction;
	}

	/**
	 * @param direction the direction to set
	 */
	public void setDirection(int axis) {
		this.direction = axis;
	}

	/**
	 * @return the axisAngle
	 */
	public AxisAngle4d getAxisAngle() {
		return axisAngle;
	}

	/**
	 * @param axisAngle the axisAngle to set
	 */
	public void setAxisAngle(AxisAngle4d axisAngle) {
		this.axisAngle = axisAngle;
	}

	/**
	 * Returns the number of starts if this rotation represents a helical rotation
	 */
	public int getNStart() {
		int nStart = 0;
		for (int i: permutation) {
			if (i == -1) {
				nStart++;
			}
		}
		return nStart;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("fold       : ");
		sb.append(fold);
		sb.append("/n");
		sb.append("orientation: ");
		sb.append("direction  : ");
		sb.append("/n");
		sb.append("axisAngle  : ");
		sb.append(axisAngle);
		sb.append("/n");
		sb.append("permutation: ");
		sb.append(permutation);
		sb.append(scores);
		return sb.toString();
	}
}
