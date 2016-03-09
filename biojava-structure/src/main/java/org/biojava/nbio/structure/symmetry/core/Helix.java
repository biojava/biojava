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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Peter
 */
public class Helix {
	private QuatSymmetryScores scores = new QuatSymmetryScores();
	private List<Integer> permutation;
	private List<List<Integer>> repeatUnits;
	private Matrix4d transformation;
	private double rise;
	private int nStart;
	private int fold;
	private int contacts;

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

	public List<List<Integer>> getRepeatUnits() {
		return repeatUnits;
	}

	public void setRepeatUnits(List<List<Integer>> repeatUnits) {
		this.repeatUnits = repeatUnits;
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

	public double getRise() {
		return rise;
	}

	public void setRise(double rise) {
		this.rise = rise;
	}

	/**
	 * Returns the pitch angle of the helix
	 * @param transformation helix transformation
	 * @return
	 */
	public double getAngle() {
		return getAxisAngle().angle;
	}

	/**
	 * Returns the AxisAngle of the helix transformation
	 * @param transformation helix transformation
	 * @return
	 */
	public AxisAngle4d getAxisAngle() {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(this.transformation);
		return axis;
	}

	public int getnStart() {
		return nStart;
	}

	public void setnStart(int nStart) {
		this.nStart = nStart;
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

	public int getContacts() {
		return contacts;
	}

	public void setContacts(int contacts) {
		this.contacts = contacts;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Permutation   : " + getPermutation() + "\n");
		sb.append("Repeat units  : " + getRepeatUnits() + "\n");
		sb.append("Rise          : " + getRise() + "\n");
		sb.append("Angle         : " + Math.toDegrees(getAngle()) +"\n");
		sb.append("Fold          : " + getFold() + "\n");
		return sb.toString();
	}

	public List<List<Integer>> getLayerLines() {
		List<List<Integer>> layerLines = new ArrayList<List<Integer>>();

		createLineSegments(permutation, layerLines);

//		System.out.println("Line segments: " + layerLines.size());
//		for (List<Integer> lineSegment: layerLines) {
//			System.out.println(lineSegment);
//		}

		int count = layerLines.size();

		// iteratively join line segments
		do {
			count = layerLines.size();
			joinLineSegments(layerLines);
			// after joining line segments, get rid of the empty line segments left behind
			trimEmptyLineSegments(layerLines);

//			System.out.println("Line segments: " + count);
//			for (List<Integer> lineSegment: layerLines) {
//				System.out.println(lineSegment);
//			}
		} while (layerLines.size() < count);

		return layerLines;
	}

	private static void createLineSegments(List<Integer> permutation,
			List<List<Integer>> layerLines) {
		for (int i = 0; i < permutation.size(); i++) {
			if (permutation.get(i) != -1 ) {
				List<Integer> lineSegment = new ArrayList<Integer>();
				lineSegment.add(i);
				lineSegment.add(permutation.get(i));
				layerLines.add(lineSegment);
			}
		}
	}

	private static void joinLineSegments(List<List<Integer>> layerLines) {
		for (int i = 0; i < layerLines.size()-1; i++) {
			List<Integer> lineSegmentI = layerLines.get(i);
			if (! lineSegmentI.isEmpty()) {
				for (int j = i + 1; j < layerLines.size(); j++) {
					List<Integer> lineSegmentJ = layerLines.get(j);
					if (! lineSegmentJ.isEmpty()) {
						if (lineSegmentI.get(lineSegmentI.size()-1).equals(lineSegmentJ.get(0))) {
//							System.out.println("join right: " + lineSegmentI + " - " + lineSegmentJ);
							lineSegmentI.addAll(lineSegmentJ.subList(1,  lineSegmentJ.size()));
//							System.out.println("joned segment: " + lineSegmentI);
							lineSegmentJ.clear();
						} else if ((lineSegmentI.get(0).equals(lineSegmentJ.get(lineSegmentJ.size()-1)))) {
							lineSegmentI.addAll(0, lineSegmentJ.subList(0,  lineSegmentJ.size()-1));
//							System.out.println("join left: " + lineSegmentJ + " - " + lineSegmentI);
//							System.out.println("joned segment: " + lineSegmentI);
							lineSegmentJ.clear();
						}
					}
				}
			}
		}
	}

	private static void trimEmptyLineSegments(List<List<Integer>> layerLines) {
		for (Iterator<List<Integer>> iter = layerLines.iterator(); iter.hasNext();) {
			if (iter.next().isEmpty()) {
				iter.remove();
			}
		}
	}
}
