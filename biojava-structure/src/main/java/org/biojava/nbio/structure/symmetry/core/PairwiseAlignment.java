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

import java.util.Arrays;

public class PairwiseAlignment {
	private SequenceAlignmentCluster cluster1 = null;
	private SequenceAlignmentCluster cluster2 = null;
	private double alignmentLengthFraction = 0;
	private double sequenceIdentity = 0;
	private double rmsd = 0;
	private int[][][] alignment = null;

	public PairwiseAlignment(SequenceAlignmentCluster cluster1, SequenceAlignmentCluster cluster2) {
		this.cluster1 = cluster1;
		this.cluster2 = cluster2;
	}

	public SequenceAlignmentCluster getCluster1() {
		return cluster1;
	}

	public SequenceAlignmentCluster getCluster2() {
		return cluster2;
	}

	public double getAlignmentLengthFraction() {
		return alignmentLengthFraction;
	}

	public double getSequenceIdentity() {
		return sequenceIdentity;
	}

	public double getRmsd() {
		return rmsd;
	}

	public int[][][] getAlignment() {
		return alignment;
	}

	public void setAlignmentLengthFraction(double alignmentLengthFraction) {
		this.alignmentLengthFraction = alignmentLengthFraction;
	}

	public void setSequenceIdentity(double sequenceIdentity) {
		this.sequenceIdentity = sequenceIdentity;
	}

	public void setRmsd(double rmsd) {
		this.rmsd = rmsd;
	}

	public void setAlignment(int[][][] alignment) {
		this.alignment = alignment;
	}

	@Override
	public String toString() {
		return new StringBuffer()
				.append("cluster1:")
				.append("\n")
				.append(cluster1)
				.append("\n")
				.append("cluster2:")
				.append("\n")
				.append(cluster2)
				.append("\n")
				.append("sequence identity: " + sequenceIdentity)
				.append("\n")
				.append("alignment fraction: " + alignmentLengthFraction)
				.append("\n")
				.append("rmsd: " + rmsd)
				.append("\n")
				.append("aligment1: " + Arrays.toString(alignment[0][0]))
				.append("\n")
				.append("aligment2: " + Arrays.toString(alignment[0][1]))
				.append("\n").toString();
	}
}
