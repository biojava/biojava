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
package org.biojava.nbio.structure.cluster;

import java.io.Serializable;
import java.util.Arrays;

/**
 * The ChainClusterParameters specify the options used for the clustering of the
 * chains in structures using the {@link ChainClusterer}.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class ChainClusterParameters implements Serializable {

	private static final long serialVersionUID = 1L;

	private int minimumSequenceLength = 20;
	private int absoluteMinimumSequenceLength = 5;
	// if the shortest sequence length is >= 0.75 * the median sequence length,
	// then the minimum sequence length is set to shortest sequence length,
	// but not shorter than the absoluteMinimumSequenceLength.
	// This adaptive feature allows the consideration of very short chains, such
	// as collagen
	private double minimumSequenceLengthFraction = 0.75;
	private double[] sequenceIdentityThresholds = { 0.0, 0.95 };
	private double sequencePseudoSymmetryThreshold = 0.95;
	private double alignmentFractionThreshold = 0.9;
	private boolean onTheFly = true;
	private boolean verbose = false;

	private static final String n = System.getProperty("line.separator");

	/**
	 * @return the minimumSequenceLength
	 */
	public int getMinimumSequenceLength() {
		return minimumSequenceLength;
	}

	/**
	 * @param minimumSequenceLength
	 *            the minimumSequenceLength to set
	 */
	public void setMinimumSequenceLength(int minimumSequenceLength) {
		this.minimumSequenceLength = minimumSequenceLength;
	}

	/**
	 * @return the absoluteMinimumSequenceLength
	 */
	public int getAbsoluteMinimumSequenceLength() {
		return absoluteMinimumSequenceLength;
	}

	/**
	 * @param absoluteMinimumSequenceLength
	 *            the absoluteMinimumSequenceLength to set
	 */
	public void setAbsoluteMinimumSequenceLength(
			int absoluteMinimumSequenceLength) {
		this.absoluteMinimumSequenceLength = absoluteMinimumSequenceLength;
	}

	/**
	 * @return the minimumSequenceLengthFraction
	 */
	public double getMinimumSequenceLengthFraction() {
		return minimumSequenceLengthFraction;
	}

	/**
	 * @param minimumSequenceLengthFraction
	 *            the minimumSequenceLengthFraction to set
	 */
	public void setMinimumSequenceLengthFraction(
			double minimumSequenceLengthFraction) {
		this.minimumSequenceLengthFraction = minimumSequenceLengthFraction;
	}

	/**
	 * @return the sequenceIdentityThreshold
	 */
	public double[] getSequenceIdentityThresholds() {
		return sequenceIdentityThresholds;
	}

	/**
	 * @param sequenceIdentityThresholds
	 *            the sequenceIdentityThresholds to set
	 */
	public void setSequenceIdentityThresholds(
			double[] sequenceIdentityThresholds) {
		this.sequenceIdentityThresholds = sequenceIdentityThresholds;
	}

	/**
	 * @return the alignmentFractionThreshold
	 */
	public double getAlignmentFractionThreshold() {
		return alignmentFractionThreshold;
	}

	/**
	 * @param alignmentFractionThreshold
	 *            the alignmentFractionThreshold to set
	 */
	public void setAlignmentFractionThreshold(double alignmentFractionThreshold) {
		this.alignmentFractionThreshold = alignmentFractionThreshold;
	}

	public double getSequencePseudoSymmetryThreshold() {
		return sequencePseudoSymmetryThreshold;
	}

	public void setSequencePseudoSymmetryThreshold(
			double sequencePseudoSymmetryThreshold) {
		this.sequencePseudoSymmetryThreshold = sequencePseudoSymmetryThreshold;
	}

	/**
	 * @return true if Jmol on the fly bioassembly generation is used
	 */
	public boolean isOnTheFly() {
		return onTheFly;
	}

	/**
	 * @param onTheFly
	 *            the onTheFly to set
	 */
	public void setOnTheFly(boolean useJmolBioAssemblies) {
		this.onTheFly = useJmolBioAssemblies;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	@Override
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("Minimum protein sequence length   : ");
		s.append(minimumSequenceLength);
		s.append(n);
		s.append("Sequence identity thresholds      : ");
		s.append(Arrays.toString(sequenceIdentityThresholds));
		s.append(n);
		s.append("Sequence pseudosymmetry threshold : ");
		s.append(sequencePseudoSymmetryThreshold);
		s.append(n);
		s.append("Alignment fraction threshold      : ");
		s.append(alignmentFractionThreshold);
		s.append(n);
		s.append("Verbose                           : ");
		s.append(verbose);
		s.append(n);
		return s.toString();
	}
}
