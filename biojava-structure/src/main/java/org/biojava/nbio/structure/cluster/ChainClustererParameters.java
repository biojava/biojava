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

/**
 * Created only to resolve compiling errors. TODO has to be removed.
 */
@Deprecated
public class ChainClustererParameters implements Serializable {

	private static final long serialVersionUID = 1L;

	private int minimumSequenceLength = 20;
	private int absoluteMinimumSequenceLength = 5;
	// if the shortest sequence length is >= 0.75 * the median sequence length,
	// then the minimum sequence length is set to shortest sequence length,
	// but not shorter than the absoluteMinimumSequenceLength.
	// This adaptive feature allows the consideration of very short chains, such as collagen
	private double minimumSequenceLengthFraction = 0.75;
	private double[] sequenceIdentityThresholds = {0.0, 0.95};
	private double sequencePseudoSymmetryThreshold = 0.95;
	private double alignmentFractionThreshold = 0.9;
	private double rmsdThreshold = 7.0;

	public double getRmsdThreshold() {
		return rmsdThreshold;
	}
	public void setRmsdThreshold(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}
	/**
	 * @return the minimumSequenceLength
	 */
	public int getMinimumSequenceLength() {
		return minimumSequenceLength;
	}
	/**
	 * @param minimumSequenceLength the minimumSequenceLength to set
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
	 * @param absoluteMinimumSequenceLength the absoluteMinimumSequenceLength to set
	 */
	public void setAbsoluteMinimumSequenceLength(int absoluteMinimumSequenceLength) {
		this.absoluteMinimumSequenceLength = absoluteMinimumSequenceLength;
	}
	/**
	 * @return the minimumSequenceLengthFraction
	 */
	public double getMinimumSequenceLengthFraction() {
		return minimumSequenceLengthFraction;
	}
	/**
	 * @param minimumSequenceLengthFraction the minimumSequenceLengthFraction to set
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
	 * @param sequenceIdentityThresholds the sequenceIdentityThresholds to set
	 */
	public void setSequenceIdentityThresholds(double[] sequenceIdentityThresholds) {
		this.sequenceIdentityThresholds = sequenceIdentityThresholds;
	}

	/**
	 * @return the alignmentFractionThreshold
	 */
	public double getAlignmentFractionThreshold() {
		return alignmentFractionThreshold;
	}
	/**
	 * @param alignmentFractionThreshold the alignmentFractionThreshold to set
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

}
