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

public class QuatSymmetryParameters {
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
	private double angleThreshold = 10.0; // max angle deviation for C2 solver
	// if a structure has both cyclic and helical symmetry (i.e., 3J4F: C2 and H), 
	// then helical symmetry is assigned if Rmsd(helical) - Rmsd(cyclic) <= helixRmsdThreshold
	// A slightly positive value gives preference to helical, if the RMSDs for the two symmetries
	// are almost identical
	private double helixRmsdThreshold = 0.05;
	private double helixRmsdToRiseRatio = 0.5; // rmsd must be < 0.5 * abs(rise) (previously set to 0.75)
	private double minimumHelixRise = 1.0;
	private double minimumHelixAngle = 5.0; // min helix angle to differentiate it from a translational repeat
	private int maximumLocalCombinations = 50000; // max number of combinations to try for local symmetry calculation
	private int maximumLocalResults = 1000;
	private int maximumLocalSubunits = 20; // maximum number of subunits for local symmetry calculations
	private boolean localSymmetry = true;
	private double localTimeLimit = 120; // time limit for local calculations in seconds
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
	/**
	 * @return the rmsdThreshold
	 */
	public double getRmsdThreshold() {
		return rmsdThreshold;
	}
	/**
	 * @param rmsdThreshold the rmsdThreshold to set
	 */
	public void setRmsdThreshold(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}
	public double getAngleThreshold() {
		return angleThreshold;
	}
	public void setAngleThreshold(double angleThreshold) {
		this.angleThreshold = angleThreshold;
	}
	
	public double getHelixRmsdThreshold() {
		return helixRmsdThreshold;
	}
	public void setHelixRmsdThreshold(double helixRmsdThreshold) {
		this.helixRmsdThreshold = helixRmsdThreshold;
	}
	/**
	 * @return the helixRmsdToRiseRatio
	 */
	public double getHelixRmsdToRiseRatio() {
		return helixRmsdToRiseRatio;
	}
	/**
	 * @param helixRmsdToRiseRatio the helixRmsdToRiseRatio to set
	 */
	public void setHelixRmsdToRiseRatio(double helixRmsdToRiseRatio) {
		this.helixRmsdToRiseRatio = helixRmsdToRiseRatio;
	}
	public double getMinimumHelixRise() {
		return minimumHelixRise;
	}
	public void setMinimumHelixRise(double minimumHelixRise) {
		this.minimumHelixRise = minimumHelixRise;
	}
	public double getMinimumHelixAngle() {
		return minimumHelixAngle;
	}
	public void setMinimumHelixAngle(double minimumHelixAngle) {
		this.minimumHelixAngle = minimumHelixAngle;
	}
	public double getSequencePseudoSymmetryThreshold() {
		return sequencePseudoSymmetryThreshold;
	}
	
	public void setSequencePseudoSymmetryThreshold(
			double sequencePseudoSymmetryThreshold) {
		this.sequencePseudoSymmetryThreshold = sequencePseudoSymmetryThreshold;
	}
	
	public int getMaximumLocalCombinations() {
		return maximumLocalCombinations;
	}
	public void setMaximumLocalCombinations(int maximumLocalCombinations) {
		this.maximumLocalCombinations = maximumLocalCombinations;
	}
	/**
	 * @return the maximumLocalResults
	 */
	public int getMaximumLocalResults() {
		return maximumLocalResults;
	}
	/**
	 * @return the maximumLocalSubunits
	 */
	public int getMaximumLocalSubunits() {
		return maximumLocalSubunits;
	}
	/**
	 * @param maximumLocalSubunits the maximumLocalSubunits to set
	 */
	public void setMaximumLocalSubunits(int maximumLocalSubunits) {
		this.maximumLocalSubunits = maximumLocalSubunits;
	}
	/**
	 * @param maximumLocalResults the maximumLocalResults to set
	 */
	public void setMaximumLocalResults(int maximumLocalResults) {
		this.maximumLocalResults = maximumLocalResults;
	}
	public boolean isLocalSymmetry() {
		return localSymmetry;
	}
	public void setLocalSymmetry(boolean localSymmetry) {
		this.localSymmetry = localSymmetry;
	}
	/**
	 * @return the localTimeLimit
	 */
	public double getLocalTimeLimit() {
		return localTimeLimit;
	}
	/**
	 * @param localTimeLimit the localTimeLimit to set
	 */
	public void setLocalTimeLimit(double localTimeLimit) {
		this.localTimeLimit = localTimeLimit;
	}
	/**
	 * @return true if Jmol on the fly bioassembly generation is used
	 */
	public boolean isOnTheFly() {
		return onTheFly;
	}
	/**
	 * @param onTheFly the onTheFly to set
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
		s.append("Angle threshold                   : ");
		s.append(angleThreshold);
		s.append(n);	
		s.append("Symmetry RMSD threshold           : ");
		s.append(rmsdThreshold);
		s.append(n);
		s.append("Local symmetry                    : ");
		s.append(localSymmetry);
		s.append(n);
		s.append("Verbose                           : ");
		s.append(verbose);
		s.append(n);
		return s.toString();
	}
}
