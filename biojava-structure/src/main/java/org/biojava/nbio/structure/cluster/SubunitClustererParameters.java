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

import org.biojava.nbio.structure.align.ce.CeMain;

import java.io.Serializable;

/**
 * The SubunitClustererParameters specifies the options used for the clustering
 * of the subunits in structures using the {@link SubunitClusterer}.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class SubunitClustererParameters implements Serializable {

	private static final long serialVersionUID = 1L;

	private int minimumSequenceLength = 20;
	private int absoluteMinimumSequenceLength = 5;
	private double minimumSequenceLengthFraction = 0.75;

	private boolean useGlobalMetrics;
	private double sequenceIdentityThreshold;
	private double sequenceCoverageThreshold = 0.75;

	private double rmsdThreshold = 3.0;
	private double structureCoverageThreshold = 0.75;
	private double tmThreshold = 0.5;

	private SubunitClustererMethod clustererMethod = SubunitClustererMethod.SEQUENCE_STRUCTURE;

	private String superpositionAlgorithm = CeMain.algorithmName;
	private boolean optimizeAlignment = true;

	private boolean useSequenceCoverage;
	private boolean useRMSD;
	private boolean useStructureCoverage;
	private boolean useTMScore;

	private boolean internalSymmetry = false;

	/**
	 * Subunits aligned with these or better scores will be considered "identical".
	 */
	private static final double hcSequenceIdentityLocal = 0.95;
	private static final double hcSequenceCoverageLocal = 0.75;
	private static final double hcSequenceIdentityGlobal = 0.85;

	/**
	 * "Local" metrics are scoring
	 * SubunitClustererMethod.SEQUENCE: sequence identity of a local alignment
	 *                                  (normalised by the number of aligned residues)
	 *                                  sequence coverage of the alignment
	 *                                  (normalised by the length of the longer sequence)
	 * SubunitClustererMethod.STRUCTURE: RMSD of the aligned substructures
	 *                                   and structure coverage of the alignment
	 *                                   (normalised by the length of the larger structure)
	 * Two thresholds for each method are required.
	 *
	 * "Global" metrics are scoring
	 * SubunitClustererMethod.SEQUENCE: sequence identity of a global alignment
	 *                                  (normalised by the length of the alignment)
	 * SubunitClustererMethod.STRUCTURE: TMScore of the aligned structures
	 *                                  (normalised by the length of the larger structure)
	 * One threshold for each method is required.
	 *
	 */
	public SubunitClustererParameters(boolean useGlobalMetrics) {
		this.useGlobalMetrics = useGlobalMetrics;

		if (useGlobalMetrics) {
			sequenceIdentityThreshold = hcSequenceIdentityGlobal;
			useSequenceCoverage = false;
			useRMSD = false;
			useStructureCoverage = false;
			useTMScore = true;
		} else {
			sequenceIdentityThreshold = hcSequenceIdentityLocal;
			useSequenceCoverage = true;
			useRMSD = true;
			useStructureCoverage = true;
			useTMScore = false;
		}
	}

	/**
	 * Initialize with "local" metrics by default.
	 */
	public SubunitClustererParameters() {
		this(false);
	}

	/**
	 * Get the minimum number of residues of a subunits to be considered in the
	 * clusters.
	 * 
	 * @return minimumSequenceLength
	 */
	public int getMinimumSequenceLength() {
		return minimumSequenceLength;
	}

	/**
	 * Set the minimum number of residues of a subunits to be considered in the
	 * clusters.
	 * 
	 * @param minimumSequenceLength
	 */
	public void setMinimumSequenceLength(int minimumSequenceLength) {
		this.minimumSequenceLength = minimumSequenceLength;
	}

	/**
	 * If the shortest subunit sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median subunit sequence length,
	 * then the minimumSequenceLength is set to shortest subunit sequence
	 * length, but not shorter than the absoluteMinimumSequenceLength.
	 * <p>
	 * This adaptive feature allows the consideration of structures mainly
	 * constructed by very short chains, such as collagen (1A3I)
	 * 
	 * @return the absoluteMinimumSequenceLength
	 */
	public int getAbsoluteMinimumSequenceLength() {
		return absoluteMinimumSequenceLength;
	}

	/**
	 * If the shortest subunit sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median subunit sequence length,
	 * then the minimumSequenceLength is set to shortest subunit sequence
	 * length, but not shorter than the absoluteMinimumSequenceLength.
	 * <p>
	 * This adaptive feature allows the consideration of structures mainly
	 * constructed by very short chains, such as collagen (1A3I)
	 * 
	 * @param absoluteMinimumSequenceLength
	 */
	public void setAbsoluteMinimumSequenceLength(
			int absoluteMinimumSequenceLength) {
		this.absoluteMinimumSequenceLength = absoluteMinimumSequenceLength;
	}

	/**
	 * If the shortest subunit sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median subunit sequence length,
	 * then the minimumSequenceLength is set to shortest subunit sequence
	 * length, but not shorter than the absoluteMinimumSequenceLength.
	 * <p>
	 * This adaptive feature allows the consideration of structures mainly
	 * constructed by very short chains, such as collagen (1A3I)
	 * 
	 * @return the minimumSequenceLengthFraction
	 */
	public double getMinimumSequenceLengthFraction() {
		return minimumSequenceLengthFraction;
	}

	/**
	 * If the shortest subunit sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median subunit sequence length,
	 * then the minimumSequenceLength is set to shortest subunit sequence
	 * length, but not shorter than the absoluteMinimumSequenceLength.
	 * <p>
	 * This adaptive feature allows the consideration of structures mainly
	 * constructed by very short chains, such as collagen (1A3I)
	 * 
	 * @param minimumSequenceLengthFraction
	 */
	public void setMinimumSequenceLengthFraction(
			double minimumSequenceLengthFraction) {
		this.minimumSequenceLengthFraction = minimumSequenceLengthFraction;
	}

	/**
	 * Sequence identity threshold to consider for the subunits clustering.
	 * <p>
	 * Two subunits with sequence identity equal or higher than the threshold
	 * will be clustered together.
	 * 
	 * @return sequenceIdentityThreshold
	 */
	public double getSequenceIdentityThreshold() {
		return sequenceIdentityThreshold;
	}

	/**
	 * Sequence identity threshold to consider for the sequence subunit
	 * clustering.
	 * <p>
	 * Two subunits with sequence identity equal or higher than the threshold
	 * will be clustered together.
	 * 
	 * @param sequenceIdentityThreshold
	 */
	public void setSequenceIdentityThreshold(double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two subunits to be
	 * clustered together.
	 *
	 * @return sequenceCoverageThreshold
	 */
	public double getSequenceCoverageThreshold() {
		return sequenceCoverageThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two subunits to be
	 * clustered together.
	 *
	 * @param sequenceCoverageThreshold
	 */
	public void setSequenceCoverageThreshold(double sequenceCoverageThreshold) {
		this.sequenceCoverageThreshold = sequenceCoverageThreshold;
	}

	/**
	 * Structure similarity threshold (measured with RMSD) to consider for the
	 * structural subunit clustering.
	 * 
	 * @return rmsdThreshold
	 */
	public double getRMSDThreshold() {
		return rmsdThreshold;
	}

	/**
	 * Structure similarity threshold (measured with RMSD) to consider for the
	 * structural subunit clustering.
	 * 
	 * @param rmsdThreshold
	 */
	public void setRMSDThreshold(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}

	/**
	 * Structure similarity threshold (measured with TMScore) to consider for the
	 * structural subunit clustering.
	 *
	 * @return tmThreshold
	 */
	public double getTMThreshold() {
		return tmThreshold;
	}

	/**
	 * Structure similarity threshold (measured with TMScore) to consider for the
	 * structural subunit clustering.
	 *
	 * @param tmThreshold
	 */
	public void setTMThreshold(double tmThreshold) {
		this.tmThreshold = tmThreshold;
	}

	/**
	 * The minimum coverage of the structure alignment between two subunits to be
	 * clustered together.
	 *
	 * @return structureCoverageThreshold
	 */
	public double getStructureCoverageThreshold() {
		return structureCoverageThreshold;
	}

	/**
	 * The minimum coverage of the structure alignment between two subunits to be
	 * clustered together.
	 *
	 * @param structureCoverageThreshold
	 */
	public void setStructureCoverageThreshold(double structureCoverageThreshold) {
		this.structureCoverageThreshold = structureCoverageThreshold;
	}

	/**
	 * Method to cluster subunits.
	 *
	 * @return clustererMethod
	 */
	public SubunitClustererMethod getClustererMethod() {
		return clustererMethod;
	}

	/**
	 * Method to cluster subunits.
	 * 
	 * @param method
	 */
	public void setClustererMethod(SubunitClustererMethod method) {
		this.clustererMethod = method;
	}

	/**
	 * The internal symmetry option divides each {@link Subunit} of each
	 * {@link SubunitCluster} into its internally symmetric repeats.
	 * <p>
	 * The {@link SubunitClustererMethod#STRUCTURE} must be chosen to consider
	 * internal symmetry, otherwise this parameter will be ignored.
	 * 
	 * @return true if internal symmetry is considered, false otherwise
	 */
	public boolean isInternalSymmetry() {
		return internalSymmetry;
	}

	/**
	 * The internal symmetry option divides each {@link Subunit} of each
	 * {@link SubunitCluster} into its internally symmetric repeats.
	 * <p>
	 * The {@link SubunitClustererMethod#STRUCTURE} must be chosen to consider
	 * internal symmetry, otherwise this parameter will be ignored.
	 * 
	 * @param internalSymmetry
	 *            true if internal symmetry is considered, false otherwise
	 */
	public void setInternalSymmetry(boolean internalSymmetry) {
		this.internalSymmetry = internalSymmetry;
	}

	@Override
	public String toString() {
		return "SubunitClustererParameters [minimumSequenceLength="
				+ minimumSequenceLength + ", absoluteMinimumSequenceLength="
				+ absoluteMinimumSequenceLength
				+ ", minimumSequenceLengthFraction="
				+ minimumSequenceLengthFraction
				+ ", sequenceIdentityThreshold=" + sequenceIdentityThreshold
				+ ", rmsdThreshold=" + rmsdThreshold + ", coverageThreshold="
				+ sequenceCoverageThreshold + ", clustererMethod=" + clustererMethod
				+ ", internalSymmetry=" + internalSymmetry + "]";
	}

	/**
	 * Method to superpose subunits (i.e., structural aligner).
	 *
	 * @return superpositionAlgorithm
	 */
	public String getSuperpositionAlgorithm() {
		return superpositionAlgorithm;
	}

	/**
	 * Method to cluster subunits.
	 *
	 * @param superpositionAlgorithm
	 */
	public void setSuperpositionAlgorithm(String superpositionAlgorithm) {
		this.superpositionAlgorithm = superpositionAlgorithm;
	}

	/**
	 * Whether the alignment algorithm should try its best to optimize the alignment,
	 * or we are happy with a quick and dirty result. Effect depends on implementation
	 * of the specific algorithm's method.	 *
	 *
	 * @return optimizeAlignment
	 */
	public boolean isOptimizeAlignment() {
		return optimizeAlignment;
	}

	/**
	 * Whether the alignment algorithm should try its best to optimize the alignment,
	 * or we are happy with a quick and dirty result. Effect depends on implementation
	 * of the specific algorithm's method.	 *
	 *
	 * @param optimizeAlignment
	 */
	public void setOptimizeAlignment(boolean optimizeAlignment) {
		this.optimizeAlignment = optimizeAlignment;
	}

	/**
	 * Use RMSD for evaluating structure similarity
	 *
	 * @return useRMSD
	 */
	public boolean isUseRMSD() { return useRMSD; }

	/**
	 * Use RMSD for evaluating structure similarity
	 *
	 * @param useRMSD
	 */
	public void setUseRMSD(boolean useRMSD) {
		this.useRMSD = useRMSD;
	}

	/**
	 * Use TMScore for evaluating structure similarity
	 *
	 * @return useTMScore
	 */
	public boolean isUseTMScore() {
		return useTMScore;
	}

	/**
	 * Use TMScore for evaluating structure similarity
	 *
	 * @param useTMScore
	 */
	public void setUseTMScore(boolean useTMScore) {
		this.useTMScore = useTMScore;
	}

	/**
	 * Use sequence coverage for evaluating sequence similarity
	 *
	 * @return useSequenceCoverage
	 */
	public boolean isUseSequenceCoverage() {
		return useSequenceCoverage;
	}

	/**
	 * Use sequence coverage for evaluating sequence similarity
	 *
	 * @param useSequenceCoverage
	 */
	public void setUseSequenceCoverage(boolean useSequenceCoverage) {
		this.useSequenceCoverage = useSequenceCoverage;
	}

	/**
	 * Use structure coverage for evaluating sequence similarity
	 *
	 * @return useStructureCoverage
	 */
	public boolean isUseStructureCoverage() {
		return useStructureCoverage;
	}

	/**
	 * Use structure coverage for evaluating sequence similarity
	 *
	 * @param useStructureCoverage
	 */
	public void setUseStructureCoverage(boolean useStructureCoverage) {
		this.useStructureCoverage = useStructureCoverage;
	}

	/**
	 * Use metrics calculated relative to the whole sequence or structure,
	 * rather than the aligned part only
	 *
	 * @return useGlobalMetrics
	 */
	public boolean isUseGlobalMetrics() {
		return useGlobalMetrics;
	}

	/**
	 * Use metrics calculated relative to the whole sequence or structure,
	 * rather than the aligned part only
	 *
	 * @param useGlobalMetrics
	 */
	public void setUseGlobalMetrics(boolean useGlobalMetrics) {
		this.useGlobalMetrics = useGlobalMetrics;
	}

	/**
	 * Whether the subunits can be considered "identical" by sequence alignment.
	 * For local sequence alignment (normalized by the number of aligned pairs)
	 * this means 0.95 or higher identity and 0.75 or higher coverage.
	 * For global sequence alignment (normalised by the alignment length)
	 * this means 0.85 or higher sequence identity.
	 *
	 * @param sequenceIdentity
	 * @param sequenceCoverage
	 * @return true if the sequence alignment scores are equal to
	 * or better than the "high confidence" scores, false otherwise.
	 */
	public boolean isHighConfidenceScores(double sequenceIdentity, double sequenceCoverage) {
		if (useGlobalMetrics)
			return sequenceIdentity>=hcSequenceIdentityGlobal;
		else
			return sequenceIdentity>=hcSequenceIdentityLocal && sequenceCoverage >= hcSequenceCoverageLocal;
	}


}
