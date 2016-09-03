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

	private double sequenceIdentityThreshold = 0.95;
	private double rmsdThreshold = 3.0;
	private double coverageThreshold = 0.75;

	private SubunitClustererMethod clustererMethod = SubunitClustererMethod.STRUCTURE;
	private boolean internalSymmetry = false;

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
	 * Structure similarity threshold (measured with RMSD) to consider for the
	 * structural subunit clustering.
	 * 
	 * @return rmsdThreshold
	 */
	public double getRmsdThreshold() {
		return rmsdThreshold;
	}

	/**
	 * Structure similarity threshold (measured with RMSD) to consider for the
	 * structural subunit clustering.
	 * 
	 * @param rmsdThreshold
	 */
	public void setRmsdThreshold(double rmsdThreshold) {
		this.rmsdThreshold = rmsdThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two subunits to be
	 * clustered together.
	 * 
	 * @return coverageThreshold
	 */
	public double getCoverageThreshold() {
		return coverageThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two subunits to be
	 * clustered together.
	 * 
	 * @param coverageThreshold
	 */
	public void setCoverageThreshold(double coverageThreshold) {
		this.coverageThreshold = coverageThreshold;
	}

	/**
	 * Method to cluster subunits.
	 * 
	 * @return ChainClustererMethod
	 */
	public SubunitClustererMethod getClustererMethod() {
		return clustererMethod;
	}

	/**
	 * Method to cluster subunits.
	 * 
	 * @param SubunitClustererMethod
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
				+ coverageThreshold + ", clustererMethod=" + clustererMethod
				+ ", internalSymmetry=" + internalSymmetry + "]";
	}

}
