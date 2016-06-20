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
 * The ChainClustererParameters specifies the options used for the clustering of
 * the chains in structures using the {@link ChainClusterer}.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class ChainClustererParameters implements Serializable {

	private static final long serialVersionUID = 1L;

	private int minimumSequenceLength = 20;
	private int absoluteMinimumSequenceLength = 5;
	private double minimumSequenceLengthFraction = 0.75;

	private double sequenceIdentityThreshold = 0.95;
	private double tmScoreThreshold = 0.5;
	private double alignmentFractionThreshold = 0.9;

	private ChainClustererMethod clustererMethod = ChainClustererMethod.SEQUENCE;

	/**
	 * Get the minimum number of residues of a chain to be considered in the
	 * clusters.
	 * 
	 * @return minimumSequenceLength
	 */
	public int getMinimumSequenceLength() {
		return minimumSequenceLength;
	}

	/**
	 * Set the minimum number of residues of a chain to be considered in the
	 * clusters.
	 * 
	 * @param minimumSequenceLength
	 */
	public void setMinimumSequenceLength(int minimumSequenceLength) {
		this.minimumSequenceLength = minimumSequenceLength;
	}

	/**
	 * If the shortest chain sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median chain sequence length,
	 * then the minimumSequenceLength is set to shortest chain sequence length,
	 * but not shorter than the absoluteMinimumSequenceLength.
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
	 * If the shortest chain sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median chain sequence length,
	 * then the minimumSequenceLength is set to shortest chain sequence length,
	 * but not shorter than the absoluteMinimumSequenceLength.
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
	 * If the shortest chain sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median chain sequence length,
	 * then the minimumSequenceLength is set to shortest chain sequence length,
	 * but not shorter than the absoluteMinimumSequenceLength.
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
	 * If the shortest chain sequence length is higher or equal the
	 * minimumSequenceLengthFraction times the median chain sequence length,
	 * then the minimumSequenceLength is set to shortest chain sequence length,
	 * but not shorter than the absoluteMinimumSequenceLength.
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
	 * Sequence identity threshold to consider for the chain clustering.
	 * <p>
	 * Two chains with sequence identity equal or higher than the threshold will
	 * be clustered together.
	 * 
	 * @return sequenceIdentityThreshold
	 */
	public double getSequenceIdentityThreshold() {
		return sequenceIdentityThreshold;
	}

	/**
	 * Sequence identity threshold to consider for the sequence chain
	 * clustering.
	 * <p>
	 * Two chains with sequence identity equal or higher than the threshold will
	 * be clustered together.
	 * 
	 * @param sequenceIdentityThreshold
	 */
	public void setSequenceIdentityThreshold(double sequenceIdentityThreshold) {
		this.sequenceIdentityThreshold = sequenceIdentityThreshold;
	}

	/**
	 * Structure similarity threshold (measured with TM-Score) to consider for
	 * the structural chain clustering.
	 * 
	 * @return tmScoreThreshold
	 */
	public double getTmScoreThreshold() {
		return tmScoreThreshold;
	}

	/**
	 * Structure similarity threshold (measured with TM-Score) to consider for
	 * the structural chain clustering.
	 * 
	 * @param tmScoreThreshold
	 */
	public void setTmScoreThreshold(double tmScoreThreshold) {
		this.tmScoreThreshold = tmScoreThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two chains to be
	 * clustered together.
	 * 
	 * @return alignmentFractionThreshold
	 */
	public double getAlignmentFractionThreshold() {
		return alignmentFractionThreshold;
	}

	/**
	 * The minimum coverage of the sequence alignment between two chains to be
	 * clustered together.
	 * 
	 * @param alignmentFractionThreshold
	 */
	public void setAlignmentFractionThreshold(double alignmentFractionThreshold) {
		this.alignmentFractionThreshold = alignmentFractionThreshold;
	}

	/**
	 * Method to cluster chains.
	 * 
	 * @return ChainClustererMethod
	 */
	public ChainClustererMethod getClustererMethod() {
		return clustererMethod;
	}

	/**
	 * Method to cluster chains.
	 * 
	 * @param ChainClustererMethod
	 */
	public void setClustererMethod(ChainClustererMethod method) {
		this.clustererMethod = method;
	}

	@Override
	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("Minimum chain sequence length  : ");
		s.append(minimumSequenceLength);
		s.append(System.getProperty("line.separator"));
		s.append("Sequence identity threshold    : ");
		s.append(sequenceIdentityThreshold);
		s.append(System.getProperty("line.separator"));
		s.append("TM-Score threshold             : ");
		s.append(tmScoreThreshold);
		s.append(System.getProperty("line.separator"));
		s.append("Alignment fraction threshold   : ");
		s.append(alignmentFractionThreshold);
		s.append(System.getProperty("line.separator"));
		s.append("Chain clusterer method         : ");
		s.append(clustererMethod);
		s.append(System.getProperty("line.separator"));
		return s.toString();
	}
}
