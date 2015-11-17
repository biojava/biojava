package org.biojava.nbio.phylo;

import java.util.List;

import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.forester.evoinference.distance.PairwiseDistanceCalculator;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.msa.Msa;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The DistanceMatrixCalculator methods generate a {@link DistanceMatrix} from a
 * {@link MultipleSequenceAlignment}.
 * <p>
 * The implementations differ in the information required to calculate the
 * distances. Thus, the difference resides in their constructor.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class DistanceMatrixCalculator {

	private static final Logger logger = LoggerFactory
			.getLogger(DistanceMatrixCalculator.class);

	/** Prevent instantiation */
	private DistanceMatrixCalculator() {
	}

	/**
	 * The fractional dissimilarity (D) is defined as the percentage of sites
	 * that differ between two aligned sequences. The percentage of identity
	 * (PID) is the fraction of identical sites between two aligned sequences.
	 * 
	 * <pre>
	 * D = 1 - PID
	 * </pre>
	 * 
	 * The gapped positons in the alignment are ignored in the calculation. This
	 * method is a wrapper to the forester implementation of the calculation:
	 * {@link PairwiseDistanceCalculator#calcFractionalDissimilarities(Msa)}
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix fractionalDissimilarity(
			MultipleSequenceAlignment<C, D> msa) {

		Msa fMsa = null;
		try {
			fMsa = ForesterWrapper.convert(msa);
		} catch (Exception e) {
			logger.warn("Could not convert BioJava MSA to forester MSA", e);
		}

		DistanceMatrix DM = PairwiseDistanceCalculator
				.calcFractionalDissimilarities(fMsa);

		return DM;
	}

	/**
	 * The Poisson (correction) evolutionary distance (d) is a function of the
	 * fractional dissimilarity (D), given by:
	 * 
	 * <pre>
	 * d = -log(1 - D)
	 * </pre>
	 * 
	 * The gapped positons in the alignment are ignored in the calculation. This
	 * method is a wrapper to the forester implementation of the calculation:
	 * {@link PairwiseDistanceCalculator#calcPoissonDistances(Msa)}
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix poissonDistance(
			MultipleSequenceAlignment<C, D> msa) {

		Msa fMsa = null;
		try {
			fMsa = ForesterWrapper.convert(msa);
		} catch (Exception e) {
			logger.warn("Could not convert BioJava MSA to forester MSA", e);
		}

		DistanceMatrix DM = PairwiseDistanceCalculator
				.calcPoissonDistances(fMsa);

		return DM;
	}

	/**
	 * The Kimura evolutionary distance (d) is a correction of the fractional
	 * dissimilarity (D) specially needed for large evolutionary distances. It
	 * is given by:
	 * 
	 * <pre>
	 * d = -log(1 - D - 0.2 * D<sup>2</sup>)
	 * </pre>
	 * 
	 * The equation is derived by fitting the relationship between the
	 * evolutionary distance (d) and the fractional dissimilarity (D) according
	 * to the PAM model of evolution (it is an empirical approximation for the
	 * method {@link #pamDistance(MultipleSequenceAlignment}). The gapped
	 * positons in the alignment are ignored in the calculation. This method is
	 * a wrapper to the forester implementation of the calculation:
	 * {@link PairwiseDistanceCalculator#calcKimuraDistances(Msa)}.
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix kimuraDistance(
			MultipleSequenceAlignment<C, D> msa) {

		Msa fMsa = null;
		try {
			fMsa = ForesterWrapper.convert(msa);
		} catch (Exception e) {
			logger.warn("Could not convert BioJava MSA to forester MSA", e);
		}

		DistanceMatrix DM = PairwiseDistanceCalculator
				.calcPoissonDistances(fMsa);

		return DM;
	}

	/**
	 * BioJava implementation for percentage of identity (PID). Although the
	 * name of the method is percentage of identity, the DistanceMatrix contains
	 * the fractional dissimilarity (D), computed as D = 1 - PID.
	 * <p>
	 * It is recommended to use the method
	 * {@link DistanceMatrixCalculator#fractionalDissimilarity(MultipleSequenceAlignment)}
	 * instead of this one.
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix percentageIdentity(
			MultipleSequenceAlignment<C, D> msa) {

		logger.info("{}:{}", "Determing Distances", 0);
		int n = msa.getSize();
		String[] sequenceString = new String[n];
		for (int i = 0; i < n; i++) {
			sequenceString[i] = msa.getAlignedSequence(i + 1)
					.getSequenceAsString();
		}

		DistanceMatrix distance = new BasicSymmetricalDistanceMatrix(n);
		int totalloopcount = (n / 2) * (n + 1);

		int loopcount = 0;
		for (int i = 0; i < (n - 1); i++) {
			logger.info("{}:{}", "Determining Distances", (loopcount * 100)
					/ totalloopcount);
			distance.setIdentifier(i, msa.getAlignedSequence(i + 1)
					.getAccession().getID());

			for (int j = i; j < n; j++) {
				loopcount++;
				if (j == i) {
					distance.setValue(i, j, 0);
				} else {
					distance.setValue(i, j, 100 - Comparison.PID(
							sequenceString[i], sequenceString[j]));
					distance.setValue(j, i, distance.getValue(i, j));
				}
			}
		}
		logger.info("{}:{}", "Determining Distances", 100);

		return distance;
	}

	/**
	 * The fractional dissimilarity score (Ds) is a relative measure of the
	 * dissimilarity between two aligned sequences. It is calculated as:
	 * 
	 * <pre>
	 * Ds = sum( max(M) - M<sub>ai,bi</sub> ) / (max(M)-min(M)) ) / L
	 * </pre>
	 * 
	 * Where the sum through i runs for all the alignment positions, ai and bi
	 * are the AA at position i in the first and second aligned sequences,
	 * respectively, and L is the total length of the alignment (normalization).
	 * <p>
	 * The fractional dissimilarity score (Ds) is in the interval [0, 1], where
	 * 0 means that the sequences are identical and 1 that the sequences are
	 * completely different.
	 * <p>
	 * Gaps do not have a contribution to the similarity score calculation (gap
	 * penalty = 0)
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @param M
	 *            SubstitutionMatrix for similarity scoring
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix fractionalDissimilarityScore(
			MultipleSequenceAlignment<C, D> msa, SubstitutionMatrix<D> M) {

		// Calculate the similarity scores using the alignment package
		logger.info("{}:{}", "Determing Distances", 0);

		int n = msa.getSize();
		DistanceMatrix DM = new BasicSymmetricalDistanceMatrix(n);
		int totalloopcount = (n / 2) * (n + 1);
		int end = msa.getLength();

		String[] sequenceString = new String[n];
		for (int i = 0; i < n; i++) {
			sequenceString[i] = msa.getAlignedSequence(i + 1)
					.getSequenceAsString();
		}
		List<C> seqs = msa.getAlignedSequences();

		int loopcount = 0;
		for (int i = 0; i < (n - 1); i++) {
			logger.info("{}:{}", "Determining Distances", (loopcount * 100)
					/ totalloopcount);

			// Obtain the similarity scores
			for (int j = i; j < n; j++) {

				double score = 0;
				loopcount++;

				for (int k = 0; k < end; k++) {
					if (Comparison.isGap(sequenceString[i].charAt(k))
							|| Comparison.isGap(sequenceString[j].charAt(k)))
						continue;
					score += M.getValue(seqs.get(i).getCompoundAt(k + 1), seqs
							.get(j).getCompoundAt(k + 1));
				}

				if (i == j)
					DM.setValue(i, j, 0.0);
				else {
					double dS = (M.getMaxValue() - score / msa.getLength())
							/ (M.getMaxValue() - M.getMinValue());
					
					DM.setValue(i, j, dS);
					DM.setValue(j, i, dS);
				}
			}
		}
		return DM;
	}

	/**
	 * The dissimilarity score is the additive inverse of the similarity score
	 * (sum of scores) between two aligned sequences using a substitution model
	 * (Substitution Matrix). The maximum dissimilarity score is taken to be the
	 * maximum similarity score between self-alignments (each sequence against
	 * itself). Calculation of the score is as follows:
	 * 
	 * <pre>
	 * Ds = maxScore - sum<sub>i</sub>(M<sub>ai,bi</sub>)
	 * </pre>
	 * 
	 * It is recommended to use the method
	 * {@link #fractionalDissimilarityScore(MultipleSequenceAlignment, SubstitutionMatrix)}
	 * , since the maximum similarity score is not relative to the data set, but
	 * relative to the Substitution Matrix, and the score is normalized by the
	 * alignment length (fractional).
	 * <p>
	 * Gaps do not have a contribution to the similarity score calculation (gap
	 * penalty = 0).
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @param M
	 *            SubstitutionMatrix for similarity scoring
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix dissimilarityScore(
			MultipleSequenceAlignment<C, D> msa, SubstitutionMatrix<D> M) {

		logger.info("{}:{}", "Determing Distances", 0);

		int n = msa.getSize();
		String[] sequenceString = new String[n];
		for (int i = 0; i < n; i++) {
			sequenceString[i] = msa.getAlignedSequence(i + 1)
					.getSequenceAsString();
		}
		List<C> seqs = msa.getAlignedSequences();

		DistanceMatrix DM = new BasicSymmetricalDistanceMatrix(n);
		int totalloopcount = (n / 2) * (n + 1);

		double maxscore = 0;
		int end = msa.getLength();
		int loopcount = 0;

		for (int i = 0; i < (n - 1); i++) {
			logger.info("{}:{}", "Determining Distances", (loopcount * 100)
					/ totalloopcount);

			// Obtain the similarity scores
			for (int j = i; j < n; j++) {

				double score = 0;
				loopcount++;

				for (int k = 0; k < end; k++) {
					if (Comparison.isGap(sequenceString[i].charAt(k))
							|| Comparison.isGap(sequenceString[j].charAt(k)))
						continue;
					score += M.getValue(seqs.get(i).getCompoundAt(k + 1), seqs
							.get(j).getCompoundAt(k + 1));
				}

				if (i != j)
					DM.setValue(i, j, score);

				if (score > maxscore) {
					maxscore = score;
				}
			}
		}

		for (int i = 0; i < n; i++) {
			DM.setIdentifier(i, msa.getAlignedSequence(i + 1).getAccession()
					.getID());

			for (int j = i; j < n; j++) {
				if (i == j)
					DM.setValue(i, j, 0.0);
				else {
					double dS = maxscore - DM.getValue(i, j);
					DM.setValue(i, j, dS);
					DM.setValue(j, i, dS);
				}
			}
		}

		logger.info("{}:{}", "Determining Distances", 100);
		return DM;
	}

	/**
	 * The PAM (Point Accepted Mutations) distance is a measure of evolutionary
	 * distance in protein sequences. The PAM unit represents an average
	 * substitution rate of 1% per site. The fractional dissimilarity (D) of two
	 * aligned sequences is related with the PAM distance (d) by the equation:
	 * 
	 * <pre>
	 * D = sum(fi * (1 - M<sub>ii</sub><sup>d</sup>))
	 * </pre>
	 * 
	 * Where the sum is for all 20 AA, fi denotes the natural fraction of the
	 * given AA and M is the substitution matrix (in this case the PAM1 matrix).
	 * <p>
	 * To calculate the PAM distance between two aligned sequences the maximum
	 * likelihood (ML) approach is used, which consists in finding d that
	 * maximazies the function:
	 * 
	 * <pre>
	 * L(d) = product(f<sub>ai</sub> * (1 - M<sub>ai,bi</sub><sup>d</sup>))
	 * </pre>
	 * 
	 * Where the product is for every position i in the alignment, and ai and bi
	 * are the AA at position i in the first and second aligned sequences,
	 * respectively.
	 * 
	 * @param msa
	 *            MultipleSequenceAlignment
	 * @return
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix pamDistance(
			MultipleSequenceAlignment<C, D> msa) {

		// Need to import PAM1 matrix to biojava TODO
		SubstitutionMatrix<AminoAcidCompound> PAM1 = SubstitutionMatrixHelper
				.getPAM250();

		return null;
	}

	/**
	 * The structural distance (d<sub>S</sub>) uses the structural similarity
	 * (or dissimilarity) from a the structural alignment of two protein
	 * strutures. It is based on the diffusive model for protein fold evolution
	 * (Grishin 1995). The structural deviations are captured as RMS deviations.
	 * 
	 * <pre>
	 * d<sub>Sij</sub> = (rmsd<sub>max</sub><sup>2</sup> / alpha<sup>2</sup>) * 
	 *        ln( (rmsd<sub>max</sub><sup>2</sup> - rmsd<sub>0</sub><sup>2</sup>) / 
	 *        (rmsd<sub>max</sub><sup>2</sup> - (rmsd<sub>ij</sub><sup>2</sup>) )
	 * </pre>
	 * 
	 * @param rmsdMat
	 *            RMSD matrix for all structure pairs (symmetric matrix)
	 * @param alpha
	 *            change in CA positions introduced by a single AA substitution
	 *            (Grishin 1995: 1 A)
	 * @param rmsdMax
	 *            estimated RMSD between proteins of the same fold when the
	 *            percentage of identity is infinitely low (the maximum allowed
	 *            RMSD of proteins with the same fold). (Grishin 1995: 5 A)
	 * @param rmsd0
	 *            arithmetical mean of squares of the RMSD for identical
	 *            proteins (Grishin 1995: 0.4 A)
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix structuralDistance(
			double[][] rmsdMat, double alpha, double rmsdMax, double rmsd0) {

		int n = rmsdMat.length;
		DistanceMatrix DM = new BasicSymmetricalDistanceMatrix(n);

		// Transform raw RMSD into distances and create matrix
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (i == j)
					DM.setValue(i, j, 0.0);
				else {
					double d = (rmsdMax * rmsdMax)
							/ (alpha * alpha)
							* Math.log((rmsdMax * rmsdMax - rmsd0 * rmsd0)
									/ (rmsdMax * rmsdMax - rmsdMat[i][j]
											* rmsdMat[i][j]));
					DM.setValue(i, j, d);
					DM.setValue(j, i, d);
				}
			}
		}

		return DM;
	}

	/**
	 * The joint sequence-structure distance (d<sub>SS</sub>) is a combination
	 * of the sequence-based and the structure-based distances.
	 * 
	 * @param rmsdMat
	 *            RMSD matrix for all structure pairs (symmetric matrix)
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix jointSeqStrucDistance(
			double[][] rmsdMat) {
		// TODO
		return null;
	}

}
