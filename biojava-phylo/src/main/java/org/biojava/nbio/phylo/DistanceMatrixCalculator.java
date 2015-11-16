package org.biojava.nbio.phylo;

import java.util.List;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
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
	private DistanceMatrixCalculator() {}

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
	 * The fractional dissimilarity score (Ds) is the average character
	 * dissimilarity between two aligned sequences. It is calculated as:
	 * 
	 * <pre>
	 * Ds = sum(max(M) - M<sub>ai,bi</sub>) / L
	 * </pre>
	 * 
	 * Where the sum through i runs for all the alignment positions, ai and bi
	 * are the AA at position i in the first and second aligned sequences,
	 * respectively, and L is the total length of the alignment (normalization).
	 * <p>
	 * The fractional dissimilarity score (Ds) is defined in the interval [0,
	 * 1], where 0 means that the sequences are identical and 1 that the
	 * sequences are completely different (similarity score of 0). Note that a
	 * value higher than one can be obtained, if the similarity score of the
	 * alignment is negative (rare case).
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

		int n = msa.getSize();
		DistanceMatrix DM = new BasicSymmetricalDistanceMatrix(n);

		// Calculate the similarity scores using the alignment package
		double[] scores = Alignments.getAllPairsScores(
				msa.getAlignedSequences(),
				PairwiseSequenceScorerType.GLOBAL_SIMILARITIES,
				new SimpleGapPenalty(0, 0), M);

		// The maximum score that can be obtained for this alignment
		int maxScore = M.getMaxValue() * msa.getLength();

		// Convert and copy the dissimilarity scores to the matrix
		int indx = 0;
		for (int i = 0; i < n; i++) {
			DM.setIdentifier(i, msa.getAlignedSequence(i + 1)
					.getAccession().getID());
			
			for (int j = i; j < n; j++) {
				if (i == j)
					DM.setValue(i, j, 0.0);
				else {
					double dS = (maxScore - scores[indx]) / msa.getLength();
					DM.setValue(i, j, dS);
					DM.setValue(j, i, dS);
					indx++;
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
	 * Ds = maxScore - sum(M<sub>ai,bi</sub>)
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
					D from = seqs.get(i).getCompoundAt(k);
					D to = seqs.get(i).getCompoundAt(k);
					if (Comparison.isGap(from.getShortName().charAt(0))
							|| Comparison.isGap(to.getShortName().charAt(0)))
						continue;
					score += M.getValue(from, to);
				}
				DM.setValue(i, j, score);

				if (score > maxscore) {
					maxscore = score;
				}
			}
		}

		for (int i = 0; i < n; i++) {
			DM.setIdentifier(i, msa.getAlignedSequence(i + 1)
					.getAccession().getID());
			
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
		// TODO
		return null;
	}

	/**
	 * The structural distance (d<sub>S</sub>) uses the structural similarity
	 * (or dissimilarity) from a the structural alignment of two protein
	 * strutures. It is based on the diffusive model for protein fold evolution.
	 * The structural deviations are captured as RMS deviations.
	 * 
	 * @param rmsdMat
	 *            RMSD matrix for all structure pairs (symmetric matrix)
	 * @return DistanceMatrix
	 */
	public static <C extends Sequence<D>, D extends Compound> DistanceMatrix structuralDistance(
			double[][] rmsdMat) {
		// TODO
		return null;
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
