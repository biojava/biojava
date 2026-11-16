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

package org.biojava.nbio.alignment;

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

/**
 * {@link AllPairsScores} against the scorers it replaces: {@link NeedlemanWunsch} and {@link SmithWaterman} run
 * through {@link Alignments#getAllPairsScorers} and {@link Alignments#runPairwiseScorers}.
 */
public class AllPairsScoresTest {

	private static final String AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV";

	/** Gap penalties as (open, extension); an open penalty of 0 makes a LINEAR SimpleGapPenalty. */
	private static final int[][] GAPS = { { 10, 1 }, { 2, 1 }, { 11, 2 }, { 0, 1 }, { 0, 4 }, { 1, 0 } };

	private static List<ProteinSequence> randomProteins(long seed, int count, int minLength, int maxLength,
			String alphabet) throws CompoundNotFoundException {
		Random random = new Random(seed);
		List<ProteinSequence> proteins = new ArrayList<>();
		for (int i = 0; i < count; i++) {
			int length = minLength + random.nextInt(maxLength - minLength + 1);
			StringBuilder sb = new StringBuilder(length);
			for (int p = 0; p < length; p++) {
				sb.append(alphabet.charAt(random.nextInt(alphabet.length())));
			}
			proteins.add(new ProteinSequence(sb.toString()));
		}
		return proteins;
	}

	/** Mutated copies of one sequence, so that the pairs are related (high scores, long local alignments). */
	private static List<ProteinSequence> family(long seed, int count, int length) throws CompoundNotFoundException {
		Random random = new Random(seed);
		String root = randomProteins(seed, 1, length, length, AMINO_ACIDS).get(0).getSequenceAsString();
		List<ProteinSequence> proteins = new ArrayList<>();
		for (int i = 0; i < count; i++) {
			StringBuilder sb = new StringBuilder();
			for (char c : root.toCharArray()) {
				double r = random.nextDouble();
				if (r < 0.05) {
					continue; // deletion
				} else if (r < 0.10) {
					sb.append(AMINO_ACIDS.charAt(random.nextInt(AMINO_ACIDS.length()))); // insertion
				} else if (r < 0.25) {
					c = AMINO_ACIDS.charAt(random.nextInt(AMINO_ACIDS.length())); // substitution
				}
				sb.append(c);
			}
			proteins.add(new ProteinSequence(sb.toString()));
		}
		return proteins;
	}

	private static <S extends Sequence<C>, C extends Compound> void assertSameAsAligners(List<S> sequences,
			GapPenalty gaps, SubstitutionMatrix<C> matrix) {
		for (PairwiseSequenceScorerType type : new PairwiseSequenceScorerType[] { PairwiseSequenceScorerType.GLOBAL,
				PairwiseSequenceScorerType.LOCAL }) {
			assertTrue(AllPairsScores.isSupported(sequences, type, gaps, matrix));
			List<PairwiseSequenceScorer<S, C>> expected = Alignments.getAllPairsScorers(sequences, type, gaps, matrix);
			String what = type + " " + gaps.getType() + " gaps " + gaps.getOpenPenalty() + "/" + gaps.getExtensionPenalty();
			double[] scores = Alignments.getAllPairsScores(sequences, type, gaps, matrix);
			List<PairwiseSequenceScorer<S, C>> actual = AllPairsScores.getScorers(sequences, type, gaps, matrix);
			assertEquals(what, expected.size(), scores.length);
			assertEquals(what, expected.size(), actual.size());
			for (int k = 0; k < expected.size(); k++) {
				PairwiseSequenceScorer<S, C> e = expected.get(k);
				PairwiseSequenceScorer<S, C> a = actual.get(k);
				assertSame(e.getQuery(), a.getQuery());
				assertSame(e.getTarget(), a.getTarget());
				double expectedScore;
				try {
					expectedScore = e.getScore();
				} catch (ArrayIndexOutOfBoundsException noPositiveCell) {
					// SmithWaterman with a linear gap penalty fails on a pair with no positive-scoring cell, whose
					// local score is 0 (see noPositiveLocalCell)
					assertEquals(what, PairwiseSequenceScorerType.LOCAL, type);
					assertEquals(what, GapPenalty.Type.LINEAR, gaps.getType());
					assertEquals(what + " pair " + k, 0.0, scores[k], 0.0);
					assertEquals(what + " pair " + k, 0.0, a.getScore(), 0.0);
					continue;
				}
				assertEquals(what + " pair " + k, expectedScore, scores[k], 0.0);
				assertEquals(what + " pair " + k, expectedScore, a.getScore(), 0.0);
				assertEquals(what + " pair " + k, e.getMaxScore(), a.getMaxScore(), 0.0);
				assertEquals(what + " pair " + k, e.getMinScore(), a.getMinScore(), 0.0);
				assertEquals(what + " pair " + k, e.getDistance(), a.getDistance(), 0.0);
				assertEquals(what + " pair " + k, e.getSimilarity(), a.getSimilarity(), 0.0);
			}
		}
	}

	/**
	 * A local alignment of sequences with no positive substitution score between them is empty, with score 0.
	 * {@link SmithWaterman} gives 0 with an affine gap penalty, but fails with a linear one; the scores computed
	 * here are 0 for both.
	 */
	@Test
	public void noPositiveLocalCell() throws CompoundNotFoundException {
		SubstitutionMatrix<AminoAcidCompound> blosum62 = SubstitutionMatrixHelper.getBlosum62();
		List<ProteinSequence> pair = Arrays.asList(new ProteinSequence("W"), new ProteinSequence("MK"));
		GapPenalty affine = new SimpleGapPenalty(10, 1);
		GapPenalty linear = new SimpleGapPenalty(0, 1);
		assertEquals(0.0, new SmithWaterman<>(pair.get(0), pair.get(1), affine, blosum62).getScore(), 0.0);
		assertArrayEquals(new double[] { 0.0 }, Alignments.getAllPairsScores(pair, PairwiseSequenceScorerType.LOCAL, affine, blosum62), 0.0);
		assertArrayEquals(new double[] { 0.0 }, Alignments.getAllPairsScores(pair, PairwiseSequenceScorerType.LOCAL, linear, blosum62), 0.0);
	}

	@Test
	public void randomProteinsWithBlosum62() throws CompoundNotFoundException {
		SubstitutionMatrix<AminoAcidCompound> blosum62 = SubstitutionMatrixHelper.getBlosum62();
		List<ProteinSequence> proteins = randomProteins(42, 24, 1, 120, AMINO_ACIDS);
		for (int[] g : GAPS) {
			assertSameAsAligners(proteins, new SimpleGapPenalty(g[0], g[1]), blosum62);
		}
	}

	@Test
	public void relatedProteinsWithOtherMatrices() throws CompoundNotFoundException {
		List<ProteinSequence> proteins = family(7, 16, 150);
		for (SubstitutionMatrix<AminoAcidCompound> matrix : Arrays.asList(SubstitutionMatrixHelper.getBlosum62(),
				SubstitutionMatrixHelper.getPAM250(), SubstitutionMatrixHelper.getBlosum30())) {
			for (int[] g : GAPS) {
				assertSameAsAligners(proteins, new SimpleGapPenalty(g[0], g[1]), matrix);
			}
		}
	}

	/** Compounds outside the 20 amino acids (X, B, Z, U, stop) score as the matrix's getValue scores them. */
	@Test
	public void ambiguousAndRareCompounds() throws CompoundNotFoundException {
		List<ProteinSequence> proteins = randomProteins(3, 12, 5, 60, AMINO_ACIDS + "XBZU*");
		assertSameAsAligners(proteins, new SimpleGapPenalty(10, 1), SubstitutionMatrixHelper.getBlosum62());
	}

	@Test
	public void identicalAndVeryDifferentLengths() throws CompoundNotFoundException {
		List<ProteinSequence> proteins = new ArrayList<>(randomProteins(5, 3, 300, 300, AMINO_ACIDS));
		proteins.add(new ProteinSequence(proteins.get(0).getSequenceAsString()));
		proteins.add(new ProteinSequence("W"));
		proteins.add(new ProteinSequence("MK"));
		for (int[] g : GAPS) {
			assertSameAsAligners(proteins, new SimpleGapPenalty(g[0], g[1]), SubstitutionMatrixHelper.getBlosum62());
		}
	}

	@Test
	public void dnaWithNuc44() throws CompoundNotFoundException {
		Random random = new Random(11);
		List<DNASequence> dna = new ArrayList<>();
		for (int i = 0; i < 16; i++) {
			StringBuilder sb = new StringBuilder();
			for (int p = 0, n = 20 + random.nextInt(200); p < n; p++) {
				sb.append("ACGT".charAt(random.nextInt(4)));
			}
			dna.add(new DNASequence(sb.toString()));
		}
		SubstitutionMatrix<NucleotideCompound> nuc44 = SubstitutionMatrixHelper.getNuc4_4();
		for (int[] g : GAPS) {
			assertSameAsAligners(dna, new SimpleGapPenalty(g[0], g[1]), nuc44);
		}
	}

	@Test
	public void unsupportedInputsTakeTheAlignerPath() throws CompoundNotFoundException {
		SubstitutionMatrix<AminoAcidCompound> blosum62 = SubstitutionMatrixHelper.getBlosum62();
		GapPenalty gaps = new SimpleGapPenalty();
		List<ProteinSequence> proteins = randomProteins(1, 4, 5, 20, AMINO_ACIDS);
		assertFalse(AllPairsScores.isSupported(proteins, PairwiseSequenceScorerType.GLOBAL_IDENTITIES, gaps, blosum62));
		assertFalse(AllPairsScores.isSupported(proteins, PairwiseSequenceScorerType.KMERS, gaps, blosum62));
		assertFalse(AllPairsScores.isSupported(proteins.subList(0, 1), PairwiseSequenceScorerType.GLOBAL, gaps, blosum62));
		assertArrayEquals(Alignments.runPairwiseScorers(Alignments.getAllPairsScorers(proteins,
				PairwiseSequenceScorerType.GLOBAL_IDENTITIES, gaps, blosum62)),
				Alignments.getAllPairsScores(proteins, PairwiseSequenceScorerType.GLOBAL_IDENTITIES, gaps, blosum62), 0.0);
	}

	/** The guide tree, and so the multiple alignment, are those the aligners' scores give. */
	@Test
	public void multipleSequenceAlignmentIsUnchanged() throws CompoundNotFoundException {
		List<ProteinSequence> proteins = family(13, 12, 80);
		SubstitutionMatrix<AminoAcidCompound> blosum62 = SubstitutionMatrixHelper.getBlosum62();
		for (PairwiseSequenceScorerType type : new PairwiseSequenceScorerType[] { PairwiseSequenceScorerType.GLOBAL,
				PairwiseSequenceScorerType.LOCAL }) {
			GapPenalty gaps = new SimpleGapPenalty();
			List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> scorers = Alignments.getAllPairsScorers(
					proteins, type, gaps, blosum62);
			Alignments.runPairwiseScorers(scorers);
			GuideTree<ProteinSequence, AminoAcidCompound> tree = new GuideTree<>(proteins, scorers);
			Profile<ProteinSequence, AminoAcidCompound> expected = Alignments.getProgressiveAlignment(tree,
					ProfileProfileAlignerType.GLOBAL, gaps, blosum62);

			Profile<ProteinSequence, AminoAcidCompound> actual = Alignments.getMultipleSequenceAlignment(proteins, type,
					gaps, blosum62, ProfileProfileAlignerType.GLOBAL);
			assertEquals(type.toString(), expected.toString(), actual.toString());
		}
	}
}
