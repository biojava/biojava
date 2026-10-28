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
package org.biojava.nbio.tornadovm.align;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.GuideTree;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.tornadovm.TestStructures;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class TornadoAlignmentsTest {

	private static final SubstitutionMatrix<AminoAcidCompound> BLOSUM62 = SubstitutionMatrixHelper.getBlosum62();

	@AfterAll
	static void stopBioJavaThreadPool() {
		ConcurrencyTools.shutdown();
	}

	static List<ProteinSequence> random(int n, int minLen, int maxLen, String alphabet, long seed) throws Exception {
		Random rnd = new Random(seed);
		List<ProteinSequence> seqs = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			int len = minLen + rnd.nextInt(maxLen - minLen + 1);
			StringBuilder sb = new StringBuilder(len);
			for (int k = 0; k < len; k++) {
				sb.append(alphabet.charAt(rnd.nextInt(alphabet.length())));
			}
			seqs.add(new ProteinSequence(sb.toString()));
		}
		return seqs;
	}

	@ParameterizedTest
	@CsvSource({ "GLOBAL,10,1", "LOCAL,10,1", "GLOBAL,8,0", "LOCAL,8,0", "GLOBAL,12,2", "LOCAL,3,3", "GLOBAL,0,1", "LOCAL,0,1", "GLOBAL,0,4", "LOCAL,0,3" })
	void scoresMatchBioJava(String type, int gop, int gep) throws Exception {
		// ambiguity codes and U/O exercise compounds outside the 20 standard amino acids
		List<ProteinSequence> seqs = random(40, 1, 300, "ACDEFGHIKLMNPQRSTVWYXBZUO", 7);
		SimpleGapPenalty gaps = new SimpleGapPenalty(gop, gep);
		PairwiseSequenceScorerType t = PairwiseSequenceScorerType.valueOf(type);
		double[] expected = Alignments.getAllPairsScores(seqs, t, gaps, BLOSUM62);
		double[] actual = TestStructures.forced(() -> TornadoAlignments.getAllPairsScores(seqs, t, gaps, BLOSUM62));
		assertArrayEquals(expected, actual, 0.0);
	}

	@Test
	void identicalAndVeryDifferentLengths() throws Exception {
		List<ProteinSequence> seqs = new ArrayList<>(random(5, 1, 3, "ACDEFGHIKLMNPQRSTVWY", 1));
		seqs.addAll(random(3, 900, 1200, "ACDEFGHIKLMNPQRSTVWY", 2));
		seqs.add(new ProteinSequence(seqs.get(6).getSequenceAsString()));
		for (PairwiseSequenceScorerType t : new PairwiseSequenceScorerType[] { PairwiseSequenceScorerType.GLOBAL,
				PairwiseSequenceScorerType.LOCAL }) {
			double[] expected = Alignments.getAllPairsScores(seqs, t, new SimpleGapPenalty(), BLOSUM62);
			double[] actual = TestStructures.forced(
					() -> TornadoAlignments.getAllPairsScores(seqs, t, new SimpleGapPenalty(), BLOSUM62));
			assertArrayEquals(expected, actual, 0.0, t.toString());
		}
	}

	@Test
	void dnaWithNuc44() throws Exception {
		Random rnd = new Random(3);
		List<DNASequence> seqs = new ArrayList<>();
		for (int i = 0; i < 30; i++) {
			StringBuilder sb = new StringBuilder();
			for (int k = 0, len = 50 + rnd.nextInt(400); k < len; k++) {
				sb.append("ACGTN".charAt(rnd.nextInt(5)));
			}
			seqs.add(new DNASequence(sb.toString()));
		}
		SubstitutionMatrix<NucleotideCompound> nuc = SubstitutionMatrixHelper.getNuc4_4();
		for (PairwiseSequenceScorerType t : new PairwiseSequenceScorerType[] { PairwiseSequenceScorerType.GLOBAL,
				PairwiseSequenceScorerType.LOCAL }) {
			double[] expected = Alignments.getAllPairsScores(seqs, t, new SimpleGapPenalty(), nuc);
			double[] actual = TestStructures.forced(
					() -> TornadoAlignments.getAllPairsScores(seqs, t, new SimpleGapPenalty(), nuc));
			assertArrayEquals(expected, actual, 0.0, t.toString());
		}
	}

	@Test
	void guideTreeOfAPfamFamilyIsIdentical() throws Exception {
		List<ProteinSequence> seqs;
		try (InputStream in = getClass().getResourceAsStream("/PF00104_small.fasta")) {
			// the file is an alignment: use the ungapped sequences (BioJava scores '-' as a residue, which gives
			// scores above the aligners' maximum and negative distances), and only the first 60 of them, since
			// BioJava's own scoring of the full family takes ~10 minutes
			seqs = new ArrayList<>();
			for (ProteinSequence s : FastaReaderHelper.readFastaProteinSequence(in).values()) {
				if (seqs.size() < 60) {
					seqs.add(new ProteinSequence(s.getSequenceAsString().replace("-", "")));
				}
			}
		}
		SimpleGapPenalty gaps = new SimpleGapPenalty();
		List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> expected = Alignments.getAllPairsScorers(seqs,
				PairwiseSequenceScorerType.GLOBAL, gaps, BLOSUM62);
		Alignments.runPairwiseScorers(expected);
		List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> actual = TestStructures.forced(
				() -> TornadoAlignments.getAllPairsScorers(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, BLOSUM62));

		assertEquals(expected.size(), actual.size());
		for (int k = 0; k < expected.size(); k++) {
			assertEquals(expected.get(k).getScore(), actual.get(k).getScore(), 0.0);
			assertEquals(expected.get(k).getMaxScore(), actual.get(k).getMaxScore(), 0.0);
			assertEquals(expected.get(k).getMinScore(), actual.get(k).getMinScore(), 0.0);
			assertEquals(expected.get(k).getDistance(), actual.get(k).getDistance(), 0.0);
		}
		assertEquals(new GuideTree<>(seqs, expected).toString(), new GuideTree<>(seqs, actual).toString());
	}
}
