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

import org.biojava.nbio.alignment.PairwiseStatisticsAligner.EncodedSequence;
import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.Assert.*;

public class PairwiseStatisticsAlignerTest {

	private static final String AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";
	private static final String NUCLEOTIDES = "ACGT";

	private static final SubstitutionMatrix<AminoAcidCompound> BLOSUM62 = SubstitutionMatrixHelper.getBlosum62();
	private static final SubstitutionMatrix<AminoAcidCompound> PAM250 = SubstitutionMatrixHelper.getPAM250();

	private static GapPenalty[] gapPenalties() {
		return new GapPenalty[] {
				new SimpleGapPenalty(11, 1), // affine
				new SimpleGapPenalty(5, 3), // affine, extension close to opening
				new SimpleGapPenalty(0, 2), // linear
				new SimpleGapPenalty(7, 0), // constant
		};
	}

	/**
	 * The statistics of the alignment that the traceback of SmithWaterman or NeedlemanWunsch returns
	 */
	private static <S extends Sequence<C>, C extends Compound> PairwiseAlignmentStatistics expected(S query, S target,
			GapPenalty gaps, SubstitutionMatrix<C> matrix, boolean local) {
		AbstractPairwiseSequenceAligner<S, C> aligner = local ? new SmithWaterman<>(query, target, gaps, matrix)
				: new NeedlemanWunsch<>(query, target, gaps, matrix);
		SequencePair<S, C> pair = aligner.getPair();
		int length = pair.getLength();
		int score = (int) aligner.getScore();
		if (!local) {
			return new PairwiseAlignmentStatistics(score, pair.getNumIdenticals(), length, 1, query.getLength(), 1,
					target.getLength());
		}
		if (length == 0) {
			return new PairwiseAlignmentStatistics(score, 0, 0, 0, 0, 0, 0);
		}
		return new PairwiseAlignmentStatistics(score, pair.getNumIdenticals(), length, pair.getIndexInQueryAt(1),
				pair.getIndexInQueryAt(length), pair.getIndexInTargetAt(1), pair.getIndexInTargetAt(length));
	}

	private static String randomSequence(Random random, int length, String alphabet) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < length; i++) {
			s.append(alphabet.charAt(random.nextInt(alphabet.length())));
		}
		return s.toString();
	}

	/**
	 * Substitutions, insertions and deletions, each at the given rate
	 */
	private static String mutate(Random random, String seq, double rate, String alphabet) {
		StringBuilder s = new StringBuilder();
		for (char c : seq.toCharArray()) {
			double p = random.nextDouble();
			if (p < rate) {
				s.append(alphabet.charAt(random.nextInt(alphabet.length())));
			} else if (p < 2 * rate) {
				// deleted
			} else if (p < 3 * rate) {
				s.append(c).append(randomSequence(random, 1 + random.nextInt(4), alphabet));
			} else {
				s.append(c);
			}
		}
		return s.toString();
	}

	/**
	 * Families of related sequences, from very close to very far, and unrelated ones
	 */
	private static List<String> sequences(Random random, String alphabet, int numFamilies, int minLength) {
		List<String> seqs = new ArrayList<>();
		for (int family = 0; family < numFamilies; family++) {
			String parent = randomSequence(random, minLength + 30 * family, alphabet);
			seqs.add(parent);
			for (double rate : new double[] {0.01, 0.05, 0.1, 0.2, 0.35}) {
				seqs.add(mutate(random, parent, rate, alphabet));
			}
			// a fragment of the parent
			seqs.add(parent.substring(family, parent.length() / 2 + family));
		}
		return seqs;
	}

	private static <S extends Sequence<C>, C extends Compound> void checkAllPairs(List<S> seqs,
			SubstitutionMatrix<C> matrix) {
		for (GapPenalty gaps : gapPenalties()) {
			for (boolean local : new boolean[] {true, false}) {
				PairwiseStatisticsAligner<C> aligner = new PairwiseStatisticsAligner<>(gaps, matrix, local);
				List<EncodedSequence<C>> encoded = seqs.stream().map(aligner::encode).collect(Collectors.toList());
				for (int i = 0; i < seqs.size(); i++) {
					for (int j = 0; j < seqs.size(); j++) {
						PairwiseAlignmentStatistics expected = expected(seqs.get(i), seqs.get(j), gaps, matrix, local);
						String message = (local ? "local" : "global") + " " + gaps.getType() + " "
								+ gaps.getOpenPenalty() + "/" + gaps.getExtensionPenalty() + ", sequences " + i + " and "
								+ j;
						assertEquals(message, expected, aligner.align(encoded.get(i), encoded.get(j)));
					}
				}
			}
		}
	}

	private static List<ProteinSequence> proteins(List<String> seqs) throws CompoundNotFoundException {
		List<ProteinSequence> proteins = new ArrayList<>();
		for (String s : seqs) {
			proteins.add(new ProteinSequence(s));
		}
		return proteins;
	}

	@Test
	public void testProteinsBlosum62() throws CompoundNotFoundException {
		List<String> seqs = sequences(new Random(12345), AMINO_ACIDS, 3, 30);
		// with letters other than the 20 standard amino acids, which BLOSUM62 scores or not
		seqs.add(seqs.get(0).substring(0, 10) + "XBZUOJ" + seqs.get(0).substring(10));
		checkAllPairs(proteins(seqs), BLOSUM62);
	}

	@Test
	public void testProteinsPam250() throws CompoundNotFoundException {
		checkAllPairs(proteins(sequences(new Random(54321), AMINO_ACIDS, 3, 25)), PAM250);
	}

	@Test
	public void testDna() throws CompoundNotFoundException {
		List<DNASequence> seqs = new ArrayList<>();
		for (String s : sequences(new Random(999), NUCLEOTIDES, 3, 30)) {
			seqs.add(new DNASequence(s));
		}
		checkAllPairs(seqs, SubstitutionMatrixHelper.getNuc4_4());
	}

	/**
	 * Low complexity sequences have many alignments with the same score: the statistics must still be those of
	 * the alignment of the traceback
	 */
	@Test
	public void testTies() throws CompoundNotFoundException {
		List<String> seqs = sequences(new Random(4242), "AG", 2, 20);
		seqs.add("AAAAAAAAAA");
		seqs.add("AAAAGAAAAA");
		seqs.add("GAGAGAGAGAGA");
		checkAllPairs(proteins(seqs), BLOSUM62);
	}

	@Test
	public void testSingleCompounds() throws CompoundNotFoundException {
		List<String> seqs = new ArrayList<>();
		for (char c : "AWPC".toCharArray()) {
			seqs.add(String.valueOf(c));
		}
		seqs.add("WP");
		checkAllPairs(proteins(seqs), BLOSUM62);
	}

	@Test
	public void testValues() throws CompoundNotFoundException {
		PairwiseStatisticsAligner<AminoAcidCompound> aligner = new PairwiseStatisticsAligner<>(
				new SimpleGapPenalty(11, 1), BLOSUM62, true);
		String seq = "MTKPVLQDGETVVCQGTHAAWYRSFNE";
		int n = seq.length();

		PairwiseAlignmentStatistics stats = aligner.align(new ProteinSequence(seq), new ProteinSequence(seq));
		assertEquals(n, stats.getNumIdenticals());
		assertEquals(n, stats.getLength());
		assertEquals(1, stats.getQueryStart());
		assertEquals(n, stats.getQueryEnd());
		assertEquals(1, stats.getTargetStart());
		assertEquals(n, stats.getTargetEnd());

		// the target within the query, with other residues around
		stats = aligner.align(new ProteinSequence("GGGGG" + seq + "GGG"), new ProteinSequence(seq));
		assertEquals(n, stats.getNumIdenticals());
		assertEquals(n, stats.getLength());
		assertEquals(6, stats.getQueryStart());
		assertEquals(5 + n, stats.getQueryEnd());
		assertEquals(1, stats.getTargetStart());
		assertEquals(n, stats.getTargetEnd());

		// 2 residues deleted in the middle of the target: 2 gap columns, costing 11 + 2 * 1
		int selfScore = aligner.align(new ProteinSequence(seq), new ProteinSequence(seq)).getScore();
		int deletedScore = BLOSUM62.getValue(aa(seq.charAt(13)), aa(seq.charAt(13)))
				+ BLOSUM62.getValue(aa(seq.charAt(14)), aa(seq.charAt(14)));
		stats = aligner.align(new ProteinSequence(seq), new ProteinSequence(seq.substring(0, 13) + seq.substring(15)));
		assertEquals(selfScore - deletedScore - 13, stats.getScore());
		assertEquals(n - 2, stats.getNumIdenticals());
		assertEquals(n, stats.getLength());
		assertEquals(1, stats.getQueryStart());
		assertEquals(n, stats.getQueryEnd());
		assertEquals(1, stats.getTargetStart());
		assertEquals(n - 2, stats.getTargetEnd());
	}

	private static AminoAcidCompound aa(char c) {
		return BLOSUM62.getCompoundSet().getCompoundForString(String.valueOf(c));
	}

	/**
	 * Nothing scores above zero, or there is nothing to align: the local alignment is empty. The linear gap case
	 * used to make SmithWaterman throw an exception.
	 */
	@Test
	public void testEmptyLocalAlignments() throws CompoundNotFoundException {
		PairwiseAlignmentStatistics empty = new PairwiseAlignmentStatistics(0, 0, 0, 0, 0, 0, 0);
		for (GapPenalty gaps : gapPenalties()) {
			PairwiseStatisticsAligner<AminoAcidCompound> aligner = new PairwiseStatisticsAligner<>(gaps, BLOSUM62, true);
			assertEquals(empty, aligner.align(new ProteinSequence("W"), new ProteinSequence("P")));
			assertEquals(empty, aligner.align(new ProteinSequence("GGGGGGGG"), new ProteinSequence("WWWWWWWW")));
			assertEquals(empty, aligner.align(new ProteinSequence(""), new ProteinSequence("ACDEF")));
			assertEquals(empty, aligner.align(new ProteinSequence("ACDEF"), new ProteinSequence("")));
			assertEquals(empty, aligner.align(new ProteinSequence(""), new ProteinSequence("")));
		}
	}

	@Test
	public void testEmptyGlobalAlignments() throws CompoundNotFoundException {
		PairwiseStatisticsAligner<AminoAcidCompound> affine = new PairwiseStatisticsAligner<>(
				new SimpleGapPenalty(11, 1), BLOSUM62, false);
		assertEquals(new PairwiseAlignmentStatistics(-16, 0, 5, 0, 0, 1, 5),
				affine.align(new ProteinSequence(""), new ProteinSequence("ACDEF")));
		assertEquals(new PairwiseAlignmentStatistics(-16, 0, 5, 1, 5, 0, 0),
				affine.align(new ProteinSequence("ACDEF"), new ProteinSequence("")));
		assertEquals(new PairwiseAlignmentStatistics(0, 0, 0, 0, 0, 0, 0),
				affine.align(new ProteinSequence(""), new ProteinSequence("")));

		PairwiseStatisticsAligner<AminoAcidCompound> linear = new PairwiseStatisticsAligner<>(
				new SimpleGapPenalty(0, 2), BLOSUM62, false);
		assertEquals(new PairwiseAlignmentStatistics(-10, 0, 5, 0, 0, 1, 5),
				linear.align(new ProteinSequence(""), new ProteinSequence("ACDEF")));
		assertEquals(new PairwiseAlignmentStatistics(-10, 0, 5, 1, 5, 0, 0),
				linear.align(new ProteinSequence("ACDEF"), new ProteinSequence("")));
	}

	@Test
	public void testLinearGapNoPositiveScore() throws CompoundNotFoundException {
		// the case that made SmithWaterman throw an exception: both must now give an empty alignment
		ProteinSequence w = new ProteinSequence("W"), p = new ProteinSequence("P");
		GapPenalty linear = new SimpleGapPenalty(0, 2);
		assertEquals(expected(w, p, linear, BLOSUM62, true),
				new PairwiseStatisticsAligner<>(linear, BLOSUM62, true).align(w, p));
	}

	@Test
	public void testThreads() throws CompoundNotFoundException {
		List<ProteinSequence> seqs = proteins(sequences(new Random(777), AMINO_ACIDS, 4, 40));
		PairwiseStatisticsAligner<AminoAcidCompound> aligner = new PairwiseStatisticsAligner<>(
				new SimpleGapPenalty(11, 1), BLOSUM62, true);
		List<EncodedSequence<AminoAcidCompound>> encoded = seqs.stream().map(aligner::encode)
				.collect(Collectors.toList());
		int n = seqs.size();
		List<PairwiseAlignmentStatistics> sequential = IntStream.range(0, n * n)
				.mapToObj(k -> aligner.align(encoded.get(k / n), encoded.get(k % n))).collect(Collectors.toList());
		List<PairwiseAlignmentStatistics> parallel = IntStream.range(0, n * n).parallel()
				.mapToObj(k -> aligner.align(encoded.get(k / n), encoded.get(k % n))).collect(Collectors.toList());
		assertEquals(sequential, parallel);
	}

	@Test(expected = IllegalArgumentException.class)
	public void testSequenceEncodedByOtherAligner() throws CompoundNotFoundException {
		GapPenalty gaps = new SimpleGapPenalty(11, 1);
		PairwiseStatisticsAligner<AminoAcidCompound> a1 = new PairwiseStatisticsAligner<>(gaps, BLOSUM62, true);
		PairwiseStatisticsAligner<AminoAcidCompound> a2 = new PairwiseStatisticsAligner<>(gaps, BLOSUM62, true);
		ProteinSequence s = new ProteinSequence("ACDEF");
		a1.align(a1.encode(s), a2.encode(s));
	}

	@Test(expected = IllegalArgumentException.class)
	public void testPositiveGapPenalty() {
		GapPenalty gaps = new SimpleGapPenalty(11, 1);
		// SimpleGapPenalty makes its penalties negative, but other implementations may not
		GapPenalty positive = new GapPenalty() {
			@Override public int getExtensionPenalty() { return 1; }
			@Override public int getOpenPenalty() { return gaps.getOpenPenalty(); }
			@Override public Type getType() { return Type.AFFINE; }
			@Override public void setExtensionPenalty(int gep) { }
			@Override public void setOpenPenalty(int gop) { }
		};
		new PairwiseStatisticsAligner<>(positive, BLOSUM62, true);
	}

	@Test
	public void testGapPenaltyIsCopied() throws CompoundNotFoundException {
		SimpleGapPenalty gaps = new SimpleGapPenalty(11, 1);
		PairwiseStatisticsAligner<AminoAcidCompound> aligner = new PairwiseStatisticsAligner<>(gaps, BLOSUM62, true);
		ProteinSequence q = new ProteinSequence("MTKPVLQDGETVVCQGTHAAWYRSFNE"),
				t = new ProteinSequence("MTKPVLQDGEAAWYRSFNE");
		PairwiseAlignmentStatistics before = aligner.align(q, t);
		gaps.setOpenPenalty(1);
		assertEquals(before, aligner.align(q, t));
	}
}
