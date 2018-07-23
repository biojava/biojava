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
package org.biojava.nbio.ronn;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.data.sequence.FastaSequence;
import org.biojava.nbio.ronn.Jronn.Range;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;


public class JronnTest {

	// Implemented with two platform checks.
	public void verifyRanges() {

	Range[]	ranges = Jronn.getDisorder(new FastaSequence("name", "LLRGRHLMNGTMIMRPWNFLNDHHFPKFFPHLIEQQAIWLADWWRKKHC" +
				"RPLPTRAPTMDQWDHFALIQKHWTANLWFLTFPFNDKWGWIWFLKDWTPGSADQAQRACTWFFCHGHDTN" +
				"CQIIFEGRNAPERADPMWTGGLNKHIIARGHFFQSNKFHFLERKFCEMAEIERPNFTCRTLDCQKFPWDDP" +
				"CSSTHSDCPKLEDLISFTETHGCSAADNADRPSQACHIGWAAMCEPTAMFMLMGSRCRCSFWPQNNAARHR" +
				"NFLIQIEMHSHLEHWIQTLHPQRPFLCNTWDDNWPICQFASQARGNSPDHHP"));
	assertEquals(4, ranges.length);
	assertEquals(53, ranges[0].from);
	assertEquals(59, ranges[0].to);

	assertEquals(190, ranges[1].from);
	assertEquals(196, ranges[1].to);

	assertEquals(210, ranges[2].from);
	assertEquals(226, ranges[2].to);

	assertEquals(305, ranges[3].from);
	assertEquals(313, ranges[3].to);
	//System.out.println(Arrays.toString(ranges));
	}

	/**
	 * Jronn was breaking on Windows platform due to the different System.getProperty("line.separator") values
	 * (CRLF vs LF).  This wraps the existing unit testing to show that it works on windows or unix.
	 */
	@Test
	public void checkJronn() {
		// Windows CRLF
		ScopedProperty lineSepProp = new ScopedProperty("line.separator", "\r\n");
		verifyRanges();
		lineSepProp.close();

		// UNIX LF
		ScopedProperty lineSepPropUnix = new ScopedProperty("line.separator", "\n");
		verifyRanges();
		lineSepPropUnix.close();
	}

	/**
	 * A scoped property helper class to check with Windows style CRLF.
	 * Credit Thomas Klambauer, but here have removed the implement of
	 * AutoCloseable (Java 7) for BioJava support of Java 6.
	 */
	public class ScopedProperty {

		private final String key;
		private final String oldValue;

		/**
		 *
		 * @param key The System.setProperty key
		 * @param value The System.setProperty value to switch to.
		 */
		public ScopedProperty(final String key, final String value) {
			this.key = key;
			oldValue = System.setProperty(key, value);
		}

		public void close() {
			// Can't use setProperty(key, null) -> Throws NullPointerException.
			if( oldValue == null ) {
				// Previously there was no entry.
				System.clearProperty(key);
			} else {
				System.setProperty(key, oldValue);
			}
		}
	}

	/** Test the user scenario when disorder scores are calculated over a sequence containing the selenocysteine (Sec) amino acid.
	 *  The user has to manually convert the stop codons symbols ("*") to the "U" symbol and Jronn is expected to handle this sequence.
	 *
	 * @throws CompoundNotFoundException
	 */

	@Test
	public void testStopCodons() throws CompoundNotFoundException {

		// gene: DIO2, NM_001007023
		String dnaString = "ATGGGCATCCTCAGCGTAGACTTGCTGATCACACTGCAAATTCTGCCAGTTTTTTTCTCCAACTGCCTCT" +
				"TCCTGGCTCTCTATGACTCGGTCATTCTGCTCAAGCACGTGGTGCTGCTGTTGAGCCGCTCCAAGTCCAC" +
				"TCGCGGAGAGTGGCGGCGCATGCTGACCTCAGAGGGACTGCGCTGCGTCTGGAAGAGCTTCCTCCTCGAT" +
				"GCCTACAAACAGCTAAATTGTCCTCCATCAGGTTTTAGCAAAGATGGACACATTTTATGACTAGTATATG" +
				"AAGCTTATAAAAGCAGACTACTGGTCTACTCACATTTGGATTTATGGATGGTGAAATTGGGTGAGGATGC" +
				"CCCCAATTCCAGTGTGGTGCATGTCTCCAGTACAGAAGGAGGTGACAACAGTGGCAATGGTACCCAGGAG" +
				"AAGATAGCTGAGGGAGCCACATGCCACCTTCTTGACTTTGCCAGCCCTGAGCGCCCACTAGTGGTCAACT" +
				"TTGGCTCAGCCACTTGACCTCCTTTCACGAGCCAGCTGCCAGCCTTCCGCAAACTGGTGGAAGAGTTCTC" +
				"CTCAGTGGCTGACTTCCTGCTGGTCTACATTGATGAGGCTCATCCATCAGATGGCTGGGCGATACCGGGG" +
				"GACTCCTCTTTGTCTTTTGAGGTGAAGAAGCACCAGAACCAGGAAGATCGATGTGCAGCAGCCCAGCAGC" +
				"TTCTGGAGCGTTTCTCCTTGCCGCCCCAGTGCCGAGTTGTGGCTGACCGCATGGACAATAACGCCAACAT" +
				"AGCTTACGGGGTAGCCTTTGAACGTGTGTGCATTGTGCAGAGACAGAAAATTGCTTATCTGGGAGGAAAG" +
				"GGCCCCTTCTCCTACAACCTTCAAGAAGTCCGGCATTGGCTGGAGAAGAATTTCAGCAAGAGATGAAAGA" +
				"AAACTAGATTAGCTGGTTAA";

		DNASequence dnaSequence = new DNASequence(dnaString);
		RNASequence mRNA = dnaSequence.getRNASequence();
		ProteinSequence protein = mRNA.getProteinSequence();

		String proteinString = protein.getSequenceAsString();
		// replace the symbol * (codon TGA) with U at amino acid
		String proteinStringU = proteinString.replaceAll("\\*", "u");
		FastaSequence fsequence = new FastaSequence("", proteinStringU);
		try {
			float[] rawProbabilityScores = Jronn.getDisorderScores(fsequence);
		}
		catch (Exception e) {
			fail("Disorder scores calculation doesn't work");
		}
	}
}
