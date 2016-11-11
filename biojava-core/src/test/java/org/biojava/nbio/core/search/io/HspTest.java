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
package org.biojava.nbio.core.search.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.function.Function;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.search.io.blast.BlastHspBuilder;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Paolo Pavan
 */

public class HspTest {
	//private static final Logger logger = LoggerFactory.getLogger(HspTest.class);

	public static Function<String, DNASequence> buildDNASeq = (seq) -> {
		try {
			return new DNASequence(seq);
		} catch (CompoundNotFoundException ex) {
			fail(ex.getMessage());
			return null;
		}
	};


	Hsp<DNASequence, NucleotideCompound> hspImpl;
	Hsp<DNASequence, NucleotideCompound> uncompleteHsp;

	public HspTest() {
	}

	@BeforeClass
	public static void setUpClass() {
	}

	@AfterClass
	public static void tearDownClass() {
	}

	@Before
	public void setUp() {
		try {
			hspImpl = new BlastHspBuilder()
			.setHspNum(1)
			.setHspBitScore(377.211)
			.setHspEvalue(8.04143e-093)
			.setHspQueryFrom(1)
			.setHspQueryTo(224)
			.setHspHitFrom(1035)
			.setHspHitTo(811)
			.setHspQueryFrame(-1)
			.setHspIdentity(213)
			.setHspPositive(213)
			.setHspGaps(5)
			.setHspAlignLen(227)
			.setHspQseq("CTGACGACAGCCATGCACCACCTGTCTCGACTTTCCCCCGAAGGGCACCTAATGTATCTCTACCTCGTTAGTCGGATGTCAAGACCTGGTAAGGTTTTTTCGCGTATCTTCGAATTAAACCACATACTCCACTGCTTGTGCGG-CCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGCCGTACTCCC-AGGTGGA-TACTTATTGTGTTAACTCCGGCACGGAAGG")
			.setHspHseq("CTGACGACAACCATGCACCACCTGTCTCAACTTTCCCC-GAAGGGCACCTAATGTATCTCTACTTCGTTAGTTGGATGTCAAGACCTGGTAAGGTT-CTTCGCGTTGCTTCGAATTAAACCACATACTCCACTGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGTGGATTACTTATTGTGTTAACTCCGGCACAGAAGG")
			.setHspIdentityString("||||||||| |||||||||||||||||| ||||||||| |||||||||||||||||||||||| |||||||| |||||||||||||||||||||||  |||||||  |||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||| ||||||||| ||||||| |||||||||||||||||||||||| |||||")
			.createBlastHsp( buildDNASeq );
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		uncompleteHsp = new BlastHspBuilder()
		.setPercentageIdentity(100.00/100)
		.setHspAlignLen(48)
		.setMismatchCount(0)
		.setHspGaps(0)
		.setHspQueryFrom(1)
		.setHspQueryTo(48)
		.setHspHitFrom(344)
		.setHspHitTo(391)
		.setHspEvalue(4e-19)
		.setHspBitScore(95.6)
		.createBlastHsp( buildDNASeq );

	}

	@After
	public void tearDown() {
	}

	/**
	 * Test of hashCode method, of class Hsp.
	 * @throws CompoundNotFoundException 
	 */
	@Test
	public void testHashCode() {
		System.out.println("hashCode");
		Hsp<DNASequence, NucleotideCompound> instance;
		int expResult;
		int result;

		instance = new BlastHspBuilder()
				.setHspNum(1)
				.setHspBitScore(377.211)
				.setHspEvalue(8.04143e-093)
				.setHspQueryFrom(1)
				.setHspQueryTo(224)
				.setHspHitFrom(1035)
				.setHspHitTo(811)
				.setHspQueryFrame(-1)
				.setHspIdentity(213)
				.setHspPositive(213)
				.setHspGaps(5)
				.setHspAlignLen(227)
				.setHspQseq("CTGACGACAGCCATGCACCACCTGTCTCGACTTTCCCCCGAAGGGCACCTAATGTATCTCTACCTCGTTAGTCGGATGTCAAGACCTGGTAAGGTTTTTTCGCGTATCTTCGAATTAAACCACATACTCCACTGCTTGTGCGG-CCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGCCGTACTCCC-AGGTGGA-TACTTATTGTGTTAACTCCGGCACGGAAGG")
				.setHspHseq("CTGACGACAACCATGCACCACCTGTCTCAACTTTCCCC-GAAGGGCACCTAATGTATCTCTACTTCGTTAGTTGGATGTCAAGACCTGGTAAGGTT-CTTCGCGTTGCTTCGAATTAAACCACATACTCCACTGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGTGGATTACTTATTGTGTTAACTCCGGCACAGAAGG")
				.setHspIdentityString("||||||||| |||||||||||||||||| ||||||||| |||||||||||||||||||||||| |||||||| |||||||||||||||||||||||  |||||||  |||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||| ||||||||| ||||||| |||||||||||||||||||||||| |||||")
				.createBlastHsp( buildDNASeq );

		expResult = hspImpl.hashCode();
		result = instance.hashCode();
		assertEquals(expResult, result);

		instance = new BlastHspBuilder()
				.setPercentageIdentity(100.00/100)
				.setHspAlignLen(48)
				.setMismatchCount(0)
				.setHspGaps(0)
				.setHspQueryFrom(1)
				.setHspQueryTo(48)
				.setHspHitFrom(344)
				.setHspHitTo(391)
				.setHspEvalue(4e-19)
				.setHspBitScore(95.6)
				.createBlastHsp( buildDNASeq );

		expResult = uncompleteHsp.hashCode();
		result = instance.hashCode();
		assertEquals(expResult, result);

		Hsp<DNASequence, NucleotideCompound> uncompleteHsp2 = new BlastHspBuilder()
				.setPercentageIdentity(100.00/100)
				.setHspAlignLen(48)
				.setMismatchCount(0)
				.setHspGaps(0)
				.setHspQueryFrom(1)
				.setHspQueryTo(48)
				.setHspHitFrom(344)
				.setHspHitTo(391)
				.setHspEvalue(4e-19)
				.setHspBitScore(95.6)
				.createBlastHsp( buildDNASeq );

		assertEquals(uncompleteHsp.hashCode(), uncompleteHsp2.hashCode());
	}

	/**
	 * Test of equals method, of class Hsp.
	 * @throws CompoundNotFoundException 
	 */
	@Test
	public void testEquals() {
		System.out.println("equals");
		Hsp<DNASequence, NucleotideCompound> o;
		o = new BlastHspBuilder()
				.setHspNum(1)
				.setHspBitScore(377.211)
				.setHspEvalue(8.04143e-093)
				.setHspQueryFrom(1)
				.setHspQueryTo(224)
				.setHspHitFrom(1035)
				.setHspHitTo(811)
				.setHspQueryFrame(-1)
				.setHspIdentity(213)
				.setHspPositive(213)
				.setHspGaps(5)
				.setHspAlignLen(227)
				.setHspQseq("CTGACGACAGCCATGCACCACCTGTCTCGACTTTCCCCCGAAGGGCACCTAATGTATCTCTACCTCGTTAGTCGGATGTCAAGACCTGGTAAGGTTTTTTCGCGTATCTTCGAATTAAACCACATACTCCACTGCTTGTGCGG-CCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGCCGTACTCCC-AGGTGGA-TACTTATTGTGTTAACTCCGGCACGGAAGG")
				.setHspHseq("CTGACGACAACCATGCACCACCTGTCTCAACTTTCCCC-GAAGGGCACCTAATGTATCTCTACTTCGTTAGTTGGATGTCAAGACCTGGTAAGGTT-CTTCGCGTTGCTTCGAATTAAACCACATACTCCACTGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGTGGATTACTTATTGTGTTAACTCCGGCACAGAAGG")
				.setHspIdentityString("||||||||| |||||||||||||||||| ||||||||| |||||||||||||||||||||||| |||||||| |||||||||||||||||||||||  |||||||  |||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||||||||||| ||||||||| ||||||| |||||||||||||||||||||||| |||||")
				.createBlastHsp( buildDNASeq );
		Hsp<DNASequence, NucleotideCompound> instance = hspImpl;

		assertEquals(o, instance);

		// example of Hsp retrieved from uncomplete report.
		// (Those HSP may come from a tabular format, for example)
		o = new BlastHspBuilder()
				.setPercentageIdentity(100.00/100)
				.setHspAlignLen(48)
				.setMismatchCount(0)
				.setHspGaps(0)
				.setHspQueryFrom(1)
				.setHspQueryTo(48)
				.setHspHitFrom(344)
				.setHspHitTo(391)
				.setHspEvalue(4e-19)
				.setHspBitScore(95.6)
				.createBlastHsp( buildDNASeq );

		assertEquals(uncompleteHsp, o);
	}
}
