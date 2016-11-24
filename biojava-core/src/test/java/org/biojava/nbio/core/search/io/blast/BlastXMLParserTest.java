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
package org.biojava.nbio.core.search.io.blast;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.search.io.Hit;
import org.biojava.nbio.core.search.io.Hsp;
import org.biojava.nbio.core.search.io.Result;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Paolo Pavan
 */
public class BlastXMLParserTest {
	private static final Logger logger = LoggerFactory.getLogger(BlastXMLParserTest.class);

	protected  File getFileForResource(String resource){
		URL resourceURL = this.getClass().getResource(resource);
		String filepath = resourceURL.getFile();
		filepath = filepath.replaceAll("%20"," ");

		File file = new File(filepath);

		return file;
	}

	/**
	 * Test of createObjects method, of class BlastXMLParser.
	 */
	@Test
	public void testCreateObjects() throws Exception {
		System.out.println("createObjects");

		String resource = "/org/biojava/nbio/core/search/io/blast/small-blastreport.blastxml";


		File file = getFileForResource(resource);
		
//		Function<String, Sequence<?>> buildSeq = (seq) -> SearchIO.getSequence(seq);
		Function<String, DNASequence> buildSeq = (seq) -> {
			try {
				return new DNASequence(seq);
			} catch (CompoundNotFoundException ex) {
				logger.error("Unexpected error, could not find compound when creating Sequence object from Hsp", ex);
				return null;
			}
		};
		List<Result<DNASequence, NucleotideCompound>> results = getResults(file,buildSeq);


		// test with random manual selected results
		BlastHsp<DNASequence, NucleotideCompound> hsp1hit1res1 = new BlastHspBuilder()
				.setHspNum(1)
				.setHspBitScore(2894.82)
				.setHspScore(1567)
				.setHspEvalue(0)
				.setHspQueryFrom(1)
				.setHspQueryTo(1567)
				.setHspHitFrom(616309)
				.setHspHitTo(617875)
				.setHspQueryFrame(1)
				.setHspHitFrame(1)
				.setHspIdentity(1567)
				.setHspPositive(1567)
				.setHspGaps(0)
				.setHspAlignLen(1567)
				.setHspQseq("TTAAATTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGTGGCGTGCCTAATACATGCAAGTCGTACGCTAGCCGCTGAATTGATCCTTCGGGTGAAGTGAGGCAATGACTAGAGTGGCGAACTGGTGAGTAACACGTAAGAAACCTGCCCTTTAGTGGGGGATAACATTTGGAAACAGATGCTAATACCGCGTAACAACAAATCACACATGTGATCTGTTTGAAAGGTCCTTTTGGATCGCTAGAGGATGGTCTTGCGGCGTATTAGCTTGTTGGTAGGGTAGAAGCCTACCAAGGCAATGATGCGTAGCCGAGTTGAGAGACTGGCCGGCCACATTGGGACTGAGACACTGCCCAAACTCCTACGGGAGGCTGCAGTAGGGAATTTTCCGCAATGCACGAAAGTGTGACGGAGCGACGCCGCGTGTGTGATGAAGGCTTTCGGGTCGTAAAGCACTGTTGTAAGGGAAGAATAACTGAATTCAGAGAAAGTTTTCAGCTTGACGGTACCTTACCAGAAAGGGATGGCTAAATACGTGCCAGCAGCCGCGGTAATACGTATGTCCCGAGCGTTATCCGGATTTATTGGGCGTAAAGCGAGCGCAGACGGTTTATTAAGTCTGATGTGAAATCCCGAGGCCCAACCTCGGAACTGCATTGGAAACTGATTTACTTGAGTGCGATAGAGGCAAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATGTGGAAGAACACCAGTGGCGAAAGCGGCTTGCTAGATCGTAACTGACGTTGAGGCTCGAAAGTATGGGTAGCAAACGGGATTAGATACCCCGGTAGTCCATACCGTAAACGATGGGTGCTAGTTGTTAAGAGGTTTCCGCCTCCTAGTGACGTAGCAAACGCATTAAGCACCCCGCCTGAGGAGTACGGCCGCAAGGCTAAAACTTAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACGCGAAAAACCTTACCAGGTCTTGACATACCAATGATCGCTTTTGTAATGAAAGCTTTTCTTCGGAACATTGGATACAGGTGGTGCATGGTCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTATTAGTTGCCAGCATTTAGTTGGGCACTCTAATGAGACTGCCGGTGATAAACCGGAGGAAGGTGGGGACGACGTCAGATCATCATGCCCCTTATGACCTGGGCAACACACGTGCTACAATGGGAAGTACAACGAGTCGCAAACCGGCGACGGTAAGCTAATCTCTTAAAACTTCTCTCAGTTCGGACTGGAGTCTGCAACTCGACTCCACGAAGGCGGAATCGCTAGTAATCGCGAATCAGCATGTCGCGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCAAATCATGGGAGTCGGAAGTACCCAAAGTCGCTTGGCTAACTTTTAGAGGCCGGTGCCTAAGGTAAAATCGATGACTGGGATTAAGTCGTAACAAGGTAGCCGTAGGAGAACCTGCGGCTGGATCACCTCCTTTCT")
				.setHspHseq("TTAAATTGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGTGGCGTGCCTAATACATGCAAGTCGTACGCTAGCCGCTGAATTGATCCTTCGGGTGAAGTGAGGCAATGACTAGAGTGGCGAACTGGTGAGTAACACGTAAGAAACCTGCCCTTTAGTGGGGGATAACATTTGGAAACAGATGCTAATACCGCGTAACAACAAATCACACATGTGATCTGTTTGAAAGGTCCTTTTGGATCGCTAGAGGATGGTCTTGCGGCGTATTAGCTTGTTGGTAGGGTAGAAGCCTACCAAGGCAATGATGCGTAGCCGAGTTGAGAGACTGGCCGGCCACATTGGGACTGAGACACTGCCCAAACTCCTACGGGAGGCTGCAGTAGGGAATTTTCCGCAATGCACGAAAGTGTGACGGAGCGACGCCGCGTGTGTGATGAAGGCTTTCGGGTCGTAAAGCACTGTTGTAAGGGAAGAATAACTGAATTCAGAGAAAGTTTTCAGCTTGACGGTACCTTACCAGAAAGGGATGGCTAAATACGTGCCAGCAGCCGCGGTAATACGTATGTCCCGAGCGTTATCCGGATTTATTGGGCGTAAAGCGAGCGCAGACGGTTTATTAAGTCTGATGTGAAATCCCGAGGCCCAACCTCGGAACTGCATTGGAAACTGATTTACTTGAGTGCGATAGAGGCAAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATGTGGAAGAACACCAGTGGCGAAAGCGGCTTGCTAGATCGTAACTGACGTTGAGGCTCGAAAGTATGGGTAGCAAACGGGATTAGATACCCCGGTAGTCCATACCGTAAACGATGGGTGCTAGTTGTTAAGAGGTTTCCGCCTCCTAGTGACGTAGCAAACGCATTAAGCACCCCGCCTGAGGAGTACGGCCGCAAGGCTAAAACTTAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGATACGCGAAAAACCTTACCAGGTCTTGACATACCAATGATCGCTTTTGTAATGAAAGCTTTTCTTCGGAACATTGGATACAGGTGGTGCATGGTCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTATTAGTTGCCAGCATTTAGTTGGGCACTCTAATGAGACTGCCGGTGATAAACCGGAGGAAGGTGGGGACGACGTCAGATCATCATGCCCCTTATGACCTGGGCAACACACGTGCTACAATGGGAAGTACAACGAGTCGCAAACCGGCGACGGTAAGCTAATCTCTTAAAACTTCTCTCAGTTCGGACTGGAGTCTGCAACTCGACTCCACGAAGGCGGAATCGCTAGTAATCGCGAATCAGCATGTCGCGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCAAATCATGGGAGTCGGAAGTACCCAAAGTCGCTTGGCTAACTTTTAGAGGCCGGTGCCTAAGGTAAAATCGATGACTGGGATTAAGTCGTAACAAGGTAGCCGTAGGAGAACCTGCGGCTGGATCACCTCCTTTCT")
				.setHspIdentityString("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
				.createBlastHsp(buildSeq);
		List<Hsp<DNASequence, NucleotideCompound>> hsplist = new ArrayList<>();
		hsplist.add(hsp1hit1res1);
		hsplist.add(hsp1hit1res1);

		BlastHit<DNASequence, NucleotideCompound> hit1res1 = new BlastHitBuilder<DNASequence, NucleotideCompound>()
				.setHitNum(1)
				.setHitId("gnl|BL_ORD_ID|2006")
				.setHitDef("CP000411 Oenococcus oeni PSU-1, complete genome")
				.setHitAccession("0")
				.setHitLen(1780517)
				.setHsps(hsplist)
				.createBlastHit();

		List<Hit<DNASequence, NucleotideCompound>> hitlist = new ArrayList<>();
		hitlist.add(hit1res1);

		BlastResult<DNASequence, NucleotideCompound> res1 = new BlastResultBuilder<DNASequence, NucleotideCompound>()
				.setProgram("blastn")
				.setVersion("BLASTN 2.2.29+")
				.setReference("Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), &quot;A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.")
				.setQueryID("Query_1")
				.setQueryDef("CP000411_-_16S_rRNA Oenococcus oeni PSU-1, complete genome")
				.setQueryLength(1567)
				.createBlastResult();

		Result<DNASequence, NucleotideCompound> expRes1 = results.get(0);
		Hit<DNASequence, NucleotideCompound> expHit1res1 = expRes1.iterator().next();
		Hsp<DNASequence, NucleotideCompound> expHsp1hit1res1 = expHit1res1.iterator().next();

		// result not testable without all hits and hsp
		//assertEquals(expRes1, res1);

		// hit test
		assertEquals(expHit1res1, hit1res1);

		// hsp test
		assertEquals(expHsp1hit1res1, hsp1hit1res1);
	}
	
	private <S extends Sequence<C>, C extends Compound>
	List<Result<S, C>> getResults(File file, Function<String,S> buildSeq) throws IOException, ParseException {
		BlastXMLParser<S,C> instance = new BlastXMLParser<>(buildSeq);
		instance.setFile(file);
		List<Result<S,C>> result = instance.createObjects(1e-10);
		return result;
		
	}
}
