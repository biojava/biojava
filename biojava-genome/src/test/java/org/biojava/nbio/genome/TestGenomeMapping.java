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
package org.biojava.nbio.genome;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import org.biojava.nbio.core.util.Download;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePosition;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePositionParser;
import org.biojava.nbio.genome.util.ChromosomeMappingTools;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by andreas on 7/19/16.
 */
public class TestGenomeMapping {

	private static final String geneChromosomeFile = "http://cdn.rcsb.org/gene/hg38/geneChromosome38.tsf.gz";

	private static final List<GeneChromosomePosition> gcps;

	static {
		List<GeneChromosomePosition> g;
		try {
			g = GeneChromosomePositionParser.getChromosomeMappings( Download.stream(new URL(geneChromosomeFile)) );
		} catch (IOException e) {
			e.printStackTrace();
			g = null;
		}
		gcps = g;
	}

//	@Before
//	public void setUp() throws Exception {
//
//	}

	@Test
	public void testAK1() {
		String geneName = "AK1";

		Assert.assertNotNull(gcps);
		Assert.assertTrue("Problems with downloading refFlat file from UCSC browser ", gcps.size() > 100);

		int uniProtLength = 194;

		try {

			for (GeneChromosomePosition pos : gcps) {

				//System.out.println(pos.getGeneName());
				if (!pos.getGeneName().equals(geneName))
					continue;

				/// there are three alternative transcripts for AK1.
				// we are just testing one here:

				if ( ! pos.getGenebankId().equals("NM_000476"))
					continue;

				Assert.assertEquals(pos.getGeneName(), geneName);
				Assert.assertEquals('-', (char) pos.getOrientation());
				Assert.assertEquals("chr9", pos.getChromosome());

				List<Range<Integer>> cdsranges = ChromosomeMappingTools.getCDSExonRanges(pos);

				validateExon(0,0,7, cdsranges  );
				validateExon(1,7,43, cdsranges  );
				validateExon(2,43,207, cdsranges  );
				validateExon(3,207,324, cdsranges  );
				validateExon(4,324,516, cdsranges  );
				validateExon(5,516,585, cdsranges  );


				int cdslength = ChromosomeMappingTools.getCDSLength(pos, 1);

				Assert.assertEquals("CDS length should be 582, but is " + cdslength, cdslength, (uniProtLength * 3));

				List<Range<Integer>> chromranges = ChromosomeMappingTools.getChromosomalRangesForCDS(pos);

				// we are reverse strand. reverse the order
				chromranges = Lists.reverse(chromranges);

				Assert.assertEquals(6, chromranges.size());

				// compare with https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS6881
				validateExon(0,127868008,127868076, chromranges  );
				validateExon(1,127868320,127868512, chromranges  );
				validateExon(2,127871822,127871939, chromranges  );
				validateExon(3,127872689,127872853, chromranges  );
				validateExon(4,127873025,127873061, chromranges  );
				validateExon(5,127874610,127874617, chromranges  );

			}
		} catch (Exception e) {
			Assert.fail(e.getMessage());
		}
	}

	@Test
	public void testHBA(){

		String geneName = "HBA1";
		Assert.assertNotNull(gcps);

		Assert.assertTrue("Problems with downloading refFlat file from UCSC browser ", gcps.size() > 100);

		try {

			for ( GeneChromosomePosition pos : gcps){

				//System.out.println(pos.getGeneName());
				if ( ! pos.getGeneName().equals(geneName))
					continue;

				Assert.assertEquals("HBA1", pos.getGeneName());
				Assert.assertEquals("NM_000558", pos.getGenebankId());
				Assert.assertEquals("chr16", pos.getChromosome());
				Assert.assertEquals(176650, (int) pos.getTranscriptionStart());
				Assert.assertEquals(177522, (int) pos.getTranscriptionEnd());
				Assert.assertEquals('+', (char) pos.getOrientation());

				List<Range<Integer>> cdsranges = ChromosomeMappingTools.getCDSExonRanges(pos);

				Assert.assertEquals(3, cdsranges.size());

				validateExon(0,0,95,cdsranges);
				validateExon(1,95,300,cdsranges);
				validateExon(2,300,429,cdsranges);


				List<Range<Integer>> chromranges = ChromosomeMappingTools.getChromosomalRangesForCDS(pos);

				validateExon(0,176716,176811, chromranges  );
				validateExon(1,176928,177133, chromranges  );
				validateExon(2,177282,177411, chromranges  );


			}
		} catch (Exception e){
			Assert.fail(e.getMessage());
		}


	}

	private void validateExon(int exonNr, int start, int stop, List<Range<Integer>> cdsranges) {

		Range<Integer> exon = cdsranges.get(exonNr);
		Assert.assertEquals("Exon " + exonNr + " boundary " + exon.lowerEndpoint() + " does not match " + start, (int) exon.lowerEndpoint(), start);
		Assert.assertEquals("Exon " + exonNr + " boundary " + exon.upperEndpoint() + " does not match " + stop, (int) exon.upperEndpoint(), stop);

	}

	/** Get the position of the nucleotide base corresponding to the position of that base on the mRNA sequence
	 * for a gene living on the reverse DNA strand.
	 *
	 * @author Yana Valasatava
	 */
	private int getPositionInmRNA(String geneName, String genebankId, int posChrom, int base) {
		for (GeneChromosomePosition gcp : gcps) {
			if ( gcp.getGeneName().equals(geneName) ) {
				if ( gcp.getGenebankId().equals(genebankId) ) {
					return ChromosomeMappingTools.getCDSPosForChromosomeCoordinate(posChrom, gcp, base);
				}
			}
		}
		return -1;
	}
	private int getPositionInmRNA(String geneName, String genebankId, int posChrom) {
		return getPositionInmRNA(geneName, genebankId, posChrom, 1);
	}

	/** Make sure the mapping tool correctly retrieves the mRNA position for a gene
	 * living on the forward DNA strand for different chromosome positions.
	 *
	 * @author Yana Valasatava
	 */
	@Test
	public void testForwardMappingPositions() {

		String geneName = "HORMAD2"; // gene on the forward DNA strand
		String genebankId = "NM_152510"; // GeneBank ID for the transcript used for testing (ENST00000336726)

		List<String> scenarios = Arrays.asList("first1exon", "last1exon", "last3exon");

		int cds;
		int posExonStart;
		int posInmRNA;
		for (String scenario : scenarios) {

			switch (scenario) {

				case "first1exon":
			    	posExonStart = 30093953; // ending position of the last exon coding region (on forward strand)
			    	posInmRNA = 1; // base 1 position in mRNA sequence
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;

				case "last1exon":
			    	posExonStart = 30094003; // starting position of the last exon coding region (on forward strand)
			    	posInmRNA = 51; // position in mRNA sequence equals to the length of the exon
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;

				case "last3exon":
					posExonStart = 30103500; // starting position of the first base in a coding region (3rd exon)
					posInmRNA = 257; // position in mRNA sequence equals to the sum length of the 3 last exons
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;
			}
		}
	}

	/** Make sure the mapping tool correctly retrieves the mRNA position for a gene
	 * living on the reverse DNA strand for different chromosome positions.
	 *
	 * @author Yana Valasatava
	 */
	@Test
	public void testReverseMappingPositions() {

		String geneName = "BCL11B"; // gene on the reverse DNA strand
		String genebankId = "NM_138576"; // GeneBank ID for the transcript used for testing (ENST00000357195)

		List<String> scenarios = Arrays.asList("first1exon", "last1exon", "last3exon");

		int cds;
		int posExonStart;
		int posInmRNA;
		for (String scenario : scenarios) {

			switch (scenario) {

				case "first1exon":
			    	posExonStart = 99271218; // ending position of the last exon coding region (on forward strand)
			    	posInmRNA = 1; // base 1 position in mRNA sequence
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;

				case "last1exon":
			    	posExonStart = 99271161; // starting position of the last exon coding region (on forward strand)
			    	posInmRNA = 58; // position in mRNA sequence equals to the length of the exon
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;

				case "last3exon":
					posExonStart = 99231345; // starting position of the first base in a coding region (3rd exon)
					posInmRNA = 640; // position in mRNA sequence equals to the sum length of the 3 last exons
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					Assert.assertEquals(cds, posInmRNA);
					break;
			}
		}
	}

	/** Test to make sure the mapping tool correctly identify that position falls outside the coding region
	 * for a gene living on the forward DNA strand.
	 *
	 * @author Yana Valasatava
	 */
	@Test
	public void testForwardMappingForExonBoundaries() {

		String geneName = "HBA1"; // gene on the reverse DNA strand
		String genebankId = "NM_000558"; // GeneBank ID for the transcript used for testing (ENST00000320868)

		int posExonStart = 176717; // starting position of the first base in a coding region (1st exon)
		int posExonEnd = 176811; // ending position of the first base in a coding region (1st exon)

		int cdsSE = getPositionInmRNA(geneName, genebankId, posExonStart-1);
		Assert.assertEquals(cdsSE, -1);

		int cdsEE = getPositionInmRNA(geneName, genebankId, posExonEnd+1);
		Assert.assertEquals(cdsEE, -1);
	}

	/** Test to make sure the mapping tool correctly identify that position falls outside the coding region
	 * for a gene living on the reverse DNA strand.
	 *
	 * @author Yana Valasatava
	 */
	@Test
	public void testReverseMappingForExonBoundaries() {

		String geneName = "BCL11B"; // gene on the reverse DNA strand
		String genebankId = "NM_138576"; // GeneBank ID for the transcript used for testing (ENST00000357195)

		int posExonStart = 99174151; // starting position of the first base in a coding region (1st exon)
		int posExonEnd = 99176195; // ending position of the first base in a coding region (1st exon)

		int cdsSE = getPositionInmRNA(geneName, genebankId, posExonStart-1);
		Assert.assertEquals(cdsSE, -1);

		int cdsEE = getPositionInmRNA(geneName, genebankId, posExonEnd+1);
		Assert.assertEquals(cdsEE, -1);
	}

	/** Test to make sure the mapping tool correctly converts the genetic position to a position on mRNA
	 * when multiple UTR regions are consecutive.
	 *
	 * @author Yana Valasatava
	 */
	@Test
	public void testMappingCromosomePosTomRNAMultiUTRs() {

		String geneName = "ILK"; // gene on the reverse DNA strand
		String genebankId = "NM_001278442"; // GeneBank ID for the transcript used for testing (ENST00000532063)

		int chromPos = 6608760;
		int mRNAPos = 16;

		int cds = getPositionInmRNA(geneName, genebankId, chromPos);
		Assert.assertEquals(cds, mRNAPos);

	}

	@Test
	public void testGenomeMappingToolGetCDSRanges(){

		List<Integer> lst1 = new ArrayList<>(Arrays.asList(86346823, 86352858, 86354529));
		List<Integer> lst2 = new ArrayList<>(Arrays.asList(86348878, 86352984, 86354692));

		int cdsStart=86348749, cdsEnd=86387027;

		List<Range<Integer>> result = ChromosomeMappingTools.getCDSRegions(lst1,lst2,cdsStart,cdsEnd);

		// makes sure the first list does not get  changed;
		Assert.assertEquals(86346823, (int) lst1.get(0));


		Assert.assertEquals(86348749, (int) result.get(0).lowerEndpoint());
		Assert.assertEquals(86352858, (int) result.get(1).lowerEndpoint());
		Assert.assertEquals(86354529, (int) result.get(2).lowerEndpoint());

		Assert.assertEquals(86348878, (int) result.get(0).upperEndpoint());
		Assert.assertEquals(86352984, (int) result.get(1).upperEndpoint());
		Assert.assertEquals(86387027, (int) result.get(2).upperEndpoint());

	}

	@Test
	public void testGenomeMappingToolGetCDSRangesSERINC2(){

		List<Integer> lst1 = new ArrayList<>(Arrays.asList(31413812, 31415872, 31423692));
		List<Integer> lst2 = new ArrayList<>(Arrays.asList(31414777, 31415907, 31423854));

		int cdsStart=31423818, cdsEnd=31434199;

		List<Range<Integer>> result = ChromosomeMappingTools.getCDSRegions(lst1,lst2,cdsStart,cdsEnd);

		// makes sure the first list does not get  changed;
		Assert.assertEquals(31423818, (int) result.get(0).lowerEndpoint());

	}
}

