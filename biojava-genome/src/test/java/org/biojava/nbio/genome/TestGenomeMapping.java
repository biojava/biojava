package org.biojava.nbio.genome;

import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import junit.framework.TestCase;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePosition;
import org.biojava.nbio.genome.parsers.genename.GeneChromosomePositionParser;
import org.biojava.nbio.genome.util.ChromosomeMappingTools;
import org.junit.Test;

import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Created by andreas on 7/19/16.
 */
public class TestGenomeMapping extends TestCase {

	private static final String geneChromosomeFile = "http://cdn.rcsb.org/gene/hg38/geneChromosome38.tsf.gz";

	private List<GeneChromosomePosition> gcps = null;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		InputStream input = new GZIPInputStream(new URL(geneChromosomeFile).openStream());
		gcps = GeneChromosomePositionParser.getChromosomeMappings(input);
	}

	@Test
	public void testAK1() {
		String geneName = "AK1";

		assertNotNull(gcps);
		assertTrue("Problems with downloading refFlat file from UCSC browser ", gcps.size() > 100);

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

				assertTrue(pos.getGeneName().equals(geneName));
				assertTrue(pos.getOrientation().equals('-'));
				assertTrue(pos.getChromosome().equals("chr9"));

				List<Range<Integer>> cdsranges = ChromosomeMappingTools.getCDSExonRanges(pos);

				validateExon(0,0,7, cdsranges  );
				validateExon(1,7,43, cdsranges  );
				validateExon(2,43,207, cdsranges  );
				validateExon(3,207,324, cdsranges  );
				validateExon(4,324,516, cdsranges  );
				validateExon(5,516,585, cdsranges  );


				int cdslength = ChromosomeMappingTools.getCDSLength(pos);

				assertTrue("CDS length should be 582, but is " + cdslength, cdslength == (uniProtLength *3));

				List<Range<Integer>> chromranges = ChromosomeMappingTools.getChromosomalRangesForCDS(pos);

				// we are reverse strand. reverse the order
				chromranges = Lists.reverse(chromranges);

				assertTrue(chromranges.size() == 6);

				// compare with https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS6881
				validateExon(0,127868008,127868076, chromranges  );
				validateExon(1,127868320,127868512, chromranges  );
				validateExon(2,127871822,127871939, chromranges  );
				validateExon(3,127872689,127872853, chromranges  );
				validateExon(4,127873025,127873061, chromranges  );
				validateExon(5,127874610,127874617, chromranges  );

			}
		} catch (Exception e) {
			fail(e.getMessage());
		}
	}

	@Test
	public void testHBA(){

		String geneName = "HBA1";
		assertNotNull(gcps);

		assertTrue("Problems with downloading refFlat file from UCSC browser ", gcps.size() > 100);

		try {

			for ( GeneChromosomePosition pos : gcps){

				//System.out.println(pos.getGeneName());
				if ( ! pos.getGeneName().equals(geneName))
					continue;

				assertTrue(pos.getGeneName().equals("HBA1"));
				assertTrue(pos.getGenebankId().equals("NM_000558"));
				assertTrue(pos.getChromosome().equals("chr16"));
				assertTrue(pos.getTranscriptionStart().equals(176650));
				assertTrue(pos.getTranscriptionEnd().equals(177522));
				assertTrue(pos.getOrientation().equals('+'));

				List<Range<Integer>> cdsranges = ChromosomeMappingTools.getCDSExonRanges(pos);

				assertTrue(cdsranges.size() == 3);

				validateExon(0,0,95,cdsranges);
				validateExon(1,95,300,cdsranges);
				validateExon(2,300,429,cdsranges);


				List<Range<Integer>> chromranges = ChromosomeMappingTools.getChromosomalRangesForCDS(pos);

				validateExon(0,176716,176811, chromranges  );
				validateExon(1,176928,177133, chromranges  );
				validateExon(2,177282,177411, chromranges  );


			}
		} catch (Exception e){
			fail(e.getMessage());
		}


	}

	private void validateExon(int exonNr, int start, int stop, List<Range<Integer>> cdsranges) {

		Range<Integer> exon = cdsranges.get(exonNr);
		assertTrue("Exon " + exonNr + " boundary "+ exon.lowerEndpoint()  + " does not match " +start , exon.lowerEndpoint().equals(start));
		assertTrue("Exon " + exonNr + " boundary " + exon.upperEndpoint() + " does not match " + stop, exon.upperEndpoint().equals(stop));

	}

	/** Get the position of the nucleotide base corresponding to the position of that base on the mRNA sequence 
	 * for a gene living on the reverse DNA strand. 
	 * 
	 * @author Yana Valasatava
	 */
	private int getPositionInmRNA(String geneName, String genebankId, int posChrom) {
		for (GeneChromosomePosition gcp : gcps) {	
			if ( gcp.getGeneName().equals(geneName) ) {
				if ( gcp.getGenebankId().equals(genebankId) ) {
					return ChromosomeMappingTools.getCDSPosForChromosomeCoordinate(posChrom, gcp);
				}
			}	
		}
		return -1;
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
					assertEquals(cds, posInmRNA);
					break;
					
				case "last1exon":
			    	posExonStart = 30094003; // starting position of the last exon coding region (on forward strand)
			    	posInmRNA = 51; // position in mRNA sequence equals to the length of the exon
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					assertEquals(cds, posInmRNA);		
					break;
				
				case "last3exon":
					posExonStart = 30103500; // starting position of the first base in a coding region (3rd exon)
					posInmRNA = 257; // position in mRNA sequence equals to the sum length of the 3 last exons 
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					assertEquals(cds, posInmRNA);
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
					assertEquals(cds, posInmRNA);
					break;
					
				case "last1exon":
			    	posExonStart = 99271161; // starting position of the last exon coding region (on forward strand)
			    	posInmRNA = 58; // position in mRNA sequence equals to the length of the exon
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					assertEquals(cds, posInmRNA);		
					break;
				
				case "last3exon":
					posExonStart = 99231345; // starting position of the first base in a coding region (3rd exon)
					posInmRNA = 640; // position in mRNA sequence equals to the sum length of the 3 last exons 
					cds = getPositionInmRNA(geneName, genebankId, posExonStart);
					assertEquals(cds, posInmRNA);
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
		assertEquals(cdsSE, -1);

		int cdsEE = getPositionInmRNA(geneName, genebankId, posExonEnd+1);
		assertEquals(cdsEE, -1);
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
		assertEquals(cdsSE, -1);
		
		int cdsEE = getPositionInmRNA(geneName, genebankId, posExonEnd+1);
		assertEquals(cdsEE, -1);	
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
		assertEquals(cds, mRNAPos);
		
	}
}

