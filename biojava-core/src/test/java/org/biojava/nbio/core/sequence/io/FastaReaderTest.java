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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.sequence.io;

import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.logging.Level;

import junit.framework.TestCase;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReaderTest extends TestCase {

	private final static Logger logger = LoggerFactory.getLogger(FastaReaderTest.class);

	public FastaReaderTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@Override
	@Before
	public void setUp() {
	}

	@Override
	@After
	public void tearDown() {
	}

	/**
	 * Test of process method, of class FastaReader.
	 */
	@Test
	public void testProcess() throws Exception {
		logger.info("process");
		InputStream inStream = this.getClass().getResourceAsStream("/PF00104_small.fasta");
		assertNotNull(inStream);

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		LinkedHashMap<String, ProteinSequence> proteinSequences = fastaReader.process();
		inStream.close();

		// Should have 282 sequences
		// logger.debug("Expecting 283 got " + proteinSequences.size());
		assertEquals(proteinSequences.size(), 283);

		int seqNum = 0;
		for (String id : proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch (seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(), "A2D504_ATEGE/1-46");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
				break;
			case 281:
				// logger.debug("Get Accession: {}", proteinSequence.getAccession());
				// logger.debug("Get Protein Sequence: {}", proteinSequence.getSequenceAsString());
				assertEquals(proteinSequence.getAccession().getID(), "Q9PU76_CRONI/141-323");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
				break;
			case 282:
				assertEquals(proteinSequence.getAccession().getID(), "Q98SJ1_CHICK/15-61");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum, 283);
	}

	@Test
	public void testProcessInt() throws Exception {
		logger.info("process(int)");
		InputStream inStream = this.getClass().getResourceAsStream("/PF00104_small.fasta");
		assertNotNull(inStream);
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		LinkedHashMap<String, ProteinSequence> proteinSequences = fastaReader.process(200);

		// Should have 200 sequences
		logger.info("Expecting full chunk size [200], got [{}]", proteinSequences.size());
		assertEquals(proteinSequences.size(), 200);

		int seqNum = 0;
		for (String id : proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch (seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(), "A2D504_ATEGE/1-46");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
				break;
			case 199:
				assertEquals(proteinSequence.getAccession().getID(), "Q5F0P7_HUMAN/248-428");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"DRELVVIIGWAKHI-PGFS-S----LS-LGD----------------Q----MSL--LQS-----------AWM----------------------EI-LILG------IVYRSLP---YDDKLVYAEDYIMD-------------------------------------------------------------------------------------------------------------------------------------------------------------EEHSRLAGLLELYRAILQLV-RRY-KKL-K-VEKEEF--------------------VTLKA-LALA-NSDSMY-------------------------------------------IEDL--EAVQKLQDLLHEALQD-Y-----ELS---------QR---HE----------------------------EP------W--------RTG-KLLLTLPL-LRQTA-----------------AKA-VQHF--YSVKLQGKV--PMH--KLF-------LEM---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum, 200);

		// Should have 83 sequences
		proteinSequences = fastaReader.process(200);
		logger.info("Expecting incomplete chunk size [83], got [{}]", proteinSequences.size());
		assertEquals(proteinSequences.size(), 83);
		seqNum = 0;
		for (String id : proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch (seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(), "RARA_CANFA/233-413");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"TKCIIKTVEFAKQL-PGFT-T----LT-IAD----------------Q----ITL--LKA-----------ACL----------------------DI-LILR------ICTRYTP---EQDTMTFSEGLTLN-------------------------------------------------------------------------------------------------------------------------------------------------------------RTQMHKAGFGPLTDLVFAFA-NQL-LPL-E-MDDAET--------------------GLLSA-ICLI-CGDRQD-------------------------------------------LEQP--DRVDMLQEPLLEALKV-Y-----VRK---------RR---PS----------------------------RP------H--------MFP-KMLMKITD-LRSIS-----------------AKG-AERV--ITLKMEIPG--SMP--PLI-------QEM---");
				break;
			case 81:
				// logger.debug(proteinSequence.getAccession());
				// logger.debug(proteinSequence.getSequenceAsString());
				assertEquals(proteinSequence.getAccession().getID(), "Q9PU76_CRONI/141-323");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
				break;
			case 82:
				assertEquals(proteinSequence.getAccession().getID(), "Q98SJ1_CHICK/15-61");
				assertEquals(
						proteinSequence.getSequenceAsString(),
						"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum, 83);
		fastaReader.close();
		inStream.close();
	}

	@Test
	public void testSmallFastaWithBlankLines() {
		logger.info("testSmallFastaWithBlankLines");

		try {
			InputStream inStream = this.getClass().getResourceAsStream("/test.fasta");

			FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
					new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));

			LinkedHashMap<String, ProteinSequence> b;

			int nrSeq = 0;

			do {
				b = fastaReader.process(10);

				if (nrSeq == 0) {
					logger.info("Expecting full chunk size [10] before end, got [{}]", b.size());
					assertEquals(b.size(), 10);
				} else {
					logger.info("Expecting empty chunk [0] and end, got [{}]", b.size());
					assertEquals(b.size(), 0);
				}

				for (String key : b.keySet()) {
					nrSeq++;

					// #282 would result in an endless loop
					// this makes sure it has been fixed.
					assertTrue("Looks like there is a problem with termination of processing of the FASTA file!",
							nrSeq < 11);
				}

			} while (b.size() == 10);
		} catch (Exception ex) {
			ex.printStackTrace();
			java.util.logging.Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);

			fail(ex.getMessage());
		}
	}

	@Test
	public void testSmallFasta2() {

		try {
			InputStream inStream = this.getClass().getResourceAsStream("/test.fasta");

			FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
					inStream, new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
					new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));

			LinkedHashMap<String, ProteinSequence> b = fastaReader.process();

			assertNotNull(b);

			// #282 make sure that process() still works

			assertTrue(b.keySet().size() == 10);

		} catch (Exception ex) {
			ex.printStackTrace();
			java.util.logging.Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);

			fail(ex.getMessage());
		}
	}
}
