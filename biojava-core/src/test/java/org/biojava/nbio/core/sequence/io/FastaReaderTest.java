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
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.logging.Level;

/**
 *
 * @author Scooter Willis 
 */
public class FastaReaderTest {

	private final static Logger logger = LoggerFactory.getLogger(FastaReaderTest.class);

	public FastaReaderTest() {
	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@Before
	public void setUp() {
	}

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
		Assert.assertNotNull(inStream);


		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		Map<String,ProteinSequence> proteinSequences = fastaReader.process();
		inStream.close();

		//Should have 282 sequences
		//logger.debug("Expecting 283 got " + proteinSequences.size());
		Assert.assertEquals(proteinSequences.size(), 283);

		int seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
				case 0:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "A2D504_ATEGE/1-46");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
					break;
				case 281:
					//logger.debug("Get Accession: {}", proteinSequence.getAccession());
					//logger.debug("Get Protein Sequence: {}", proteinSequence.getSequenceAsString());
					Assert.assertEquals(proteinSequence.getAccession().getID(), "Q9PU76_CRONI/141-323");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
					break;
				case 282:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "Q98SJ1_CHICK/15-61");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
					break;
			}
			seqNum++;
		}
		Assert.assertEquals(seqNum, 283);
	}

	@Test
	public void processIntTest() throws Exception {
		logger.info("process(int)");
		InputStream inStream = this.getClass().getResourceAsStream("/PF00104_small.fasta");
		Assert.assertNotNull(inStream);
		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		Map<String,ProteinSequence> proteinSequences = fastaReader.process(200);

		//Should have 200 sequences
		//logger.debug("Expecting 200 got " + proteinSequences.size());
		Assert.assertEquals(proteinSequences.size(), 200);

		int seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
				case 0:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "A2D504_ATEGE/1-46");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
					break;
				case 199:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "Q5F0P7_HUMAN/248-428");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "DRELVVIIGWAKHI-PGFS-S----LS-LGD----------------Q----MSL--LQS-----------AWM----------------------EI-LILG------IVYRSLP---YDDKLVYAEDYIMD-------------------------------------------------------------------------------------------------------------------------------------------------------------EEHSRLAGLLELYRAILQLV-RRY-KKL-K-VEKEEF--------------------VTLKA-LALA-NSDSMY-------------------------------------------IEDL--EAVQKLQDLLHEALQD-Y-----ELS---------QR---HE----------------------------EP------W--------RTG-KLLLTLPL-LRQTA-----------------AKA-VQHF--YSVKLQGKV--PMH--KLF-------LEM---");
					break;
			}
			seqNum++;
		}
		Assert.assertEquals(seqNum, 200);

		//Should have 83 sequences
		proteinSequences = fastaReader.process(200);
		Assert.assertEquals(proteinSequences.size(), 83);
		seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
				case 0:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "RARA_CANFA/233-413");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "TKCIIKTVEFAKQL-PGFT-T----LT-IAD----------------Q----ITL--LKA-----------ACL----------------------DI-LILR------ICTRYTP---EQDTMTFSEGLTLN-------------------------------------------------------------------------------------------------------------------------------------------------------------RTQMHKAGFGPLTDLVFAFA-NQL-LPL-E-MDDAET--------------------GLLSA-ICLI-CGDRQD-------------------------------------------LEQP--DRVDMLQEPLLEALKV-Y-----VRK---------RR---PS----------------------------RP------H--------MFP-KMLMKITD-LRSIS-----------------AKG-AERV--ITLKMEIPG--SMP--PLI-------QEM---");
					break;
				case 81:
					//logger.debug(proteinSequence.getAccession());
					//logger.debug(proteinSequence.getSequenceAsString());
					Assert.assertEquals(proteinSequence.getAccession().getID(), "Q9PU76_CRONI/141-323");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
					break;
				case 82:
					Assert.assertEquals(proteinSequence.getAccession().getID(), "Q98SJ1_CHICK/15-61");
					Assert.assertEquals(proteinSequence.getSequenceAsString(), "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
					break;
			}
			seqNum++;
		}
		Assert.assertEquals(seqNum, 83);
		fastaReader.close();
		inStream.close();
	}

	@Test
	public void testSmallFasta() throws IOException {

		InputStream inStream = this.getClass().getResourceAsStream("/test.fasta");

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));

		Map<String, ProteinSequence> b;

		int nrSeq = 0;

		while ((b = fastaReader.process(10)) != null) {
			for (String key : b.keySet()) {
				nrSeq++;

				// #282 would result in an endless loop
				// this makes sure it has been fixed.
				Assert.assertTrue("Looks like there is a problem with termination of processing of the FASTA file!", nrSeq < 15);
			}

		}

	}


	@Test
	public void testSmallFasta2() throws IOException {

		InputStream inStream = this.getClass().getResourceAsStream("/test.fasta");

		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
				inStream,
				new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));


		int nrSeq = 0;

		Map<String, ProteinSequence> b = fastaReader.process();

		Assert.assertNotNull(b);

		// #282 make sure that process() still works

		Assert.assertEquals(10, b.keySet().size());


	}
}
