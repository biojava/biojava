/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.InputStream;
import java.util.LinkedHashMap;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReaderTest {

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
		System.out.println("process");
		InputStream inStream = this.getClass().getResourceAsStream("/PF00104_small.fasta");
		assertNotNull(inStream);


		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		LinkedHashMap<String,ProteinSequence> proteinSequences = fastaReader.process();
		inStream.close();

		//Should have 282 sequences
		//System.out.println("Expecting 283 got " + proteinSequences.size());
		assertEquals(proteinSequences.size() ,  283 );

		int seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(),"A2D504_ATEGE/1-46");
				assertEquals(proteinSequence.getSequenceAsString(),"-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
				break;
			case 281:
				//System.out.println(proteinSequence.getAccession());
				//System.out.println(proteinSequence.getSequenceAsString());
				assertEquals(proteinSequence.getAccession().getID(),"Q9PU76_CRONI/141-323");
				assertEquals(proteinSequence.getSequenceAsString(),"VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
				break;
			case 282:
				assertEquals(proteinSequence.getAccession().getID(),"Q98SJ1_CHICK/15-61");
				assertEquals(proteinSequence.getSequenceAsString(),"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum,283);
	}
	
	@Test
	public void processIntTest() throws Exception {
		System.out.println("process(int)");
		InputStream inStream = this.getClass().getResourceAsStream("/PF00104_small.fasta");
		assertNotNull(inStream);
		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(inStream, new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		LinkedHashMap<String,ProteinSequence> proteinSequences = fastaReader.process(200);

		//Should have 200 sequences
		//System.out.println("Expecting 200 got " + proteinSequences.size());
		assertEquals(proteinSequences.size() ,  200 );
		
		int seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(),"A2D504_ATEGE/1-46");
				assertEquals(proteinSequence.getSequenceAsString(),"-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
				break;
			case 199:
				assertEquals(proteinSequence.getAccession().getID(),"Q5F0P7_HUMAN/248-428");
				assertEquals(proteinSequence.getSequenceAsString(),"DRELVVIIGWAKHI-PGFS-S----LS-LGD----------------Q----MSL--LQS-----------AWM----------------------EI-LILG------IVYRSLP---YDDKLVYAEDYIMD-------------------------------------------------------------------------------------------------------------------------------------------------------------EEHSRLAGLLELYRAILQLV-RRY-KKL-K-VEKEEF--------------------VTLKA-LALA-NSDSMY-------------------------------------------IEDL--EAVQKLQDLLHEALQD-Y-----ELS---------QR---HE----------------------------EP------W--------RTG-KLLLTLPL-LRQTA-----------------AKA-VQHF--YSVKLQGKV--PMH--KLF-------LEM---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum,200);
		
		//Should have 83 sequences
		proteinSequences = fastaReader.process(200);
		assertEquals(proteinSequences.size() , 83 );
		seqNum = 0;
		for(String id:proteinSequences.keySet()) {
			ProteinSequence proteinSequence = proteinSequences.get(id);
			switch(seqNum) {
			case 0:
				assertEquals(proteinSequence.getAccession().getID(),"RARA_CANFA/233-413");
				assertEquals(proteinSequence.getSequenceAsString(), "TKCIIKTVEFAKQL-PGFT-T----LT-IAD----------------Q----ITL--LKA-----------ACL----------------------DI-LILR------ICTRYTP---EQDTMTFSEGLTLN-------------------------------------------------------------------------------------------------------------------------------------------------------------RTQMHKAGFGPLTDLVFAFA-NQL-LPL-E-MDDAET--------------------GLLSA-ICLI-CGDRQD-------------------------------------------LEQP--DRVDMLQEPLLEALKV-Y-----VRK---------RR---PS----------------------------RP------H--------MFP-KMLMKITD-LRSIS-----------------AKG-AERV--ITLKMEIPG--SMP--PLI-------QEM---");
				break;
			case 81:
				//System.out.println(proteinSequence.getAccession());
				//System.out.println(proteinSequence.getSequenceAsString());
				assertEquals(proteinSequence.getAccession().getID(),"Q9PU76_CRONI/141-323");
				assertEquals(proteinSequence.getSequenceAsString(),"VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");
				break;
			case 82:
				assertEquals(proteinSequence.getAccession().getID(),"Q98SJ1_CHICK/15-61");
				assertEquals(proteinSequence.getSequenceAsString(),"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");
				break;
			}
			seqNum++;
		}
		assertEquals(seqNum,83);
		fastaReader.close();
		inStream.close();
	}
}
