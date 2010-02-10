/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.InputStream;
import java.util.List;

import org.biojava3.core.sequence.ProteinSequence;
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


        FastaReader<ProteinSequence> fastaReader = new FastaReader<ProteinSequence>(inStream, new GenericFastaHeaderParser(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
        List<ProteinSequence> proteinSequences = fastaReader.process();
        inStream.close();


        //Should have 282 sequences
        System.out.println("Expecting 283 got " + proteinSequences.size());
        assertEquals(proteinSequences.size() ,  283 );

        ProteinSequence proteinSequence = proteinSequences.get(0);
        assertEquals(proteinSequence.getAccession().getID(),"A2D504_ATEGE/1-46");
        assertEquals(proteinSequence.getSequenceAsString(),"-----------------FK-N----LP-LED----------------Q----ITL--IQY-----------SWM----------------------CL-SSFA------LSWRSYK---HTNSQFLYFAPDLVF-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

        proteinSequence = proteinSequences.get(281);
        System.out.println(proteinSequence.getAccession());
        System.out.println(proteinSequence.getSequenceAsString());
        assertEquals(proteinSequence.getAccession().getID(),"Q9PU76_CRONI/141-323");
        assertEquals(proteinSequence.getSequenceAsString(),"VETVTELTEFAKSI-PGFS-N----LD-LND----------------Q----VTL--LKY-----------GVY----------------------EA-IFAM------LASVMNK---DGMPVAYGNGFITRE------------------------------------------------------------------------------------------------------------------------------------------------------------FLKSLRKPFCDIMEPKFDFA-MKF-NSL-E-LDDSDI--------------------SLFVA-AIIC-CGDRPG-------------------------------------------LVNV--GHIEKMQESIVHVLKL-H-----LQN---------NH---PD----------------------------DI------F--------LFP-KLLQKMAD-LRQLV-----------------TEH-AQLV--QIIKK---TESDAHLHPLL-------QEI---");

        proteinSequence = proteinSequences.get(282);
        assertEquals(proteinSequence.getAccession().getID(),"Q98SJ1_CHICK/15-61");
        assertEquals(proteinSequence.getSequenceAsString(),"---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------Q-----------------NW------Q--------RFY-QLTKLLDS-MHDVV-----------------ENL-LSFC--FQTFLDKSM--SIEFPEML-------AEI---");

        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }
}
