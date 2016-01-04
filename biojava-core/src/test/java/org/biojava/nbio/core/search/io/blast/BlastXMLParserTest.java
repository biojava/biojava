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

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.core.search.io.Hit;
import org.biojava.nbio.core.search.io.Hsp;
import org.biojava.nbio.core.search.io.Result;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Paolo Pavan
 */
public class BlastXMLParserTest {
    
    public BlastXMLParserTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of setFile method, of class BlastXMLParser.
     */
    @Test
    public void testSetFile() {
        System.out.println("setFile");
        File f = null;
        BlastXMLParser instance = new BlastXMLParser();
        instance.setFile(f);
    }

    /**
     * Test of createObjects method, of class BlastXMLParser.
     */
    @Test
    public void testCreateObjects() throws Exception {
        System.out.println("createObjects");
        
        String resource = "/org/biojava/nbio/core/search/io/blast/small-blastreport.blastxml";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        BlastXMLParser instance = new BlastXMLParser();
        instance.setFile(file);
        
        //instance.setQueryReferences(null);
        //instance.setDatabaseReferences(null);
        List<Result> result = instance.createObjects(1e-10);
        
        // test with random manual selected results
        BlastHsp hsp1hit1res1 = new BlastHspBuilder()
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
                .createBlastHsp();
        List<Hsp> hsplist = new ArrayList<Hsp>();
        hsplist.add(hsp1hit1res1);
        hsplist.add(hsp1hit1res1);
        
        BlastHit hit1res1 = new BlastHitBuilder()
                .setHitNum(1)
                .setHitId("gnl|BL_ORD_ID|2006")
                .setHitDef("CP000411 Oenococcus oeni PSU-1, complete genome")
                .setHitAccession("0")
                .setHitLen(1780517)
                .setHsps(hsplist)
                .createBlastHit();
        
        List<Hit> hitlist = new ArrayList<Hit>();
        hitlist.add(hit1res1);
        
        BlastResult res1 = new BlastResultBuilder()
                .setProgram("blastn")
                .setVersion("BLASTN 2.2.29+")
                .setReference("Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), &quot;A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.")
                .setQueryID("Query_1")
                .setQueryDef("CP000411_-_16S_rRNA Oenococcus oeni PSU-1, complete genome")
                .setQueryLength(1567)
                .createBlastResult();
        
        Result expRes1 = result.get(0);
        Hit expHit1res1 = expRes1.iterator().next();
        Hsp expHsp1hit1res1 = expHit1res1.iterator().next();
        
        // result not testable without all hits and hsp
        //assertEquals(expRes1, res1);
        
        // hit test
        assertEquals(expHit1res1, hit1res1);
        
        // hsp test
        assertEquals(expHsp1hit1res1, hsp1hit1res1);
    }

    /**
     * Test of getFileExtensions method, of class BlastXMLParser.
     */
    @Test
    public void testGetFileExtensions() {
        System.out.println("getFileExtensions");
        BlastXMLParser instance = new BlastXMLParser();
        List<String> result = instance.getFileExtensions();
        assertTrue(result.contains("blastxml"));
    }

    /**
     * Test of setQueryReferences method, of class BlastXMLParser.
     */
    @Test
    @Ignore public void testSetQueryReferences() {
        System.out.println("setQueryReferences");
        List<org.biojava.nbio.core.sequence.template.Sequence> sequences = null;
        BlastXMLParser instance = new BlastXMLParser();
        instance.setQueryReferences(sequences);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setDatabaseReferences method, of class BlastXMLParser.
     */
    @Test
    @Ignore public void testSetDatabaseReferences() {
        System.out.println("setDatabaseReferences");
        List<org.biojava.nbio.core.sequence.template.Sequence> sequences = null;
        BlastXMLParser instance = new BlastXMLParser();
        instance.setDatabaseReferences(sequences);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of storeObjects method, of class BlastXMLParser.
     */
    @Test
    public void testStoreObjects() throws Exception {
        // not implemented yet
    }
}
