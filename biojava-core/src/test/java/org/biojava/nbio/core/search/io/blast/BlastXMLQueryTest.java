/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.search.io.blast;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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
 * @author pavanpa
 */
public class BlastXMLQueryTest {
    
    public BlastXMLQueryTest() {
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
     * Test of setFile method, of class BlastXMLQuery.
     */
    @Test
    public void testSetFile() {
        System.out.println("setFile");
        File f = null;
        BlastXMLQuery instance = new BlastXMLQuery();
        instance.setFile(f);
    }

    /**
     * Test of createObjects method, of class BlastXMLQuery.
     */
    @Test
    public void testCreateObjects() throws Exception {
        System.out.println("createObjects");
        
        String resource = "/org/biojava/nbio/core/search/io/blast/testBlastReport.xml";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        BlastXMLQuery instance = new BlastXMLQuery();
        instance.setFile(file);
        
        //instance.setQueryReferences(null);
        //instance.setDatabaseReferences(null);
        ArrayList<Result> result = instance.createObjects(1e-10);
        
        // test with random manual selected results
        BlastHsp hsp1hit1res1 = new BlastHspBuilder()
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
                .createBlastHsp();
        BlastHit hit1res1 = new BlastHitBuilder()
                .setHitNum(1)
                .setHitId("gnl|BL_ORD_ID|2006")
                .setHitDef("gi|265679047|ref|NR_029355.1| Clostridium methylpentosum DSM 5476 strain R2 16S ribosomal RNA, partial sequence")
                .setHitAccession("2006")
                .setHitLen(1435)
                .createBlastHit();
        
        BlastResult res1 = new BlastResultBuilder()
                .setProgram("blastn")
                .setVersion("BLASTN 2.2.29+")
                .setReference("Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), &quot;A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.")
                .setQueryID("Query_1")
                .setQueryDef("42DKN:00022:00047")
                .setQueryLength(226)
                .createBlastResult();
        
        Result expRes1 = result.get(0);
        Hit expHit1res1 = expRes1.iterator().next();
        Hsp expHsp1hit1res1 = expHit1res1.iterator().next();
        
        // result test
        assertEquals(expRes1, res1);
        
        // hit test
        assertEquals(expHit1res1, hit1res1);
        
        // hsp test
        assertEquals(expHsp1hit1res1, hsp1hit1res1);
    }

    /**
     * Test of getFileExtensions method, of class BlastXMLQuery.
     */
    @Test
    @Ignore public void testGetFileExtensions() {
        System.out.println("getFileExtensions");
        BlastXMLQuery instance = new BlastXMLQuery();
        List<String> expResult = null;
        List<String> result = instance.getFileExtensions();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setQueryReferences method, of class BlastXMLQuery.
     */
    @Test
    @Ignore public void testSetQueryReferences() {
        System.out.println("setQueryReferences");
        List sequences = null;
        BlastXMLQuery instance = new BlastXMLQuery();
        instance.setQueryReferences(sequences);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setDatabaseReferences method, of class BlastXMLQuery.
     */
    @Test
    @Ignore public void testSetDatabaseReferences() {
        System.out.println("setDatabaseReferences");
        List sequences = null;
        BlastXMLQuery instance = new BlastXMLQuery();
        instance.setDatabaseReferences(sequences);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of storeObjects method, of class BlastXMLQuery.
     */
    @Test
    public void testStoreObjects() throws Exception {
        // not implemented yet
    }
    
    private File unpackResourceFile(String resource) {
        
        File file = null;
        URL res = getClass().getResource(resource);
        if (res == null) {
            throw new RuntimeException(" not supported");
        }
        if (res.toString().startsWith("jar:")) {
            try {
                InputStream input = getClass().getResourceAsStream(resource);
                file = File.createTempFile("javawrapped", ".exe");
                //if (!DEBUG) file.deleteOnExit();
                file.setExecutable(true, true);
                OutputStream out = new FileOutputStream(file);
                int read;
                byte[] bytes = new byte[1024];

                while ((read = input.read(bytes)) != -1) {
                    out.write(bytes, 0, read);
                }
                out.close();
                out.flush();
            } catch (IOException ex) {
                throw new RuntimeException("Resource not available for this OS!");
            } 
        } else {
            //this will probably work in your IDE, but not from a JAR
            file = new File(res.getFile());
        }

        if (file == null) {
            throw new RuntimeException("Error: File " + file + " not found!");
        }
        return file;
    }
    
}
