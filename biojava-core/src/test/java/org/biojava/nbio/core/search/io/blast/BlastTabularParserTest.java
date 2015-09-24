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

/**
 *
 * @author Paolo Pavan
 */
public class BlastTabularParserTest {
    
    public BlastTabularParserTest() {
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
     * Test of getFileExtensions method, of class BlastTabularParser.
     */
    @Test
    public void testGetFileExtensions() {
        System.out.println("getFileExtensions");
        BlastTabularParser instance = new BlastTabularParser();
        
        List<String> expResult = new ArrayList();
        expResult.add("blasttabular");
        expResult.add("blasttxt");
        
        List<String> result = instance.getFileExtensions();
        assertEquals(expResult, result);
    }

    /**
     * Test of setFile method, of class BlastTabularParser.
     */
    @Test
    public void testSetFile() {
        System.out.println("setFile");
        File f = null;
        BlastTabularParser instance = new BlastTabularParser();
        instance.setFile(f);
    }

    /**
     * Test of createObjects method, of class BlastTabularParser.
     */
    @Test
    public void testCreateObjects() throws Exception {
        System.out.println("createObjects");
        Result expRes1;
        Hit expHit1res1;
        Hsp expHsp1hit1res1;
        
        String resource = "/org/biojava/nbio/core/search/io/blast/test.two-query.blasttxt";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        BlastTabularParser instance = new BlastTabularParser();
        instance.setFile(file);
        
        List<Result> results = instance.createObjects(1e-10);
        
        BlastHsp hsp1Hit1Res1 = new BlastHspBuilder()
                .setHspNum(1)
                .setPercentageIdentity(97.40/100)
                .setHspAlignLen(77)
                .setMismatchCount(2)
                .setHspGaps(0)
                .setHspQueryFrom(774)
                .setHspQueryTo(850)
                .setHspHitFrom(45396566)
                .setHspQueryTo(45396336)
                .setHspEvalue(1e-46)
                .setHspBitScore(157)
                .createBlastHsp();
        
        BlastHit hit1Res1 = new BlastHitBuilder()
                .setHitDef("chr15")
                .createBlastHit();
        
        BlastResult res1 = new BlastResultBuilder()
                .setQueryID("Query_1")
                .setQueryDef("Dual oxidase (DUOX1_RAT)")
                .createBlastResult();
        
        expRes1 = results.get(0);
        expHit1res1 = expRes1.iterator().next();
        expHsp1hit1res1 = expHit1res1.iterator().next();
        
        // results test
        assertEquals(expRes1, res1);
        // hit test
        assertEquals(expHit1res1, hit1Res1);
        // hsp test
        assertEquals(expHsp1hit1res1, hsp1Hit1Res1);
        
        
        String resource2 = "/org/biojava/nbio/core/search/io/blast/testBlastTabularReport.txt";
        URL resourceURL2 = getClass().getResource(resource2);
        File file2 = new File(resourceURL2.getFile());
        
        BlastTabularParser instance2 = new BlastTabularParser();
        instance2.setFile(file2);
        
        List<Result> results2 = instance2.createObjects(1e-10);
        expRes1 = results2.get(0);
        expHit1res1 = expRes1.iterator().next();
        expHsp1hit1res1 = expHit1res1.iterator().next();
        
        hsp1Hit1Res1 = new BlastHspBuilder()
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
                .createBlastHsp();
        
        hit1Res1 = new BlastHitBuilder()
                .setHitDef("KF031625.1.1775")
                .createBlastHit();
        
        res1 = new BlastResultBuilder()
                .setQueryID("1_759_906_F3")
                .setQueryDef("1_759_906_F3")
                .createBlastResult();
        
        // results test
        assertEquals(expRes1, res1);
        // hit test
        assertEquals(expHit1res1, hit1Res1);
        // hsp test
        assertEquals(expHsp1hit1res1, hsp1Hit1Res1);
    }

    
    
}
