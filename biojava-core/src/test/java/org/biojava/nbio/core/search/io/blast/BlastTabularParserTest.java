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
        
        List<String> expResult = new ArrayList<String>();
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
        
        String resource = "/org/biojava/nbio/core/search/io/blast/small-blastreport.blasttxt";
        URL resourceURL = getClass().getResource(resource);
        File file = new File(resourceURL.getFile());
        
        BlastTabularParser instance = new BlastTabularParser();
        instance.setFile(file);
        
        List<Result> results = instance.createObjects(1e-10);
        
        BlastHsp hsp1Hit1Res1 = new BlastHspBuilder()
                .setHspNum(1)
                .setPercentageIdentity(100.0/100)
                .setHspAlignLen(1567)
                .setMismatchCount(0)
                .setHspGaps(0)
                .setHspQueryFrom(1)
                .setHspQueryTo(1567)
                .setHspHitFrom(616309)
                .setHspQueryTo(617875)
                .setHspEvalue(0)
                .setHspBitScore(2894)
                .createBlastHsp();
        
        BlastHsp hsp1Hit1Res2 = new BlastHspBuilder()
                .setHspNum(1)
                .setPercentageIdentity(100.0/100)
                .setHspAlignLen(1567)
                .setMismatchCount(0)
                .setHspGaps(0)
                .setHspQueryFrom(1)
                .setHspQueryTo(1567)
                .setHspHitFrom(1278699)
                .setHspQueryTo(1277133)
                .setHspEvalue(0)
                .setHspBitScore(2894)
                .createBlastHsp();
        
        List<Hsp> hsplist = new ArrayList<Hsp>();
        hsplist.add(hsp1Hit1Res1);
        hsplist.add(hsp1Hit1Res2);
        
        BlastHit hit1Res1 = new BlastHitBuilder()
                .setHitDef("CP000411")
                .setHsps(hsplist)
                .createBlastHit();
        List<Hit> hitlist = new ArrayList<Hit>();
        hitlist.add(hit1Res1);
        
        BlastResult res1 = new BlastResultBuilder()
                .setQueryID("CP000411_-_16S_rRNA")
                .setQueryDef("CP000411_-_16S_rRNA")
                .setHits(hitlist)
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
        
        // results test
        assertEquals(expRes1.getQueryID(), "1_759_906_F3");
        assertEquals(results2.size(), 298);
        // only one hsp test
        assertEquals(expHsp1hit1res1, hsp1Hit1Res1);
    }

    
    
}
