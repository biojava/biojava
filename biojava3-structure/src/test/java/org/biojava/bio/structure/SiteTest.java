/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.bio.structure;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import junit.framework.TestCase;

import org.biojava.bio.structure.io.PDBParseException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;


/**
 * Tests functionality of Site class.
 * @author Jules Jacobsen <jacobsen@ebi.ac.uk>
 */
public class SiteTest extends TestCase {

    public static final String newline = System.getProperty("line.separator");

    private static Site catSite;
    private static Site bindingSite;

    private static List<Group> bindingSiteGroups;

    public SiteTest() {
        try {
            // sites CAT and AC1 from PDB entry 1a4w
            //SITE     1 CAT  3 HIS H  57  ASP H 102  SER H 195
            //SITE     1 AC1  6 ARG H 221A LYS H 224  HOH H 403  HOH H 460
            //SITE     2 AC1  6 HOH H 464  HOH H 497
            //groups for site CAT
            Chain chain = new ChainImpl();
            chain.setChainID("H");
            Group his57 = new AminoAcidImpl();
//            his57.setPDBCode("57");
            his57.setResidueNumber("H", 57, ' ');
            his57.setPDBName("HIS");
            his57.setChain(chain);
            Group asp102 = new AminoAcidImpl();
            asp102.setResidueNumber("H", 102, ' ');
            asp102.setPDBName("ASP");
            asp102.setChain(chain);
            Group ser195 = new AminoAcidImpl();
            ser195.setResidueNumber("H", 195, ' ');
            ser195.setPDBName("SER");
            ser195.setChain(chain);
            List<Group> catSiteGroups = new ArrayList<Group>();
            catSiteGroups.add(his57);
            catSiteGroups.add(asp102);
            catSiteGroups.add(ser195);
            //make the catalytic site CAT
            catSite = new Site();
            catSite.setSiteID("CAT");
            catSite.setGroups(catSiteGroups);
            catSite.setEvCode("UNKNOWN");
            catSite.setDescription("ACTIVE SITE");
            //groups for site AC1
            Group arg221a = new AminoAcidImpl();
            arg221a.setResidueNumber("H", 221, 'A');
            arg221a.setPDBName("ARG");
            arg221a.setChain(chain);
            Group lys224 = new AminoAcidImpl();
            lys224.setResidueNumber("H", 224, ' ');
            lys224.setPDBName("LYS");
            lys224.setChain(chain);
            Group hoh403 = new HetatomImpl();
            hoh403.setResidueNumber("H", 403, ' ');
            hoh403.setPDBName("HOH");
            hoh403.setChain(chain);
            Group hoh460 = new HetatomImpl();
            hoh460.setResidueNumber("H", 460, ' ');
            hoh460.setPDBName("HOH");
            hoh460.setChain(chain);
            Group hoh464 = new HetatomImpl();
            hoh464.setResidueNumber("H", 464, ' ');
            hoh464.setPDBName("HOH");
            hoh464.setChain(chain);
            Group hoh497 = new HetatomImpl();
            hoh497.setResidueNumber("H", 497, ' ');
            hoh497.setPDBName("HOH");
            hoh497.setChain(chain);
            
            bindingSiteGroups = new ArrayList<Group>();
            
            bindingSiteGroups.add(arg221a);
            bindingSiteGroups.add(lys224);
            bindingSiteGroups.add(hoh403);
            bindingSiteGroups.add(hoh460);
            bindingSiteGroups.add(hoh464);
            bindingSiteGroups.add(hoh497);
            //make the binding site AC1
            bindingSite = new Site();
            bindingSite.setSiteID("AC1");
            bindingSite.setGroups(bindingSiteGroups);
            bindingSite.setEvCode("SOFTWARE");
            bindingSite.setDescription("BINDING SITE FOR RESIDUE NA H 541");

        } catch (PDBParseException ex) {
            Logger.getLogger(SiteTest.class.getName()).log(Level.SEVERE, null, ex);
        }
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
     * Test of toPDB method, of class Site.
     */
    @Test
    public void testToPDB_0args() {
//        System.out.println("toPDB_0args");
        String expResult =  "SITE     1 AC1  6 ARG H 221A LYS H 224  HOH H 403  HOH H 460                    " + newline +
                            "SITE     2 AC1  6 HOH H 464  HOH H 497                                          "+ newline;
        String result = bindingSite.toPDB();
//        System.out.println("Expected:");
//        System.out.println(expResult);
//        System.out.println("Got:");
//        System.out.println(result);
        assertEquals(expResult, result);
    }

    /**
     * Test of toPDB method, of class Site.
     */
    @Test
    public void testToPDB_StringBuffer() {
//        System.out.println("toPDB");
        StringBuffer buf = new StringBuffer("");
        String expResult =  "SITE     1 AC1  6 ARG H 221A LYS H 224  HOH H 403  HOH H 460                    " + newline +
                            "SITE     2 AC1  6 HOH H 464  HOH H 497                                          "+ newline;
        bindingSite.toPDB(buf);
        String result = buf.toString();
//        System.out.println("Expected:");
//        System.out.println(expResult);
//        System.out.println("Got:");
//        System.out.println(result);
        assertEquals(expResult, result);
    }


        /**
     * Test of toPDB method, of class Site.
     */
    @Test
    public void testRemark800ToPDB_0args() {
//        System.out.println("remark800toPDB_0args");
        String expResult =  "REMARK 800 SITE_IDENTIFIER: CAT                                                 " + newline +
                            "REMARK 800 EVIDENCE_CODE: UNKNOWN                                               " + newline +
                            "REMARK 800 SITE_DESCRIPTION: ACTIVE SITE                                        " + newline;
        String result = catSite.remark800toPDB();
//        System.out.println(result);
        assertEquals(expResult, result);
    }

    /**
     * Test of toPDB method, of class Site.
     */
    @Test
    public void testRemark800ToPDB_StringBuffer() {
//        System.out.println("remark800toPDB");
        StringBuffer buf = new StringBuffer("");
        String expResult =  "REMARK 800 SITE_IDENTIFIER: CAT                                                 " + newline +
                            "REMARK 800 EVIDENCE_CODE: UNKNOWN                                               " + newline +
                            "REMARK 800 SITE_DESCRIPTION: ACTIVE SITE                                        " + newline;
        catSite.remark800toPDB(buf);
        String result = buf.toString();
//        System.out.println(result);
        assertEquals(expResult, result);
    }



    /**
     * Test of getSiteID method, of class Site.
     */
    @Test
    public void testGetSiteID() {
//        System.out.println("getSiteID");
        String expResult = "CAT";
        String result = catSite.getSiteID();
        assertEquals(expResult, result);
    }

    /**
     * Test of getGroups method, of class Site.
     */
    @Test
    public void testGetGroups() {
//        System.out.println("getGroups");
        List<Group> expResult = bindingSiteGroups;
        List<Group> result = bindingSite.getGroups();
        assertEquals(expResult, result);
    }

    /**
     * Test to see how the groups have been set in the Groups list
     */
    @Test
    public void testGroup() {
        List<Group> result = bindingSite.getGroups();
        Group arg221 = result.get(0);
        ResidueNumber testResNum = new ResidueNumber("H", 221, 'A');
//        testResNum.setChainId("H");
//        testResNum.setSeqNum(221);
//        testResNum.setInsCode("A");
//        System.out.println(arg221);
        assertEquals(testResNum, arg221.getResidueNumber());
        //test the chainId is also set
        assertEquals("H", arg221.getChainId());


        Group hoh403 = null;

        for (Group group : result) {
            if (group.getResidueNumber().getSeqNum() == 403) {
                hoh403 = group;
            }
        }

        ResidueNumber testResNum2 = new ResidueNumber("H", 403, ' ');
//        testResNum2.setChainId("H");
//        testResNum2.setSeqNum(403);
//        testResNum2.setInsCode("");
//        System.out.println(hoh403);
        assertEquals(testResNum2, hoh403.getResidueNumber());
        //test the chaiId is also set
        assertEquals("H", hoh403.getChainId());
    }
}