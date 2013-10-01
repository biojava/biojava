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
 * Created on Sep 30, 2013
 * Author: blivens
 *
 */

package org.biojava.bio.structure.scop;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

/**
 * Generic tests for ScopDatabases. All implementing classes should pass these tests.
 *
 * This abstract class defines the tests for the interface. Each implementation
 * should extend this class and implement the requirements for a parametic test
 * (annotate the class with @RunWith, add an @Parameters static method,
 * and provide a single constructor calling super)
 * @author blivens
 *
 */
@RunWith(Parameterized.class)
public abstract class ScopDatabaseTest {
    protected ScopDatabase scop;
    String tag;

    /**
     *
     * @param tag A short string, displayed for failed asserts
     * @param scop the database instance to test
     */
    public ScopDatabaseTest(String tag, ScopDatabase scop) {
        if( tag != null) {
            this.tag = "["+tag+"] "; // for messages
        } else {
            this.tag = "";
        }

        this.scop = scop;
    }

    /**
     * Traverse through the SCOP hierarchy
     *
     */
    @Test
    public void traverseHierarchy()
    {
        String pdbId = "4HHB";

        List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);
        assertEquals(tag+"Wrong number of domains",4,domains.size());

        // Check domains (order doesn't matter)
        assertEquals(tag+"Wrong domain","d4hhba_",domains.get(0).getScopId());
        assertEquals(tag+"Wrong domain","d4hhbb_",domains.get(2).getScopId());
        assertEquals(tag+"Wrong domain","d4hhbc_",domains.get(1).getScopId());
        assertEquals(tag+"Wrong domain","d4hhbd_",domains.get(3).getScopId());

        // Check the heirarchy
        ScopNode node = scop.getScopNode(domains.get(0).getSunid());
        ScopDescription desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,15251,node.getSunid());
        assertEquals(tag,"d4hhba_",desc.getName());
        assertEquals(tag,"4hhb A:",desc.getDescription());
        assertEquals(tag,"a.1.1.2",desc.getClassificationId());

        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46487,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"Human (Homo sapiens) [TaxId: 9606]",desc.getDescription());
        assertEquals(tag,"a.1.1.2",desc.getClassificationId());

        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46486,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"Hemoglobin, alpha-chain",desc.getDescription());
        assertEquals(tag,"a.1.1.2",desc.getClassificationId());

        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46463,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"Globins",desc.getDescription());
        assertEquals(tag,"a.1.1.2",desc.getClassificationId());

        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46458,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"Globin-like",desc.getDescription());
        assertEquals(tag,"a.1.1",desc.getClassificationId());


        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46457,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"Globin-like",desc.getDescription());
        assertEquals(tag,"a.1",desc.getClassificationId());

        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,46456,node.getSunid());
        assertEquals(tag,"-",desc.getName());
        assertEquals(tag,"All alpha proteins",desc.getDescription());
        assertEquals(tag,"a",desc.getClassificationId());

        // root node
        node = scop.getScopNode(node.getParentSunid());
        desc = scop.getScopDescriptionBySunid(node.getSunid());
        assertEquals(tag,0,node.getSunid());
        assertNull(tag+"Root should not have a description", desc);
    }

    /** Get various categories
     *
     */
    public void getCategories(){
        List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

        assertNotNull(tag,superfams);
        if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
            assertEquals(tag,2223,superfams.size());
        } else {
            // defaults for other versions
            assertFalse(tag,superfams.isEmpty());
        }

        List<ScopDescription> folds = scop.getByCategory(ScopCategory.Fold);

        assertNotNull(tag,folds);
        assertFalse(tag,folds.isEmpty());
        if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
            assertEquals(tag,1393,folds.size());
        }

    }

    public void alignSuperfamily(){
        // download SCOP if required and load into memory
        ScopDatabase scop = ScopFactory.getSCOP();
        List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

        System.out.println("Total nr. of superfamilies:" + superfams.size());


        // configure where to load PDB files from and
        // what information to load
        AtomCache cache = new AtomCache();
        FileParsingParameters fileparams = new FileParsingParameters() ;
        fileparams.setAlignSeqRes(false);
        fileparams.setLoadChemCompInfo(true);
        fileparams.setParseSecStruc(false);
        cache.setFileParsingParams(fileparams);

        // get tge first superfamily
        ScopDescription superfam1 = superfams.get(0);
        System.out.println("First superfamily: " + superfam1);

        ScopNode node = scop.getScopNode(superfam1.getSunID());
        System.out.println("scopNode for first superfamily:" + node);

        List<ScopDomain> doms4superfam1 = scop.getScopDomainsBySunid(superfam1.getSunID());
        ScopDomain dom1 = doms4superfam1.get(0);

        // align the first domain against all others members of this superfamily
        for ( int i = 1 ; i < doms4superfam1.size() ; i ++){

            ScopDomain dom2 = doms4superfam1.get(i);

            try {
                Structure s1 = cache.getStructureForDomain(dom1);
                Structure s2 = cache.getStructureForDomain(dom2);

                Atom[] ca1 = StructureTools.getAtomCAArray(s1);
                Atom[] ca2 = StructureTools.getAtomCAArray(s2);
                StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
                AFPChain afpChain = ce.align(ca1, ca2);

                //System.out.println(afpChain.toCE(ca1, ca2));

                //StructureAlignmentDisplay.display(afpChain, ca1, ca2);

                System.out.println(dom1.getScopId() + " vs. " + dom2.getScopId()+ " :" + afpChain.getProbability());
                double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
                afpChain.setTMScore(tmScore);
                System.out.println(AfpChainWriter.toScoresList(afpChain));

            } catch (Exception e){
                e.printStackTrace();
            }
        }

    }

    public void printDomainsForPDB(){
        String pdbId = "4HHB";

        // download SCOP if required and load into memory
        ScopDatabase scop = ScopFactory.getSCOP();

        List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);

        System.out.println(domains);

    }





    @Test
    public void testComments() {
        List<String> comments;

        // root node
        comments = scop.getComments(0);
        assertTrue(comments.isEmpty());

        //TODO add additional version checks, since comments change a lot

        if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) <= 0 ) {
            // Note: only tested so far with 1.75, so may need some modification for earlier versions

            comments = scop.getComments(127355);
            assertEquals(tag+"Wrong number of comments", 2, comments.size());
            assertEquals(tag+"Wrong comment", "automatically matched to d2hbia_", comments.get(0).trim());
            assertEquals(tag+"Wrong comment", "complexed with hem; mutant", comments.get(1).trim());
        }
        if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75) == 0 ) {
            comments = scop.getComments(160555);
            assertEquals(tag+"Wrong number of comments", 1, comments.size());
            assertEquals(tag+"Wrong comment", "<a href=\"http://pfam.sanger.ac.uk/family?acc=PF06262\">PF06262</a>; DUF1025; minimal zincin fold that retains 3-stranded mixed beta-sheet and contains HExxH motif in the C-terminal helix; no metal ion bound to this motif is observed in the first determined structures", comments.get(0));


        }
        if(scop.getScopVersion().compareToIgnoreCase( ScopFactory.VERSION_1_75B) >= 0 ) {
            // The following were added or modified in 1.75B

            comments = scop.getComments(127355);
            assertEquals(tag+"Wrong number of comments", 2, comments.size());
            assertEquals(tag+"Wrong comment", "automatically matched to d2hbia_", comments.get(0).trim());
            assertEquals(tag+"Wrong comment", "complexed with hem; mutant", comments.get(1).trim());

            comments = scop.getComments(160555);
            assertEquals(tag+"Wrong number of comments", 1, comments.size());
            assertEquals(tag+"Wrong comment", "PF06262; DUF1025; minimal zincin fold that retains 3-stranded mixed beta-sheet and contains HExxH motif in the C-terminal helix; no metal ion bound to this motif is observed in the first determined structures", comments.get(0).trim());

            // d3ueea_ was added in 1.75B update
            // domain
            comments = scop.getComments(190700);
            assertEquals(tag+"Wrong number of comments", 1, comments.size());
            assertEquals(tag+"Wrong comment", "not a true protein", comments.get(0));

            // fold
            comments = scop.getComments(57923);
            assertEquals(tag+"Wrong number of comments", 1, comments.size());
            assertEquals(tag+"Wrong comment", "metal(zinc)-bound alpha+beta fold", comments.get(0));
        }
    }

}
