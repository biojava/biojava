/**
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

package org.biojavax.bio.phylo.io.nexus;

//import java.io.BufferedReader;
//import java.io.ByteArrayOutputStream;
import java.io.File;
//import java.io.IOException;
//import java.io.InputStream;
//import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.Set;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.biojava.bio.seq.io.ParseException;
import org.biojavax.bio.phylo.io.nexus.NexusFile;
import org.biojavax.bio.phylo.io.nexus.NexusFileBuilder;
import org.biojavax.bio.phylo.io.nexus.NexusFileFormat;
import org.jgrapht.WeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

/**
 * JUnit test for TreesBlock parsing.
 *
 * Currently very simple, but test for parsing ability.
 *
 * @author Tiago Antao
 */
public class TreesBlockTest extends TestCase {
    NexusFile nexus1;
    NexusFile nexus2;

    public TreesBlockTest(String name) {
        super(name);
    }

    /*
     * Tests a binary tree.
     */
    public void testSimple() {
        doVertexCount(nexus1, "test1", 3);
    }

    /*
     * Tests more than binary.
     */
    public void testThreeOffspring() {
        doVertexCount(nexus1, "test2", 4);
    }

    /*
     * Tests depth >1.
     */
    public void testMoreDepth() {
        doVertexCount(nexus1, "test3", 5);
    }


    /*
     * Tests parsing distance.
     */
    public void testDistanceParsing() {
        doVertexCount(nexus1, "test4", 5);
    }

    /*
     * Tests complex tree, no distance.
     */
    public void testComplex() {
        doVertexCount(nexus1, "test5", 28);
    }

    /*
     * Tests complex tree, with distance.
     */
    public void testComplexWithDistance() {
        doVertexCount(nexus1, "test6", 28);
    }

    /*
     * Tests complex tree with name clashing, no distance.
     */
    public void testClash() {
        doVertexCount(nexus1, "test7", 28);
    }

    /*
     * Tests complex tree, name clashing, with distance.
     */
    public void testClashWithDistance() {
        doVertexCount(nexus1, "test8", 28);
    }

    /*
     * Tests resiliency to comma after last translate info.
     */
    public void testParseIncorrectTranslation() {
        doVertexCount(nexus2, "test1", 3);
    }


    protected void setUp() {
        try {
            NexusFileBuilder builder = new NexusFileBuilder();
            NexusFileFormat.parseInputStream(builder, this.getClass().getResourceAsStream("/test1.nex"));
            nexus1 = builder.getNexusFile();
            NexusFileFormat.parseInputStream(builder, this.getClass().getResourceAsStream("/test2.nex"));
            nexus2 = builder.getNexusFile();
        }
        catch (Exception e) {
            //Should not happen
        }
    }

    private void doVertexCount(NexusFile nexus, String tree, int count) {
        try {
            getTree(nexus, tree);
            assertEquals(countVertexes(), count);
        }
        catch (ParseException pe) {
            fail("ParserFailure: " + pe.getMessage());
        }
    }

    private int countVertexes() {
        Set<String> vertexes = tree.vertexSet();
        int cnt = 0;

        for (String v : vertexes) {
            cnt++;
        }
        return cnt;
    }

    private String topNode;

    private TreesBlock getTreeNode(NexusFile nexus) {
        Iterator it = nexus.blockIterator();
        NexusBlock block;
        while(it.hasNext()) {
            block = (NexusBlock)it.next();
            if (block.getBlockName().equals("TREES")) {
                return (TreesBlock)block;
            }
        }
        return null;
    }

    private WeightedGraph<String, DefaultWeightedEdge> tree;

    private void getTree(NexusFile nexus, String name)
    throws ParseException {
        TreesBlock node = getTreeNode(nexus);
        tree = node.getTreeAsWeightedJGraphT(name);
        //System.out.println(node.getTopNode());
        topNode = node.getTopNode();
    }


    // creates a suite
    public static Test suite()
    {
        TestSuite suite = new TestSuite(TreesBlockTest.class);

        return suite;
    }

    // harness for tests
    public static void main(String [] args)
    {
        junit.textui.TestRunner.run(suite());
    }

}

