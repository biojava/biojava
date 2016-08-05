package org.biojava.nbio.structure.gui;

import junit.framework.TestCase;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.seq.SmithWaterman3DParameters;
import org.biojava.nbio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Created by andreas on 8/3/16.
 */
public class TestSimilarityCalc extends TestCase{

    @Test
    public void testSimilarityDisplay(){

        String name1 = "1CDG.A";
        String name2 = "1TIM.A";

        AtomCache cache = new AtomCache();

        Structure structure1 = null;
        Structure structure2 = null;

        try {

            StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName);

            SmithWaterman3DParameters params = new SmithWaterman3DParameters();
            
            structure1 = cache.getStructure(name1);
            structure2 = cache.getStructure(name2);

            Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
            Atom[] ca2 = StructureTools.getAtomCAArray(structure2);


            AFPChain afpChain = algorithm.align(ca1, ca2, params);

            afpChain.setName1(name1);
            afpChain.setName2(name2);



            // get the scores
            int ca1Length = afpChain.getCa1Length();
            int ca2Length = afpChain.getCa2Length();

            int blockNum = afpChain.getBlockNum();

            int optLength = afpChain.getOptLength();
            double totalRmsdOpt = afpChain.getTotalRmsdOpt();

            double alignScore = afpChain.getAlignScore();
            int alnLength = afpChain.getAlnLength();

            assertTrue(afpChain.getAlnLength() == 71);
            int gapLen = afpChain.getGapLen();


            assertTrue(afpChain.getAlnLength() == 71);
            assertTrue(afpChain.getSimilarity() > .57);
            assertTrue(afpChain.getSimilarity() <= .6);

            // this calls the internal invalidate method
            afpChain.setOptAln(afpChain.getOptAln());

            assertTrue("wrong similarity score : " + afpChain.getSimilarity(), afpChain.getSimilarity() > .57);
            assertTrue("wrong similarity score : " + afpChain.getSimilarity(), afpChain.getSimilarity() <= .6);


            double similarity = afpChain.getSimilarity();
            double identity = afpChain.getIdentity();

            assertTrue(afpChain.getSimilarity() > .58);
            assertTrue(afpChain.getSimilarity() < .59);
            assertTrue("similarity score is " + afpChain.getSimilarity()  , afpChain.getSimilarity() > .46);







        } catch (Exception e){

            fail(e.getMessage());
        }

    }
}
