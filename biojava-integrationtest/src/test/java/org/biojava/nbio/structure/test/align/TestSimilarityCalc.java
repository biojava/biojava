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
package org.biojava.nbio.structure.test.align;

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

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;


/**
 * Created by andreas on 8/3/16.
 */
public class TestSimilarityCalc  {

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


            assertTrue(afpChain.getAlnLength() == 71);

            assertTrue(afpChain.getAlnLength() == 71);
            assertTrue(afpChain.getSimilarity() > .57);
            assertTrue(afpChain.getSimilarity() <= .6);

            // this calls the internal invalidate method
            afpChain.setOptAln(afpChain.getOptAln());

            assertTrue("wrong similarity score : " + afpChain.getSimilarity(), afpChain.getSimilarity() > .57);
            assertTrue("wrong similarity score : " + afpChain.getSimilarity(), afpChain.getSimilarity() <= .6);

            assertTrue(afpChain.getSimilarity() > .58);
            assertTrue(afpChain.getSimilarity() < .59);
            assertTrue("similarity score is " + afpChain.getSimilarity()  , afpChain.getSimilarity() > .46);


        } catch (Exception e){

            fail(e.getMessage());
        }

    }
}
