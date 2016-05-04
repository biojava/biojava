package org.biojava.nbio.structure.io.mmcif;

import junit.framework.TestCase;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Created by andreas on 5/3/16.
 */
public class TestParseInternalChainId extends TestCase{

    public void test2I13(){


        AtomCache cache = new AtomCache();
        cache.setUseMmCif(true);

        try {
            Structure s = cache.getStructure("2I13");

            System.out.println(s);


            assertTrue(s.getPolyChains().size() == 6);
            assertTrue(s.getNonPolyChains().size() >20);

            /** the four nucleic chains
             *
             */
            Chain asymA = s.getPolyChain("A");
            Chain asymB = s.getPolyChain("B");
            Chain asymC = s.getPolyChain("C");
            Chain asymD = s.getPolyChain("D");

            /** the protein chains
             *
             */
            Chain asymE = s.getPolyChain("E");
            Chain asymF = s.getPolyChain("F");

            Chain[] nucleicChains = new Chain[]{asymA,asymB,asymC,asymD};
            Chain[] proteinChains = new Chain[]{asymE,asymF};

            for ( Chain c : proteinChains){
                System.out.println(c);
                assertNotNull("Chain is null!",c);
            }

            for ( Chain c : nucleicChains){
                System.out.println(c);
                assertNotNull("Chain is null!", c);
            }


            assertTrue(asymA.getChainID().equals("A"));
            assertTrue(asymA.getName().equals("C"));

            assertTrue(asymB.getChainID().equals("B"));
            assertTrue(asymB.getName().equals("D"));

            assertTrue(asymC.getChainID().equals("C"));
            assertTrue(asymC.getName().equals("E"));

            assertTrue(asymD.getChainID().equals("D"));
            assertTrue(asymD.getName().equals("F"));

            assertTrue(asymE.getChainID().equals("E"));
            assertTrue(asymE.getName().equals("A"));

            assertTrue(asymF.getChainID().equals("F"));
            assertTrue(asymF.getName().equals("B"));

            Chain chainG = s.getNonPolyChain("G");
            assertTrue(chainG.getName().equals("A"));

        } catch (Exception e){
            e.printStackTrace();
            fail(e.getMessage());
        }
    }


}
