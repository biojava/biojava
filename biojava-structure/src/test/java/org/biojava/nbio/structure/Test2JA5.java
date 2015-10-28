package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;

import org.junit.Test;

import static org.junit.Assert.*;

import java.io.IOException;

/**
 * Created by andreas on 9/16/15.
 */
public class Test2JA5 {

    @Test
    public void test2JA5() throws IOException, StructureException {

        FileParsingParameters fileParsingParameters = new FileParsingParameters();
        fileParsingParameters.setStoreEmptySeqRes(true);
        fileParsingParameters.setLoadChemCompInfo(true);
        fileParsingParameters.setHeaderOnly(true);

        AtomCache cache = new AtomCache();
        cache.setUseMmCif(false);
        cache.setFileParsingParams(fileParsingParameters);

        StructureIO.setAtomCache(cache);


        Structure s1 = StructureIO.getStructure("2ja5");

        assertTrue(StructureTools.getNrAtoms(s1) == 0);

        assertTrue(s1.getChains().size() == 14);

        Chain nChain = null;
        try {
        	nChain = s1.getChainByPDB("N");
        } catch (StructureException e){
        	// this is expected here, since there is no chain N
        }
        assertNull(nChain);



        
    }
}
