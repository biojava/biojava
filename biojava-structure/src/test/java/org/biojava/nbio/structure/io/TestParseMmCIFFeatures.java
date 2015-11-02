package org.biojava.nbio.structure.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.SSBond;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * Created by larsonmattr on 10/31/2015.
 */
public class TestParseMmCIFFeatures {
    @Test
    public void testSSBond()throws IOException, StructureException {    	
        AtomCache cache = new AtomCache();

        StructureIO.setAtomCache(cache);

        cache.setUseMmCif(true);
        Structure sCif = StructureIO.getStructure("2OYA");

        assertNotNull(sCif);

        // After it has read the file, it should check that expected SSBONDs are present.
        List<SSBond> bonds = sCif.getSSBonds();

        // 2OYA has 6 ssbonds, 3 on molecule A and 3 on molecule B
        assertEquals(6, bonds.size());

        // Check the bonds
        assertEquals("A446-A507", printBond(bonds.get(0)));
        assertEquals("A459-A517", printBond(bonds.get(1)));
        assertEquals("A487-A497", printBond(bonds.get(2)));
        assertEquals("B446-B507", printBond(bonds.get(3)));
        assertEquals("B459-B517", printBond(bonds.get(4)));
        assertEquals("B487-B497", printBond(bonds.get(5)));
    }
    
    String printBond(SSBond bond) {
    	StringBuilder str = new StringBuilder();
    	
    	str.append(bond.getChainID1());
    	str.append(bond.getResnum1());
    	str.append("-");
    	str.append(bond.getChainID2());
    	str.append(bond.getResnum2());
    	
    	return str.toString();
    }
}