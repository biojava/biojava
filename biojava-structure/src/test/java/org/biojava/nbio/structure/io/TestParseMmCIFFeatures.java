package org.biojava.nbio.structure.io;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.SSBond;
import org.biojava.nbio.structure.Site;
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

	public void testSites()throws IOException, StructureException {
        Structure sCif = StructureIO.getStructure("4HHB");

        assertNotNull(sCif);

        // After it has read the file, it should check that expected SITES are present.
        List<Site> sites = sCif.getSites();

        // 4HHB has 6 sites from ligands.
        assertEquals(6, sites.size());

        // Check for each site that it has parsed all residues.
        assertEquals(1, getGroupsInSite(sCif, "AC1"));
        assertEquals(1, getGroupsInSite(sCif, "AC2"));
        assertEquals(16, getGroupsInSite(sCif, "AC3"));
        assertEquals(13, getGroupsInSite(sCif, "AC4"));
        assertEquals(15, getGroupsInSite(sCif, "AC5"));
        assertEquals(7, getGroupsInSite(sCif, "AC6"));

        // Check that struct_site parsing worked, and they have descriptions.
        assertEquals(getDescription(sCif, "AC1"), "BINDING SITE FOR RESIDUE PO4 D 147");
        assertEquals(getDescription(sCif, "AC2"), "BINDING SITE FOR RESIDUE PO4 B 147");
        assertEquals(getDescription(sCif, "AC3"), "BINDING SITE FOR RESIDUE HEM A 142");
        assertEquals(getDescription(sCif, "AC4"), "BINDING SITE FOR RESIDUE HEM B 148");
        assertEquals(getDescription(sCif, "AC5"), "BINDING SITE FOR RESIDUE HEM C 142");
        assertEquals(getDescription(sCif, "AC6"), "BINDING SITE FOR RESIDUE HEM D 148");
    }

    private int getGroupsInSite(Structure structure, String site) {
        for (Site a_site : structure.getSites()) {
            if (a_site.getSiteID().equals(site)) {
                return a_site.getGroups().size();
            }
        }
        return 0;
    }

    private String getDescription(Structure structure, String site) {
        for (Site a_site : structure.getSites()) {
            if (a_site.getSiteID().equals(site)) {
                return a_site.getDescription();
            }
        }
        return "";
    }
}
