package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Site;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

/**
 * Created by Matt on 10/31/2015.
 */
public class TestParseMmCIFFeatures {
    @Test
    public void testSites()throws IOException, StructureException {
        AtomCache cache = new AtomCache();

        StructureIO.setAtomCache(cache);

        cache.setUseMmCif(true);
        Structure sCif = StructureIO.getStructure("4HHB");

        assertNotNull(sCif);

        // After it has read the file, it should check that expected SITES are present.
        List<Site> sites = sCif.getSites();

        // 4HHB has 6 sites from ligands.
        assertEquals(6, sites.size());
    }
}
