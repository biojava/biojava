package org.biojava.nbio.structure.test.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class TestMMCIFWriting {

    /**
     * Test that a bioassembly formed with a symmetry mate gets written out without problems.
     * @throws IOException
     * @throws StructureException
     */
    @Test
    public void testBiounitWriting6a63() throws IOException, StructureException {
        Structure s = StructureIO.getBiologicalAssembly("6a63", 1);
        String mmcif = s.toMMCIF();
        assertNotNull(mmcif);
        assertTrue(mmcif.contains("A_1"));
        assertTrue(mmcif.contains("A_2"));
    }

}
