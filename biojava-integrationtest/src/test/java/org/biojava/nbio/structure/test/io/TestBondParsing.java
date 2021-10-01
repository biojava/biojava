package org.biojava.nbio.structure.test.io;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TestBondParsing {

    /**
     * Integration test for bond parsing in PDB-format, where author chain ids and asym ids differ and can cause
     * problems. See https://github.com/biojava/biojava/issues/943
     */
    @Test
    public void testIssue943() throws Exception {
        PDBFileReader reader = new PDBFileReader();
        FileParsingParameters params = new FileParsingParameters();
        params.setCreateAtomBonds(true);
        reader.setFileParsingParameters(params);
        Structure s = reader.getStructureById("1v9i");

        Group his95 = s.getPolyChain("A").getAtomGroup(94);
        Atom ne2His95 = his95.getAtom("NE2");
        assertEquals(3, ne2His95.getBonds().size());

        Group zn = s.getNonPolyChain("B").getAtomGroup(0);
        assertEquals(3, zn.getAtom("ZN").getBonds().size());

    }

    /**
     * Integration test for SS bond parsing in PDB-format, where author chain ids and asym ids differ and can cause
     * problems. See https://github.com/biojava/biojava/issues/929
     */
    @Test
    public void testIssue929() throws Exception {
        PDBFileReader reader = new PDBFileReader();
        FileParsingParameters params = new FileParsingParameters();
        params.setCreateAtomBonds(true);
        reader.setFileParsingParameters(params);
        Structure s = reader.getStructureById("1a4w");

        Group cysB = s.getPolyChain("B").getAtomGroup(118);
        Atom sgCysB = cysB.getAtom("SG");
        assertEquals(2, sgCysB.getBonds().size());

        Group cysA = s.getPolyChain("A").getAtomGroup(1);
        Atom sgCysA = cysA.getAtom("SG");
        assertEquals(2, sgCysA.getBonds().size());


    }
}
