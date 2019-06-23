package org.biojava.nbio.structure.test.io.cif;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.cif.CifFileConverter;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class CifFileSupplierIntegrationTest {
    @Test
    public void test1SMT() throws IOException {
        // an x-ray structure
        testRoundTrip("1SMT");
    }

    /**
     * MMCIF write test for an NMR structure with 2 chains
     * @throws IOException propagated
     */
    @Test
    public void test2N3J() throws IOException {
        // an NMR structure (multimodel) with 2 chains
        testRoundTrip("2N3J");
    }

    @Test
    public void test1A2C() throws IOException {
        // a structure with insertion codes
        testRoundTrip("1A2C");
    }

    private static void testRoundTrip(String pdbId) throws IOException {
        URL url = new URL("https://files.rcsb.org/download/" + pdbId + ".cif");
        Structure originalStruct = CifFileConverter.fromURL(url);

        InputStream inputStream = new ByteArrayInputStream(CifFileConverter.toText(originalStruct).getBytes());
        Structure readStruct = CifFileConverter.fromInputStream(inputStream);

        assertNotNull(readStruct);
        assertEquals(originalStruct.getChains().size(), readStruct.getChains().size());
        assertEquals(originalStruct.nrModels(), readStruct.nrModels());

        for (int i = 0; i < originalStruct.nrModels(); i++) {
            assertEquals(originalStruct.getModel(i).size(), readStruct.getModel(i).size());
        }

        for (int modelIdx = 0; modelIdx < originalStruct.nrModels(); modelIdx++) {
            for (int i = 0; i < originalStruct.getModel(modelIdx).size(); i++) {
                assertEquals(originalStruct.getChains().get(i).getAtomGroups().size(),
                        readStruct.getChains().get(i).getAtomGroups().size());

                Chain origChain = originalStruct.getModel(modelIdx).get(i);
                Chain readChain = readStruct.getModel(modelIdx).get(i);

                assertEquals(origChain.getAtomGroups().size(), readChain.getAtomGroups().size());
                assertEquals(origChain.getId(), readChain.getId());
                assertEquals(origChain.getName(), readChain.getName());

                Atom[] origAtoms = StructureTools.getAllAtomArray(origChain);
                Atom[] readAtoms = StructureTools.getAllAtomArray(readChain);

                assertEquals(origAtoms.length, readAtoms.length);

                for (int atomIdx = 0; atomIdx < origAtoms.length; atomIdx++) {
                    assertEquals("atom serials don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getPDBserial(), readAtoms[atomIdx].getPDBserial());

                    assertEquals("atom names don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getName(), readAtoms[atomIdx].getName());

                    assertEquals("atom elements don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getElement(), readAtoms[atomIdx].getElement());

                    assertEquals("x values don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getX(), readAtoms[atomIdx].getX(),0.0001);

                    assertEquals("y values don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getY(), readAtoms[atomIdx].getY(),0.0001);

                    assertEquals("z values don't match for atom " + origAtoms[atomIdx].toString(),
                            origAtoms[atomIdx].getZ(), readAtoms[atomIdx].getZ(),0.0001);
                }
            }
        }

        // Test cell and symmetry
        assertEquals(originalStruct.getCrystallographicInfo().getSpaceGroup(),
                readStruct.getCrystallographicInfo().getSpaceGroup());
    }

    /**
     * Tests that structures containing symmetry mates with modified chain identifiers
     * can be written out correctly.
     */
    @Test
    public void testBiounitWriting() throws IOException {
        Structure s = createDummyStructure();
        String mmcif = CifFileConverter.toText(s);
        String[] lines = mmcif.split("\n");
        long atomLines = Arrays.stream(lines).filter(l -> l.startsWith("ATOM")).count();
        assertNotNull(mmcif);
        assertEquals(4, atomLines);
    }

    private static Structure createDummyStructure() {
        Group g = new AminoAcidImpl();
        Atom a = getAtom("CA", Element.C, 1, 1, 1, 1);
        g.addAtom(a);
        g.setResidueNumber(new ResidueNumber("A", 1, null));
        Group altLocG = new AminoAcidImpl();
        Atom a2 = getAtom("CA", Element.C, 2, 2, 2, 2);
        altLocG.addAtom(a2);
        altLocG.setResidueNumber(new ResidueNumber("A", 1, null));

        g.addAltLoc(altLocG);

        Chain c1 = new ChainImpl();
        c1.addGroup(g);
        c1.setId("A");
        EntityInfo entityInfo = new EntityInfo();
        entityInfo.setMolId(1);
        entityInfo.addChain(c1);
        c1.setEntityInfo(entityInfo);

        Group gc2 = new AminoAcidImpl();
        Atom ac2 = getAtom("CA", Element.C, 3, 3, 3, 3);
        gc2.addAtom(ac2);
        gc2.setResidueNumber(new ResidueNumber("A_1", 1, null));

        Group altLocGc2 = new AminoAcidImpl();
        Atom ac22 = getAtom("CA", Element.C, 4, 4, 4, 4);
        altLocGc2.addAtom(ac22);
        altLocGc2.setResidueNumber(new ResidueNumber("A_1", 1, null));

        gc2.addAltLoc(altLocGc2);

        Chain c2 = new ChainImpl();
        c2.addGroup(gc2);
        c2.setId("A_1");
        c2.setEntityInfo(entityInfo);
        entityInfo.addChain(c2);

        Structure s = new StructureImpl();
        s.addChain(c1);
        s.addChain(c2);
        return s;
    }

    private static Atom getAtom(String name, Element e, int id, double x, double y, double z) {
        Atom a = new AtomImpl();
        a.setX(x);
        a.setY(y);
        a.setZ(z);
        a.setPDBserial(id);
        a.setName(name);
        a.setElement(e);
        return a;
    }
}