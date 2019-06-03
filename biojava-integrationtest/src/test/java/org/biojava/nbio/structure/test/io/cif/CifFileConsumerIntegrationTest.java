package org.biojava.nbio.structure.test.io.cif;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.DBRef;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.BcifFileReader;
import org.biojava.nbio.structure.io.CifFileReader;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import static org.junit.Assert.*;

public class CifFileConsumerIntegrationTest {
    private static boolean headerOnly;
    private static boolean binary;

    @Test
    public void testLoad() throws IOException {
        headerOnly = false;
        doTestLoad();
    }

    @Test
    public void testLoadHeaderOnly() throws IOException {
        headerOnly = true;
        doTestLoad();
    }

    @Test
    public void testLoadBinary() throws IOException {
        headerOnly = false;
        binary = true;
        doTestLoad();
    }

    @Test
    public void testLoadHeaderOnlyBinary() throws IOException {
        headerOnly = true;
        binary = true;
        doTestLoad();
    }

    private void doTestLoad() throws IOException {
        // test a simple protein
        comparePDB2cif("5pti", "A");

        // test a protein with modified residues
        comparePDB2cif("1a4w", "L");
        comparePDB2cif("1a4w", "H");
        comparePDB2cif("1a4w", "I");

        //non-standard encoded amino acid
        comparePDB2cif("1fdo", "A");

        // test a DNA binding protein
        comparePDB2cif("1j59", "A");
        comparePDB2cif("1j59", "E");

        // test a NMR protein
        comparePDB2cif("2kc9", "A");
    }

    private void comparePDB2cif(String id, String chainId) throws IOException {
        String fileName = binary ? "/" + id + ".bcif" : "/" + id + ".cif";
        System.out.println(fileName);
        InputStream inStream = getClass().getResourceAsStream(fileName);
        assertNotNull("Could not find file " + fileName + ". Config problem?", inStream);

        LocalPDBDirectory reader = binary ? new BcifFileReader() : new CifFileReader();

        FileParsingParameters params = new FileParsingParameters();
        params.setHeaderOnly(headerOnly);
        reader.setFileParsingParameters(params);

        Structure cifStructure = reader.getStructure(inStream);
        assertNotNull(cifStructure);

        // load the PDB file via the PDB parser
        Structure pdbStructure;
        InputStream pinStream = this.getClass().getResourceAsStream("/" + id + ".pdb");
        assertNotNull(inStream);

        PDBFileParser pdbParser = new PDBFileParser();
        pdbParser.setFileParsingParameters(params);

        pdbStructure = pdbParser.parsePDBFile(pinStream);

        assertNotNull(pdbStructure);

        // check NMR data
        assertEquals(id + ": the isNMR flag is not the same!",
                pdbStructure.isNmr(),
                cifStructure.isNmr());

        if (pdbStructure.isNmr()) {
            assertEquals(id + ": the nr of NMR models is not the same!",
                    pdbStructure.nrModels(),
                    pdbStructure.nrModels());
            checkNMR(pdbStructure);
            checkNMR(cifStructure);
        }

        Chain a_pdb = pdbStructure.getPolyChainByPDB(chainId);
        Chain a_cif = cifStructure.getPolyChainByPDB(chainId);

        String pdb_SEQseq = a_pdb.getSeqResSequence();
        String cif_SEQseq = a_cif.getSeqResSequence();

        assertEquals(id + ": the SEQRES sequences don't match!",
                pdb_SEQseq,
                cif_SEQseq);

        assertEquals(id + ":  The nr of ATOM groups does not match!",
                a_pdb.getAtomGroups(GroupType.AMINOACID).size(),
                a_cif.getAtomGroups(GroupType.AMINOACID).size());

        // actually this check not necessarily works, since there can be waters in PDB that we don;t deal with yet in cif...
        for (int i = 0; i < a_pdb.getAtomGroups(GroupType.AMINOACID).size(); i++) {
            Group gp = a_pdb.getAtomGroups(GroupType.AMINOACID).get(i);
            List<Group> cifGroups = a_cif.getAtomGroups(GroupType.AMINOACID);
            Group gc = cifGroups.get(i);
            checkGroups(gp, gc);
        }

        String pdb_seq = a_pdb.getAtomSequence();
        String cif_seq = a_cif.getAtomSequence();

        assertEquals("the sequences obtained from PDB and mmCif don't match!", pdb_seq, cif_seq);

        List<DBRef> pdb_dbrefs = pdbStructure.getDBRefs();
        List<DBRef> cif_dbrefs = cifStructure.getDBRefs();

        assertEquals("nr of DBrefs found does not match!", pdb_dbrefs.size(), cif_dbrefs.size());

        DBRef p = pdb_dbrefs.get(0);
        DBRef c = cif_dbrefs.get(0);

        String pdb_dbref = p.toPDB();
        String cif_dbref = c.toPDB();
        assertEquals("DBRef is not equal", pdb_dbref, cif_dbref);

        PDBHeader h1 = pdbStructure.getPDBHeader();
        PDBHeader h2 = cifStructure.getPDBHeader();

        if (!h1.toPDB().toUpperCase().equals(h2.toPDB().toUpperCase())) {
            System.err.println(h1.toPDB());
            System.err.println(h2.toPDB());
            assertEquals(h1.toPDB(), h2.toPDB());
        }
        assertEquals("the PDBHeader.toPDB representation is not equivalent",
                h1.toPDB().toUpperCase(),
                h2.toPDB().toUpperCase());
    }

    private void checkGroups(Group g1, Group g2) {
        String pdbId1 = g1.getChain().getStructure().getPDBCode();
        String pdbId2 = g1.getChain().getStructure().getPDBCode();
        assertEquals(pdbId1, pdbId2);

        assertEquals(g1.getType(), g2.getType());
        assertEquals(g1.getResidueNumber().getSeqNum(), g2.getResidueNumber().getSeqNum());
        assertEquals(g1.getResidueNumber().getInsCode(), g2.getResidueNumber().getInsCode());
        assertEquals(g1.getPDBName(), g2.getPDBName());
        assertEquals(g1.has3D(), g2.has3D());

        assertEquals(g1.hasAltLoc(), g2.hasAltLoc());
        assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2 + ":" + g2, g1.getAltLocs().size(), g2.getAltLocs().size());
        assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2 + ":" + g2, g1.getAtoms().size(), g2.getAtoms().size());

        if (g1.has3D()) {
            Atom a1 = g1.getAtom(0);
            Atom a2 = g2.getAtom(0);
            if (a1 == null) {
                fail("could not get atom for group " + g1);
            }
            if (a2 == null) {
                fail("could not get atom for group " + g2);
            }
            assertEquals(a1.getX(), a2.getX(), 0.0001);
            assertEquals(a1.getOccupancy(), a2.getOccupancy(), 0.0001);
            assertEquals(a1.getTempFactor(), a2.getTempFactor(), 0.0001);
            assertEquals(a1.getName(), a2.getName());
        }
    }

    private void checkNMR(Structure s) {
        assertTrue(s.isNmr());
        int models = s.nrModels();
        assertTrue(models > 0);
        List<Chain> model0 = s.getModel(0);

        // compare with all others
        for (int i = 1; i < models; i++) {
            List<Chain> modelX = s.getModel(i);
            assertEquals(model0.size(), modelX.size());

            // compare lengths:
            for (int j = 0; j < model0.size(); j++) {
                Chain c1 = model0.get(j);
                Chain cx = modelX.get(j);
                assertEquals(c1.getAtomLength(), cx.getAtomLength());
                assertEquals(c1.getAtomSequence(), cx.getAtomSequence());
                assertEquals(c1.getAtomGroups(GroupType.AMINOACID).size(), cx.getAtomGroups(GroupType.AMINOACID).size());
                assertEquals(c1.getAtomGroups(GroupType.NUCLEOTIDE).size(), cx.getAtomGroups(GroupType.NUCLEOTIDE).size());
                assertEquals(c1.getAtomGroups(GroupType.HETATM).size(), cx.getAtomGroups(GroupType.HETATM).size());
            }
        }
    }
}