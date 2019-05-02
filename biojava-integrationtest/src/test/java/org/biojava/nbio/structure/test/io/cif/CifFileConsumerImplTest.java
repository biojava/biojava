package org.biojava.nbio.structure.test.io.cif;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.*;
import org.biojava.nbio.structure.io.cif.CifFileConverter;
import org.junit.Ignore;
import org.junit.Test;
import org.rcsb.cif.CifReader;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

public class CifFileConsumerImplTest {
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
        comparePDB2cif("5pti","A");

        // test a protein with modified residues
        comparePDB2cif("1a4w","L");
        comparePDB2cif("1a4w","H");
        comparePDB2cif("1a4w","I");

        //non-standard encoded amino acid
        comparePDB2cif("1fdo","A");

        // test a DNA binding protein
        comparePDB2cif("1j59","A");
        comparePDB2cif("1j59","E");

        // test a NMR protein
        comparePDB2cif("2kc9","A");
    }

    private void comparePDB2cif(String id, String chainId) throws IOException {
        String fileName = binary ? "/" + id + ".bcif" : "/" + id + ".cif";
        System.out.println(fileName);
        InputStream inStream = getClass().getResourceAsStream(fileName);
        assertNotNull("Could not find file " + fileName + ". Config problem?" , inStream);

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

        if (pdbStructure.isNmr()){
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
        for (int i = 0 ; i < a_pdb.getAtomGroups(GroupType.AMINOACID).size(); i++){
            Group gp = a_pdb.getAtomGroups(GroupType.AMINOACID).get(i);
            List<Group> cifGroups = a_cif.getAtomGroups(GroupType.AMINOACID);
            Group gc = cifGroups.get(i);
            checkGroups(gp, gc);
        }

        String pdb_seq = a_pdb.getAtomSequence();
        String cif_seq = a_cif.getAtomSequence();

        assertEquals("the sequences obtained from PDB and mmCif don't match!", pdb_seq, cif_seq);

        List<DBRef> pdb_dbrefs= pdbStructure.getDBRefs();
        List<DBRef> cif_dbrefs= cifStructure.getDBRefs();

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

    private void checkGroups(Group g1, Group g2){
        String pdbId1 = g1.getChain().getStructure().getPDBCode();
        String pdbId2 = g1.getChain().getStructure().getPDBCode();
        assertEquals(pdbId1, pdbId2);

        assertEquals(g1.getType(), g2.getType());
        assertEquals(g1.getResidueNumber().getSeqNum(), g2.getResidueNumber().getSeqNum());
        assertEquals(g1.getResidueNumber().getInsCode(), g2.getResidueNumber().getInsCode());
        assertEquals(g1.getPDBName(), g2.getPDBName());
        assertEquals(g1.has3D(), g2.has3D());

        assertEquals(g1.hasAltLoc(), g2.hasAltLoc());
        assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2 + ":"+ g2, g1.getAltLocs().size(), g2.getAltLocs().size());
        assertEquals(pdbId1 + ":" + g1 + " - " + pdbId2 + ":"+ g2, g1.getAtoms().size(), g2.getAtoms().size());

        if (g1.has3D()) {
            Atom a1 = g1.getAtom(0);
            Atom a2 = g2.getAtom(0);
            if (a1 == null) {
                fail("could not get atom for group " + g1);
            }
            if (a2 == null) {
                fail("could not get atom for group " + g2);
            }
            assertEquals(a1.getX(),a2.getX(), 0.0001);
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
        for (int i = 1 ; i < models; i++){
            List<Chain> modelX = s.getModel(i);
            assertEquals(model0.size(), modelX.size());

            // compare lengths:
            for (int j = 0 ; j < model0.size(); j++){
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

    /**
     * There were issues when parsing files in parallel using SimpleDateFormat. The culprit are not the actual files,
     * but rather the parallel execution. No problems spotted by this test, run for whole archive is now flawless after
     * dropping SimpleDateFormat for date parsing.
     */
    @Test
    @Ignore("ignored for now as Bcif file source may change - currently using local files")
    public void testFailingEntries() {
        Pattern.compile(", ").splitAsStream("4he8, 1z4u, 4fp1, 1blc, 4cit, 2y2y, 4exq, 2n0f, 2d9o, 2v16, 1kqv, " +
                "1bwo, 2k2g, 1qhd, 5mhj, 2dn3, 5pq8, 5cay, 6ms1, 2vhu, 2gi0, 3swe, 3daz, 5yel, 2pxp, 4uis, 3cs1, 3in5," +
                " 1sl3, 4hjc, 3hj2, 5kpi, 1gyq, 1yq8, 4yqz, 1ox3, 2pls, 1vne, 4q02, 1dtt, 1jau, 5h3b, 5sxk, 4el7, 5q7w," +
                " 4zuz, 1n6i, 1dhg, 3dhe, 2gpo, 5if7, 5ld8, 1jhz, 4fr3, 1r6u, 3hdl, 5fse, 1iho, 1t10, 2oc6, 3czx, 3b3o," +
                " 5i6w, 2ecv, 4l2x, 441d, 2i0x, 1xq4, 3tbb, 4mmz, 1qew, 6i16, 1t8d, 5w7r, 6gm1, 1s7u, 2qp3, 1cf3, 4myb," +
                " 1omh, 1zog, 2b68, 1nqb, 1t7k")
                .parallel()
                .forEach(pdbId -> {
                    System.out.println(pdbId);
                    Structure cif = loadLocalCif(pdbId);
                    assertNumberFormat(cif);
                    Structure bcif = loadLocalBcif(pdbId);
                    assertNumberFormat(bcif);
                });
    }

    private Structure loadLocalCif(String pdbId) {
        try {
            String middle = pdbId.substring(1, 3);
            return CifFileConverter.convert(CifReader.readText(new GZIPInputStream(Files.newInputStream(Paths.get("/var/pdb/" + middle + "/" + pdbId + ".cif.gz")))));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private Structure loadLocalBcif(String pdbId) {
        try {
            String middle = pdbId.substring(1, 3);
            return CifFileConverter.convert(CifReader.readBinary(Files.newInputStream(Paths.get("/var/bcif/" + middle + "/" + pdbId + ".bcif"))));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private void assertNumberFormat(Structure structure) {
        PDBHeader header = structure.getPDBHeader();
        assertNotNull(header);
        Date relDate = header.getRelDate();
        assertNotNull(relDate);
        Date modDate = header.getModDate();
        assertNotNull(modDate);
    }

    enum Source {
        BCIF, // binary cif using ciftools
        CIF, // mmCIF using ciftools
        MMCIF, // mmCIF using BioJava impl
        MMTF // mmtf using BioJava impl
    }

    /**
     * Performance diary;
     *
     * 05/01/19 - ciftools v0.3.0, parallel, bcif, non-gzipped, 12 worker threads
     * BCIF: 918 s for 151079 structures, 6073 µs per structure, failed for 0 entries
     * CIF: 2502 s for 151079 structures, 16562 µs per structure, failed for 0 entries
     *
     * 05/02/19 -  ciftools v0.3.0, parallel, bcif, non-gzipped for bcif / gzipped for cif, 12 worker threads
     * BCIF: 35 s for 5000 structures, 6228 µs per structure, failed for 0 entries
     * CIF: 132 s for 5000 structures, 23902 µs per structure, failed for 0 entries
     * MMCIF: 310 s for 5000 structures, 52329 µs per structure, failed for 0 entries
     * MMTF: 17 s for 5000 fetched structures, 3211 µs per structure, failed for 0 entries
     *
     * CIF parsing using David's tokenizer approach takes roughly half the time. Binary parsing using ciftools takes
     * roughly double the time compared to MMTF and should be even more at a disadvantage when MMTF data is stored
     * locally. However, MMTF provides almost no metadata.
     *
     * @throws IOException propagated
     */
    @Test
    @Ignore("ignore long-running test, do run to track performance")
    public void parseEntireArchive() throws IOException {
        // TODO create an objective test case - all local, all/none gzipped, mmtf parses almost no annotation data
        AtomicInteger counter = new AtomicInteger(0);
        long start = System.nanoTime();
        int chunkSize = 250;
        List<String> failed = Collections.synchronizedList(new ArrayList<>());

        Source source = Source.CIF;
        Path archivePath = null;
        switch (source) {
            // change to your own paths
            case BCIF: archivePath = Paths.get("/var/bcif/"); break;
            case CIF: case MMCIF: archivePath = Paths.get("/var/pdb/"); break;
            // TODO obtain local MMTF archive - for now just pdbIds are obtained and then fetched from the RCSB PDB
            case MMTF: archivePath = Paths.get("/var/pdb/"); break;
        }

        Files.walk(archivePath)
                .parallel()
                // either process whole archive or limit to some number of structures
//                .limit(5000)
                .filter(path -> !Files.isDirectory(path))
                .forEach(path -> {
                    int count = counter.incrementAndGet();
                    if (count % chunkSize == 0) {
                        long end_chunk = System.nanoTime();
                        System.out.println("[" + count + "] @ " + (((end_chunk - start) /
                                1_000 / count) + " µs per structure"));
                    }

                    try {
                        // the work is to obtain the CifFile instance and convert into a BioJava structure
                        switch (source) {
                            case BCIF:
                                CifFileConverter.convert(CifReader.readBinary(Files.newInputStream(path))); break;
                            case CIF:
                                CifFileConverter.convert(CifReader.readText(new GZIPInputStream(Files.newInputStream(path)))); break;
                            case MMCIF:
                                new MMCIFFileReader().getStructure(new GZIPInputStream(Files.newInputStream(path))); break;
                            case MMTF:
                                // absolutely not comparable as files are fetched individually
                                new MMTFFileReader().getStructureById(path.toFile().getName().split("\\.")[0]); break;
                        }
                    } catch (Exception e) {
                        System.err.println("failed for " + path.toFile().getAbsolutePath());
                        e.printStackTrace();
                        failed.add(path.toFile().getName().split("\\.")[0]);
                    }
                });

        long end = System.nanoTime();
        System.out.println((end - start) / 1_000_000_000 + " s");
        System.out.println("failed for " + failed.size() + " structures");
        System.out.println("failed ids: " + failed);
    }
}