package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.*;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.Test;
import org.rcsb.cif.CifReader;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.model.Column;
import org.rcsb.cif.model.ValueKind;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

public class CifFileConsumerImplTest {
    // TODO use test resources provided by integration-test module

    @Test
    public void testEntityId() throws IOException, StructureException {
        // Set up the atom cache to parse on Internal chain id
        AtomCache cache = new AtomCache();
        cache.setUseCif(true);
        FileParsingParameters params = cache.getFileParsingParams();

        DownloadChemCompProvider cc = new DownloadChemCompProvider();
        ChemCompGroupFactory.setChemCompProvider(cc);
        cc.checkDoFirstInstall();
        cache.setFileParsingParams(params);
        StructureIO.setAtomCache(cache);
        // This is hte information we want to test against
        String[] typeInformation = new String[] {"POLYMER", "NONPOLYMER", "NONPOLYMER", "NONPOLYMER", "NONPOLYMER", "WATER"};
        String[] descriptionInformation = new String[] {"BROMODOMAIN ADJACENT TO ZINC FINGER DOMAIN PROTEIN 2B","4-Fluorobenzamidoxime",  "METHANOL", "METHANOL", "METHANOL", "water"};

        // Now some other information fields to test this data is collated correctly
        String[] geneSourceSciName = new String[] {"HOMO SAPIENS", null, null, null, null, null};
        String[] geneSourceTaxId = new String[] {"9606", null, null, null, null, null};
        String[] hostOrganismSciName = new String[] {"ESCHERICHIA COLI", null, null, null, null, null};
        String[] hostOrganismTaxId = new String[] {"469008", null, null, null, null, null};

        /// TODO GET ALL THE ENTITY INFORMATION REQUIRED FOR 4CUP
        // Get this structure
        Structure bioJavaStruct = StructureIO.getStructure("4cup");
        String[] testTypeInfo = new String[6];
        String[] testDescInfo = new String[6];

        String[] testGeneSourceSciName = new String[6];
        String[] testGeneSourceTaxId = new String[6];
        String[] testHostOrganismSciName = new String[6];
        String[] testHostOrganismTaxId = new String[6];

        // Now loop through the structure
        int chainCounter = 0;
        for (Chain c: bioJavaStruct.getChains()) {
            // Now get the entity information we want to test
            EntityInfo thisCmpd = c.getEntityInfo();
            testTypeInfo[chainCounter] = thisCmpd.getType().toString();
            testDescInfo[chainCounter] = thisCmpd.getDescription();
            testGeneSourceSciName[chainCounter] =  thisCmpd.getOrganismScientific();
            testGeneSourceTaxId[chainCounter] = thisCmpd.getOrganismTaxId();
            testHostOrganismSciName[chainCounter] = thisCmpd.getExpressionSystem();
            testHostOrganismTaxId[chainCounter] = thisCmpd.getExpressionSystemTaxId();

            chainCounter++;
        }
        // Now check they're both the same
        assertArrayEquals(descriptionInformation, testDescInfo);
        assertArrayEquals(typeInformation, testTypeInfo);
        // Now check these work too
        assertArrayEquals(geneSourceSciName, testGeneSourceSciName);
        assertArrayEquals(geneSourceTaxId, testGeneSourceTaxId);
        assertArrayEquals(hostOrganismSciName, testHostOrganismSciName);
        assertArrayEquals(hostOrganismTaxId, testHostOrganismTaxId);
    }

    /**
     * Test parsing dates from MMCIF file version 4.
     */
    @Test
    public void testDatesV4() throws IOException, ParseException {
        InputStream inputStream = getClass().getResourceAsStream("org/biojava/nbio/structure/io/mmcif/1stp_v4.cif");
        Objects.requireNonNull(inputStream, "could not acquire test resource org/biojava/nbio/structure/io/mmcif/1stp_v4.cif");
        Structure s = new CifFileReader().getStructure(inputStream);

        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd", Locale.US);

        Date modDate = dateFormat.parse("2011-07-13");
        assertEquals(modDate, s.getPDBHeader().getModDate());

        Date releaseDate = dateFormat.parse("1992-10-15");
        assertEquals(releaseDate, s.getPDBHeader().getRelDate());

        Date depositionDate = dateFormat.parse("1992-03-12");
        assertEquals(depositionDate, s.getPDBHeader().getDepDate());
    }

    /**
     * Test parsing dates from MMCIF file version 5.
     */
    @Test
    public void testDatesV5() throws IOException, ParseException {
        InputStream inputStream = getClass().getResourceAsStream("org/biojava/nbio/structure/io/mmcif/1stp_v5.cif");
        Objects.requireNonNull(inputStream, "could not acquire test resource org/biojava/nbio/structure/io/mmcif/1stp_v5.cif");
        Structure s = new CifFileReader().getStructure(inputStream);

        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd", Locale.US);

        Date modDate = dateFormat.parse("2011-07-13");
        assertEquals(modDate, s.getPDBHeader().getModDate());

        Date releaseDate = dateFormat.parse("1992-10-15");
        assertEquals(releaseDate, s.getPDBHeader().getRelDate());

        Date depositionDate = dateFormat.parse("1992-03-12");
        assertEquals(depositionDate, s.getPDBHeader().getDepDate());
    }

    /**
     * A test for reading a phenix-produced (ver 1.9_1692) mmCIF file. This is the file submitted to the PDB for
     * deposition of entry 4lup - See github issue #234
     * @throws IOException propagated
     */
    @Test
    public void testPhenixCifFile() throws IOException {
        InputStream inputStream = new GZIPInputStream(getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4lup_phenix_output.cif.gz"));

        Structure structure = new CifFileReader().getStructure(inputStream);

        assertNotNull(structure);
        assertTrue(structure.isCrystallographic());

        // all ligands are into their own chains, so we have 2 proteins, 2 nucleotide chains, 1 ligand chain and 1 purely water chain
        assertEquals(6, structure.getChains().size());

        // 4 entities: 1 protein, 1 nucleotide, 1 water, 1 ligand (EDO)
        assertEquals(4, structure.getEntityInfos().size());
        int[] counts = countEntityTypes(structure.getEntityInfos());
        assertEquals(2, counts[0]);
        assertEquals(1, counts[1]);
        assertEquals(1, counts[2]);
    }

    /**
     * This test represents a common situation for a non-deposited structure.
     * When building with common crystallography software, the user often adds new
     * ligands (or solvent) molecules as new chains.  Only prior to deposition
     * then relabel them so that they belong to the same chain as the polymeric residues.
     *
     * In this case, the ligands represent valuable information and should not be discarded.
     */
    @Test
    public void testNewLigandChain() throws IOException {
        // Test the file parsing speed when the files are already downloaded.
        InputStream pdbStream = new GZIPInputStream(getClass().getResourceAsStream("/ligandTest.pdb.gz"));
        InputStream cifStream = new GZIPInputStream(getClass().getResourceAsStream("/ligandTest.cif.gz"));

        assertNotNull(cifStream);
        assertNotNull(pdbStream);

        FileParsingParameters params = new FileParsingParameters();
        PDBFileParser pdbpars = new PDBFileParser();
        pdbpars.setFileParsingParameters(params);
        Structure s1 = pdbpars.parsePDBFile(pdbStream) ;

        // The chain B should be present with 1 ligand HEM
        Chain c1 = s1.getNonPolyChainsByPDB("B").get(0);
        assertNotNull(c1);

        int expectedNumLigands = 1;
        assertEquals(expectedNumLigands, c1.getAtomGroups().size());

        Structure s2 = new CifFileReader().getStructure(cifStream);

        // The chain B should be present with 1 ligand HEM
        Chain c2 = s2.getNonPolyChainsByPDB("B").get(0);
        assertNotNull(c2);
        assertEquals(expectedNumLigands, c2.getAtomGroups().size());

        // pdb and mmcif should have same number of chains
        assertEquals(s1.getChains().size(), s2.getChains().size());
    }

    @Test
    public void testWaterOnlyChainCif() throws IOException {

        // following file is cut-down versions of 4a10
        InputStream cifStream = new GZIPInputStream(getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4a10_short.cif.gz"));

        Structure s2 = new CifFileReader().getStructure(cifStream);

        assertEquals(2, s2.getChains().size());

        Chain c = s2.getWaterChainByPDB("F");

        assertNotNull("Got null when looking for water-only chain with author id F", c);
        assertTrue(c.getAtomGroups().size() > 0);

        // checking that compounds are linked
        assertNotNull(c.getEntityInfo());

        // checking that the water molecule was assigned an ad-hoc compound
        assertEquals(2, s2.getEntityInfos().size());

        Chain cAsymId = s2.getWaterChain("E");
        assertNotNull("Got null when looking for water-only chain with asym id E", cAsymId);
        assertTrue(cAsymId.getAtomGroups().size() > 0);
        assertSame(c, cAsymId);
    }

    private static int[] countEntityTypes(List<EntityInfo> entities) {
        int countPoly = 0;
        int countNonPoly = 0;
        int countWater = 0;
        for (EntityInfo e : entities) {
            if (e.getType() == EntityType.POLYMER) {
                countPoly++;
            }
            if (e.getType() == EntityType.NONPOLYMER) {
                countNonPoly++;
            }
            if (e.getType() == EntityType.WATER) {
                countWater++;
            }
        }
        return new int[] { countPoly, countNonPoly, countWater };
    }

    /**
     * This tests for cases where dots appear in integer fields. Unusual but it happens in some PDB entries like 1s32.
     * See issue https://github.com/biojava/biojava/issues/368
     */
    @Test
    public void specialCases() throws IOException {
        // taken from 1s32
        String mmcifStr =
                "data_\n" +
                "loop_\n" +
                "_struct_ref_seq_dif.align_id\n" +
                "_struct_ref_seq_dif.pdbx_pdb_id_code\n"+
                "_struct_ref_seq_dif.mon_id\n"+
                "_struct_ref_seq_dif.pdbx_pdb_strand_id\n"+
                "_struct_ref_seq_dif.seq_num\n"+ // integer field that contains '.'
                "_struct_ref_seq_dif.pdbx_seq_db_name\n"+
                "_struct_ref_seq_dif.pdbx_seq_db_accession_code\n"+
                "_struct_ref_seq_dif.db_mon_id\n"+
                "_struct_ref_seq_dif.pdbx_seq_db_seq_num\n"+
                "_struct_ref_seq_dif.details\n"+
                "_struct_ref_seq_dif.pdbx_auth_seq_num\n"+
                "_struct_ref_seq_dif.pdbx_pdb_ins_code\n"+
                "_struct_ref_seq_dif.pdbx_ordinal\n"+
                "1 1S32 . A . GB  30268544 MET 1 'INTIATING METHIONINE' ? ? 1\n"+
                "2 1S32 . E . GB  30268544 MET 1 'INTIATING METHIONINE' ? ? 2\n"+
                "3 1S32 . B . UNP P02304   MET 0 'INTIATING METHIONINE' ? ? 3\n"+
                "4 1S32 . F . UNP P02304   MET 0 'INTIATING METHIONINE' ? ? 4\n"+
                "5 1S32 . C . GB  30268540 MET 1 'INTIATING METHIONINE' ? ? 5\n"+
                "6 1S32 . G . GB  30268540 MET 1 'INTIATING METHIONINE' ? ? 6\n"+
                "7 1S32 . D . GB  30268542 MET 1 'INTIATING METHIONINE' ? ? 7\n"+
                "8 1S32 . H . GB  30268542 MET 1 'INTIATING METHIONINE' ? ? 8\n" +
                "#" ;
        CifFile cifFile = CifReader.parseText(new ByteArrayInputStream(mmcifStr.getBytes()));
        Column column = cifFile.getFirstBlock().getCategory("struct_ref_seq_dif").getColumn("seq_num");

        assertNotNull(column);
        assertTrue(column.isDefined());
        assertEquals(8, column.getRowCount());
        column.valueKinds().forEach(vk -> assertEquals(ValueKind.NOT_PRESENT, vk));
        column.getStringData().forEach(sd -> assertTrue(sd.isEmpty()));
    }

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
        InputStream inStream = this.getClass().getResourceAsStream(fileName);
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

        if ( pdbStructure.isNmr()){
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

        if (g1.has3D()){
            Atom a1 = g1.getAtom(0);
            Atom a2 = g2.getAtom(0);
            if ( a1 == null)
                fail("could not get atom for group " + g1);
            if (a2 == null)
                fail("could not get atom for group " + g2);
            assertEquals(a1.getX(),a2.getX(), 0.0001);
            assertEquals(a1.getOccupancy(), a2.getOccupancy(), 0.0001);
            assertEquals(a1.getTempFactor(), a2.getTempFactor(), 0.0001);
            assertEquals(a1.getName(), a2.getName());
        }
    }

    private void checkNMR(Structure s){
        assertTrue(s.isNmr());
        int models = s.nrModels();
        assertTrue(models > 0);
        List<Chain> model0 = s.getModel(0);

        // compare with all others
        for (int i = 1 ; i < models; i++){
            List<Chain> modelX = s.getModel(i);
            assertEquals(model0.size(),modelX.size());

            // compare lengths:
            for (int j=0 ; j< model0.size(); j++){
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