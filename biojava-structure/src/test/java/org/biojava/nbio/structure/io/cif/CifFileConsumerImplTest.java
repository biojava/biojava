package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.CifFileReader;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;
import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.IntColumn;
import org.rcsb.cif.model.ValueKind;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.MmCifFile;

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
    /**
     * Test parsing dates from MMCIF file version 4.
     */
    @Test
    public void testDatesV4() throws IOException, ParseException {
        InputStream inputStream = getClass().getResourceAsStream("/org/biojava/nbio/structure/io/mmcif/1stp_v4.cif");
        Objects.requireNonNull(inputStream, "could not acquire test resource /org/biojava/nbio/structure/io/mmcif/1stp_v4.cif");
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
        InputStream inputStream = getClass().getResourceAsStream("/org/biojava/nbio/structure/io/mmcif/1stp_v5.cif");
        Objects.requireNonNull(inputStream, "could not acquire test resource /org/biojava/nbio/structure/io/mmcif/1stp_v5.cif");
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
     * See issue <a href="https://github.com/biojava/biojava/issues/368">...</a>
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
        MmCifFile cifFile = CifIO.readFromInputStream(new ByteArrayInputStream(mmcifStr.getBytes())).as(StandardSchemata.MMCIF);
        IntColumn column = cifFile.getFirstBlock().getStructRefSeqDif().getSeqNum();

        assertNotNull(column);
        assertTrue(column.isDefined());
        assertEquals(8, column.getRowCount());
        column.valueKinds().forEach(vk -> assertEquals(ValueKind.NOT_PRESENT, vk));
        column.stringData().forEach(sd -> assertTrue(sd.isEmpty()));
    }

    /**
     * Testing files with atom_site that doesn't have author fields. E.g. cif files from Meta's ESM Atlas (<a href="https://esmatlas.com">...</a>)
     */
    @Test
    public void testAtomSiteWithMissingAuthFields() throws IOException {
        // taken from MGYP000911143359.cif
        String mmcifStr =
                "data_\n" +
                        "loop_\n" +
                        "_atom_site.group_PDB\n" +
                        "_atom_site.id\n" +
                        "_atom_site.type_symbol\n" +
                        "_atom_site.label_atom_id\n" +
                        "_atom_site.label_comp_id\n" +
                        "_atom_site.label_asym_id\n" +
                        "_atom_site.label_entity_id\n" +
                        "_atom_site.label_seq_id\n" +
                        "_atom_site.Cartn_x\n" +
                        "_atom_site.Cartn_y\n" +
                        "_atom_site.Cartn_z\n" +
                        "_atom_site.occupancy\n" +
                        "_atom_site.B_iso_or_equiv\n" +
                        "_atom_site.pdbx_PDB_model_num\n" +
                        "\n" +
                        "ATOM 1 N N MET A 1 1 -26.091 68.903 7.841 1.00 90.0 1\n" +
                        "ATOM 2 C CA MET A 1 1 -26.275 67.677 7.069 1.00 91.0 1\n" +
                        "ATOM 3 C C MET A 1 1 -24.933 67.025 6.755 1.00 90.0 1\n" +
                        "ATOM 4 C CB MET A 1 1 -27.033 67.967 5.773 1.00 89.0 1\n" +
                        "ATOM 5 O O MET A 1 1 -24.314 67.331 5.734 1.00 90.0 1\n" +
                        "ATOM 6 C CG MET A 1 1 -28.544 67.973 5.934 1.00 86.0 1\n" +
                        "ATOM 7 S SD MET A 1 1 -29.390 68.904 4.598 1.00 86.0 1\n" +
                        "ATOM 8 C CE MET A 1 1 -29.202 67.734 3.224 1.00 83.0 1\n" +
                        "ATOM 9 N N ASN A 1 2 -24.267 66.233 7.730 1.00 90.0 1\n" +
                        "ATOM 10 C CA ASN A 1 2 -22.897 65.827 8.029 1.00 91.0 1\n" +
                        "ATOM 11 C C ASN A 1 2 -22.600 64.427 7.500 1.00 90.0 1\n" +
                        "ATOM 12 C CB ASN A 1 2 -22.634 65.893 9.535 1.00 88.0 1\n" +
                        "ATOM 13 O O ASN A 1 2 -23.092 63.436 8.044 1.00 89.0 1\n" +
                        "ATOM 14 C CG ASN A 1 2 -22.191 67.269 9.990 1.00 86.0 1\n" +
                        "ATOM 15 N ND2 ASN A 1 2 -22.255 67.511 11.294 1.00 87.0 1\n" +
                        "ATOM 16 O OD1 ASN A 1 2 -21.795 68.108 9.177 1.00 87.0 1\n" ;
        MmCifFile cifFile = CifIO.readFromInputStream(new ByteArrayInputStream(mmcifStr.getBytes())).as(StandardSchemata.MMCIF);
        Structure s = CifStructureConverter.fromCifFile(cifFile);
        assertNotNull(s);
        assertEquals(2, s.getPolyChain("A").getAtomGroups().size());
        assertEquals(2, s.getPolyChainByPDB("A").getAtomGroups().size());
    }
}