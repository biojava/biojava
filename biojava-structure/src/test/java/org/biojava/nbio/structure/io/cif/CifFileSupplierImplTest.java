package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileParser;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

public class CifFileSupplierImplTest {

    @Test
    public void shouldReadRawPdbOutputtingCifWithEntity() throws IOException {
        InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4lup_phaser_output.pdb.gz"));

        PDBFileParser pdbpars = new PDBFileParser();
        FileParsingParameters params = new FileParsingParameters();
        params.setAlignSeqRes(true);
        pdbpars.setFileParsingParameters(params);

        Structure s = pdbpars.parsePDBFile(inStream);

        String cifText = CifStructureConverter.toText(s);
        assertTrue(cifText.contains("_entity.type"));
        assertTrue(cifText.contains("_entity_poly.pdbx_seq_one_letter_code_can"));
        assertFalse(cifText.contains("null"));
        assertTrue(cifText.contains("MSEQLTDQVLVERVQKGDQKAFNLLVVRYQHKVASLVSRYVPSGDVPDVVQEAFIKA"));

        InputStream inputStream = new ByteArrayInputStream(cifText.getBytes());
        Structure readStruct = CifStructureConverter.fromInputStream(inputStream);

        assertEquals(s.getEntityInfos().size(), readStruct.getEntityInfos().size());
        for (int i=0; i<s.getEntityInfos().size(); i++) {
            assertEquals(s.getEntityInfos().get(i).getMolId(), readStruct.getEntityInfos().get(i).getMolId());
            assertEquals(s.getEntityInfos().get(i).getType(), readStruct.getEntityInfos().get(i).getType());
        }

    }

    /**
     * The identifier must be written as a data item and not only as the name of the data block: consumers read it
     * from _entry.id or from _struct.entry_id, so writing the block header alone loses it. See issue #1143.
     */
    @Test
    public void shouldWriteEntryIdAndSurviveRoundTrip() throws IOException {
        Structure s;
        try (InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.cif.gz"))) {
            s = CifStructureConverter.fromInputStream(inStream);
        }
        assertEquals(new PdbId("4HHB"), s.getPdbId());

        String cifText = CifStructureConverter.toText(s);
        assertTrue("_entry.id must be written", cifText.contains("_entry.id"));
        assertTrue("_struct.entry_id must be written", cifText.contains("_struct.entry_id"));

        Structure readStruct = CifStructureConverter.fromInputStream(
                new ByteArrayInputStream(cifText.getBytes()));

        assertEquals(s.getPdbId(), readStruct.getPdbId());
        assertEquals(s.getPdbId(), readStruct.getPDBHeader().getPdbId());
    }

    /**
     * Structures without an identifier must not gain empty entry categories.
     */
    @Test
    public void shouldNotWriteEntryIdWhenPdbIdIsAbsent() throws IOException {
        Structure s;
        try (InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/4hhb.cif.gz"))) {
            s = CifStructureConverter.fromInputStream(inStream);
        }
        s.setPdbId(null);

        String cifText = CifStructureConverter.toText(s);
        assertFalse(cifText.contains("_entry.id"));
        assertFalse(cifText.contains("_struct.entry_id"));
    }
}
