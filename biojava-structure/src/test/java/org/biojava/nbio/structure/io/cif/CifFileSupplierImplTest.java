package org.biojava.nbio.structure.io.cif;

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

        InputStream inputStream = new ByteArrayInputStream(cifText.getBytes());
        Structure readStruct = CifStructureConverter.fromInputStream(inputStream);

        assertEquals(s.getEntityInfos().size(), readStruct.getEntityInfos().size());
        for (int i=0; i<s.getEntityInfos().size(); i++) {
            assertEquals(s.getEntityInfos().get(i).getMolId(), readStruct.getEntityInfos().get(i).getMolId());
            assertEquals(s.getEntityInfos().get(i).getType(), readStruct.getEntityInfos().get(i).getType());
        }

    }
}
