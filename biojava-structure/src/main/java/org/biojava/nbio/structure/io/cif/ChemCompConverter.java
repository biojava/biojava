package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.ChemicalComponentDictionary;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.CifFile;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.MmCifBlock;
import org.rcsb.cif.schema.mm.MmCifFile;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Convert CifFiles to chem comps.
 * @author Sebastian Bittrich
 * @since 6.0.0
 */
public class ChemCompConverter {
    /**
     * Read data from a file and convert to chem comp dictionary.
     * @param path the source of information - can be gzipped or binary or text data
     * @return the target
     */
    public static ChemicalComponentDictionary fromPath(Path path) throws IOException {
        return fromInputStream(Files.newInputStream(path));
    }

    /**
     * Get data from a URL and convert to chem comp dictionary.
     * @param url the source of information - can be gzipped or binary or text data
     * @return the target
     * @throws IOException thrown when reading fails
     */
    public static ChemicalComponentDictionary fromURL(URL url) throws IOException {
        return fromInputStream(url.openStream());
    }

    /**
     * Convert InputStream to chem comp dictionary.
     * @param inputStream the InputStream of information - can be gzipped or binary or text data
     * @return the target
     * @throws IOException thrown when reading fails
     * @see CifStructureConverter#fromInputStream(InputStream, FileParsingParameters)
     */
    public static ChemicalComponentDictionary fromInputStream(InputStream inputStream) throws IOException {
        return fromCifFile(CifIO.readFromInputStream(inputStream));
    }

    /**
     * Convert CifFile to chem comp dictionary.
     * @param cifFile the source
     * @return the target
     */
    public static ChemicalComponentDictionary fromCifFile(CifFile cifFile) {
        // initialize consumer
        ChemCompConsumer consumer = new ChemCompConsumerImpl();

        // init structure
        consumer.prepare();

        // feed individual categories to consumer
        MmCifFile mmCifFile = cifFile.as(StandardSchemata.MMCIF);
        for (MmCifBlock cifBlock : mmCifFile.getBlocks()) {
            consumer.consumeChemComp(cifBlock.getChemComp());
            consumer.consumeChemCompAtom(cifBlock.getChemCompAtom());
            consumer.consumeChemCompBond(cifBlock.getChemCompBond());
        }

        // prepare structure to be retrieved
        consumer.finish();

        return consumer.getContainer();
    }
}
