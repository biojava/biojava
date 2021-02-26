package org.biojava.nbio.structure.io.cif;

import org.biojava.nbio.structure.chem.MetalBondDistance;
import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.Block;
import org.rcsb.cif.model.CifFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.InputStream;
import java.util.List;
import java.util.Map;

/**
 * Created by andreas on 6/6/16.
 */
public class MetalBondConverter {
    private static final Logger logger = LoggerFactory.getLogger(MetalBondConverter.class);
    private static final String BONDS_FILE = "org/biojava/nbio/structure/bond_distance_limits.cif.gz";
    private static final Map<String, List<MetalBondDistance>> definitions;

    static {
        definitions = init();
    }

    public static Map<String,List<MetalBondDistance>> getMetalBondDefinitions() {
        return definitions;
    }

    private static Map<String,List<MetalBondDistance>> init() {
        InputStream inputStream = MetalBondConverter.class.getClassLoader().getResourceAsStream(BONDS_FILE);

        if (inputStream == null) {
            throw new RuntimeException("Could not find resource " + BONDS_FILE + ".  This probably means that your " +
                    "biojava.jar file is corrupt or incorrectly built.");
        }

        try {
            CifFile cifFile = CifIO.readFromInputStream(inputStream);
            // initialize consumer
            MetalBondConsumerImpl consumer = new MetalBondConsumerImpl();

            // init structure
            consumer.prepare();

            // feed individual categories to consumer
            for (Block cifBlock : cifFile.getBlocks()) {
                cifBlock.categories().forEach(consumer::consume);
            }

            // prepare structure to be retrieved
            consumer.finish();

            return consumer.getContainer();
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }
        return null;
    }
}
