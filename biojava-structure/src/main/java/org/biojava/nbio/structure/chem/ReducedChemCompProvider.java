package org.biojava.nbio.structure.chem;

import org.biojava.nbio.structure.io.cif.ChemCompConsumer;
import org.biojava.nbio.structure.io.cif.ChemCompConsumerImpl;
import org.biojava.nbio.structure.io.cif.ChemCompConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 * Unlike the {@link DownloadChemCompProvider}, this  {@link ChemCompProvider} does not download any chem comp
 * definitions. It has access to a limited set of files that are part of the biojava distribution.
 *
 * @author Andreas Prlic
 * @since 3.0
 */
public class ReducedChemCompProvider implements ChemCompProvider {
    private static final Logger logger = LoggerFactory.getLogger(ReducedChemCompProvider.class);

    public ReducedChemCompProvider(){
        logger.debug("Initialising ReducedChemCompProvider");
    }

    @Override
    public ChemComp getChemComp(String recordName) {
        String name = recordName.toUpperCase().trim();
        try (InputStream inStream = this.getClass().getResourceAsStream("/chemcomp/" + name + ".cif.gz")) {
            logger.debug("Reading chemcomp/{}.cif.gz", recordName);

            if (inStream == null) {
                //System.out.println("Could not find chem comp: " + name + " ... using generic Chem Comp");
                // could not find the chem comp definition for this in the jar file
                logger.debug("Getting empty chem comp for {}", name);
                ChemComp cc = ChemComp.getEmptyChemComp();
                cc.setId(name);
                return cc;
            }

            // The Consumer builds up the BioJava - structure object.
            // you could also hook in your own and build up you own data model.
            ChemicalComponentDictionary dict = ChemCompConverter.fromInputStream(inStream);

            return dict.getChemComp(name);
        } catch (IOException e) {
            logger.error("IOException caught while reading chem comp {}.", name, e);
        }
        logger.warn("Problem when loading chem comp {}, will use an empty chem comp for it", name);
        ChemComp cc = ChemComp.getEmptyChemComp();
        cc.setId(name);
        return cc;
    }
}

