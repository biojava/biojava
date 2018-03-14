/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.io.mmcif;

import org.biojava.nbio.structure.io.mmcif.chem.MetalBondDistance;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.InputStream;
import java.util.*;

import java.util.zip.GZIPInputStream;

/**
 * Created by andreas on 6/6/16.
 */
public class MetalBondParser  {

    private static final Logger logger = LoggerFactory.getLogger(MetalBondParser.class);

    private static final String BONDS_FILE = "org/biojava/nbio/structure/bond_distance_limits.cif.gz";


    static Map<String,List<MetalBondDistance>> definitions;

    static {
         definitions = init();
    }


    public static Map<String,List<MetalBondDistance>> getMetalBondDefinitions(){
        return definitions;

    }


    private static Map<String,List<MetalBondDistance>> init(){

        InputStream inputStream = MetalBondParser.class.getClassLoader().getResourceAsStream(BONDS_FILE);

        if (inputStream == null) {
            throw new RuntimeException("Could not find resource "+BONDS_FILE+".  This probably means that your biojava.jar file is corrupt or incorrectly built.");
        }

        try {
            GZIPInputStream gzIS = new GZIPInputStream(inputStream);

            SimpleMMcifParser parser = new SimpleMMcifParser();

            MetalBondConsumer consumer = new MetalBondConsumer();
            parser.addMMcifConsumer(consumer);

            parser.parse(gzIS);

            Map<String,List<MetalBondDistance>> defs = consumer.getDefinitions();

            return defs;

        } catch ( Exception e){
            logger.error(e.getMessage(),e);

        }
        return null;
    }



}
