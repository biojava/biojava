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
package demo;

import java.util.Map;
import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class DemoLoadSecStruc {
    public static void main(String[] args){

        try {
            FileParsingParameters params = new FileParsingParameters();
            params.setParseSecStruc(true);

            AtomCache cache = new AtomCache();
            cache.setFileParsingParams(params);

            Structure s = cache.getStructure("4hhb");

            for ( Chain c : s.getChains()) {
                for (Group g: c.getAtomGroups()){

                    if ( g instanceof AminoAcid ){

                        AminoAcid aa = (AminoAcid)g;

                        Map<String,String> sec = aa.getSecStruc();

                        System.out.println(c.getChainID() + " " + g.getResidueNumber() + " " + g.getPDBName() + " " + " " +sec);
                    }
                }
            }

        } catch (Exception e) {

            e.printStackTrace();
        }        
    }
}