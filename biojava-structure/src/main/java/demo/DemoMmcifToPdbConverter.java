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

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/** 
 * An example of how to convert mmCIF file to PDB file
 *
 * @author Jose Duarte
 *
 */
public class DemoMmcifToPdbConverter
{

	public static void main(String[] args) throws Exception {

		File inFile = new File(args[0]);
		File outFile = new File(args[1]);
		convert(inFile, outFile);
	}



	public static void convert(File inFile, File outFile) throws IOException {
				 
        MMcifParser parser = new SimpleMMcifParser();
 
        SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();       
        parser.addMMcifConsumer(consumer);
        parser.parse(new BufferedReader(new InputStreamReader(new FileInputStream(inFile))));
        
        // now get the protein structure.
        Structure cifStructure = consumer.getStructure();

        // and write it out as PDB format
        PrintWriter pr = new PrintWriter(outFile);
        for (Chain c : cifStructure.getChains()) {
        		// we can override the chain name, the mmCIF chain names might have more than 1 character
        		c.setName(c.getName().substring(0, 1));
        		pr.print(c.toPDB());
        		pr.println("TER");
        }
        
		pr.close();
		

	}
}
