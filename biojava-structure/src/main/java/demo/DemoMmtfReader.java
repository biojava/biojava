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

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;

/**
 * Class to show how to read a Biojava structure using MMTF
 * @author Anthony Bradley
 *
 */
public class DemoMmtfReader {

	/**
	 * Main function to run the demo
	 * @param args no args to specify
	 * @throws IOException 
	 * @throws StructureException
	 */
	public static void main(String[] args) throws IOException, StructureException {
		Structure structure = MmtfActions.readFromWeb("4cup");
		System.out.println(structure.getChains().size());
	}
	
}
