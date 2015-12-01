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
import java.util.List;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.secstruc.SecStrucElement;
import org.biojava.nbio.structure.secstruc.SecStrucPred;
import org.biojava.nbio.structure.secstruc.SecStrucTools;

/**
 * Demonstration on how to use the Secondary Structure Prediction (DSSP)
 * implementation in BioJava and obtain different SS representations and
 * outputs.
 * 
 * @author Aleix Lafita
 *
 */
public class DemoSecStrucPred {

	public static void main(String[] args) throws IOException,
			StructureException {

		String pdbID = "5pti";

		AtomCache cache = new AtomCache();

		// Load structure without any SS assignment
		Structure s = cache.getStructure(pdbID);

		// Predict and assign the SS of the Structure
		SecStrucPred ssp = new SecStrucPred();
		ssp.predict(s, true);

		// Print the DSSP output
		System.out.println("******DSSP output: ");
		System.out.println(ssp.printDSSP());

		// Print the FASTA sequence of SS
		System.out.println("\n******FASTA output: ");
		System.out.println(ssp.printFASTA());

		// Print the Helix Summary
		System.out.println("\n******Helix Summary: ");
		System.out.println(ssp.printHelixSummary());

		// Obtain and print the SS elements of the Structure
		List<SecStrucElement> sse = SecStrucTools.getSecStrucElements(s);
		System.out.println("\n******SecStrucElements: ");
		for (SecStrucElement e : sse)
			System.out.println(e);

	}
}
