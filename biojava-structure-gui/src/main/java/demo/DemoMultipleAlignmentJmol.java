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
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.TestSampleGenerator;

/**
 * Demo for visualizing the results of a Multiple Alignment, 
 * from a sample MultipleAlignment object.
 * 
 * @author Aleix Lafita
 * 
 */
public class DemoMultipleAlignmentJmol {

	public static void main(String[] args)
			throws StructureException, IOException{
		
		//Obtain a sample MultipleAlignment from the tests
		MultipleAlignment msa = TestSampleGenerator.testAlignment2();
		
		//Displaying a MultipleAlignment is a one line code
		MultipleAlignmentDisplay.display(msa);
	}
}
