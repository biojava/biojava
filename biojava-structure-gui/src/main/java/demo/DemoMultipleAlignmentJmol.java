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
