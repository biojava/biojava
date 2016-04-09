package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;

public class DemoMmtfRoundTrip {

	public static void main(String[] args) throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		Structure structure = StructureIO.getStructure("/Users/abradley/Downloads/4cup.cif");
		// We can do somme comparisons on the round tripped structure
		MmtfActions.roundTrip(structure);
	}
	

}
