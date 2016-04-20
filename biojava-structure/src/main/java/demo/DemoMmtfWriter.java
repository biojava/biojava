package demo;

import java.io.IOException;
import java.nio.file.Paths;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;

public class DemoMmtfWriter {

	public static void main(String[] args) throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		Structure structure = StructureIO.getStructure("4cup");
		MmtfActions.writeToFile(structure, Paths.get("/tmp/4cup.mmtf"));
	}
	
	
}
