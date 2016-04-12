package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;

/**
 * Class to show how to read a Biojava structure using MMTF
 * @author Anthony Bradley
 *
 */
public class DemoMmtfReader {

	public static void main(String[] args) throws IOException, StructureException {
		Structure structure = MmtfActions.readBiojavaStruct("/path/to/file");
		System.out.println(structure.getChains().size());
	}
	
}
