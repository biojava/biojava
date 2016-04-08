package demo;

import java.io.IOException;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;
import org.biojava.nbio.structure.io.mmtf.MmtfUtils;

public class DemoMmtfReader {

	public static void main(String[] args) throws IOException, StructureException {
		MmtfUtils.setUpBioJava();
		// TODO Read in the byte array
		byte[] inputByteArray = new byte[0];
		MmtfActions.getBiojavaStruct(inputByteArray);
	}
	
}
