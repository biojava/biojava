package org.biojava.nbio.structure.io;

import java.io.IOException;
import java.io.InputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmtf.MmtfActions;

/**
 * A class to read MMTF files and cache them locally.
 * @author Anthony Bradley
 *
 */
public class MMTFFileReader extends LocalPDBDirectory {
	
	
	public static final String[] MMTF_SPLIT_DIR    = new String[]{"data","structures","divided" ,"mmtf"};
	public static final String[] MMTF_OBSOLETE_DIR = new String[]{"data","structures","obsolete","mmtf"};
	
	public static void main(String[] args) throws Exception {

		MMTFFileReader reader = new MMTFFileReader();
		FileParsingParameters params = new FileParsingParameters();
		reader.setFileParsingParameters(params);
		Structure struc = reader.getStructureById("1m4x");
		System.out.println(struc);
	}
	
	
	@Override
	public Structure getStructure(InputStream inStream) throws IOException {
		return MmtfActions.readFromInputStream(inStream);
	}

	@Override
	protected String getFilename(String pdbId) {
		return pdbId.toLowerCase()+".mmtf.gz";
	}

	@Override
	protected String[] getSplitDirPath() {
		return MMTF_SPLIT_DIR;
	}

	@Override
	protected String[] getObsoleteDirPath() {
		return MMTF_OBSOLETE_DIR;
	}
	

}
