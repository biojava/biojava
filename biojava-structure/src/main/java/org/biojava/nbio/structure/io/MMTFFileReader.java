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
package org.biojava.nbio.structure.io;

import java.io.IOException;
import java.io.InputStream;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.UserConfiguration;
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
	
	/**
	 * Constructs a new {@link MMTFFileReader}, initializing the extensions member variable.
	 * The path is initialized in the same way as {@link UserConfiguration},
	 * i.e. to system property/environment variable {@link UserConfiguration#PDB_DIR}.
	 * Both autoFetch and splitDir are initialized to false
	 */
	public MMTFFileReader() {
		this(null);
	}
	
	/**
	 * Constructs a new {@link MMTFFileReader}, initializing the extensions member variable.
	 * The path is initialized to the given path, both autoFetch and splitDir are initialized to false.
	 */
	public MMTFFileReader(String path) {
		super(path);
		addExtension(".mmtf");
		addExtension(".mmtf.gz");	
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
