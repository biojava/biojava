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
 * Created on 16.03.2004
 * @author Andreas Prlic
 *
 *
 */
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import java.io.IOException;
import java.io.InputStream;

/**
 * <p>
 *  The wrapper class for parsing a PDB file.
 *  </p>
 *
 *
 *  <p>
 *  Several flags can be set for this class
 *  <ul>
 *
 * <li> {@link #setAutoFetch(boolean)} - if the PDB file can not be found locally, should it be fetched
 *  from the PDB ftp servers? (default:false)</li>
 *  <li> Other parameters can be set using the {@link #setFileParsingParameters(FileParsingParameters)}</li>
 *  </ul>
 *  </p>
 *
 *
 *
 *<h2>Example</h2>
 * <p>
 * Q: How can I get a Structure object from a PDB file?
 * </p>
 * <p>
 * A:
 * <pre>
 * public {@link Structure} loadStructure(String pathToPDBFile){
 * 	{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();
 *
 * 	{@link Structure} structure = null;
 * 	try{
 * 		structure = pdbreader.getStructure(pathToPDBFile);
 * 		System.out.println(structure);
 * 	} catch (IOException e) {
 * 		e.printStackTrace();
 * 	}
 * 	return structure;
 * }
 * </pre>
 *
 * Access PDB files from a directory, take care of compressed PDB files
 * <pre>
 * public {@link Structure} loadStructureById() {
 * 	String path = "/path/to/PDB/directory/";
 *
 * 	{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();
 * 	pdbreader.setPath(path);
 * 	{@link Structure} structure = null;
 * 	try {
 * 		structure = pdbreader.getStructureById("5pti");
 * 	} catch (IOException e){
 * 		e.printStackTrace();
 * 	}
 * 	return structure;
 *
 * }
 * </pre>
 *
 *
 * @author Andreas Prlic
 *
 */
public class PDBFileReader extends LocalPDBDirectory {

	//private static final Logger logger = LoggerFactory.getLogger(PDBFileReader.class);

	public static final String[] PDB_SPLIT_DIR    = new String[]{"data","structures","divided" ,"pdb"};
	public static final String[] PDB_OBSOLETE_DIR = new String[]{"data","structures","obsolete","pdb"};


	/**
	 * Constructs a new PDBFileReader, initializing the extensions member variable.
	 * The path is initialized in the same way as {@link UserConfiguration},
	 * i.e. to system property/environment variable {@link UserConfiguration#PDB_DIR}.
	 * Both autoFetch and splitDir are initialized to false
	 */
	public PDBFileReader() {
		this(null);
	}

	/**
	 * Constructs a new PDBFileReader, initializing the extensions member variable.
	 * The path is initialized to the given path, both autoFetch and splitDir are initialized to false.
	 *
	 * <p>If path is null, initialize using the system property/environment variable
	 * {@link UserConfiguration#PDB_DIR}.
	 * @param path Path to the PDB file directory
	 */
	public PDBFileReader(String path) {
		super(path);

		addExtension(".ent");
		addExtension(".pdb");
		addExtension(".ent.gz");
		addExtension(".pdb.gz");
		addExtension(".ent.Z");
		addExtension(".pdb.Z");
	}

	@Override
	protected String getFilename(String pdbId) {
		return "pdb"+pdbId.toLowerCase()+".ent.gz";
	}

	@Override
	public Structure getStructure(InputStream inStream) throws IOException {
		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setFileParsingParameters(getFileParsingParameters());

		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;
	}

	@Override
	protected String[] getSplitDirPath() {
		return PDB_SPLIT_DIR;
	}

	@Override
	protected String[] getObsoleteDirPath() {
		return PDB_OBSOLETE_DIR;
	}



}
