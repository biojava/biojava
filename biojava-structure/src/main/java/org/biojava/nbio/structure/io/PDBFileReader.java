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


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Compound;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

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
 public {@link Structure} loadStructure(String pathToPDBFile){
		{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();

		{@link Structure} structure = null;
		try{
			structure = pdbreader.getStructure(pathToPDBFile);
			System.out.println(structure);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return structure;
	}
 </pre>
 *
 * Access PDB files from a directory, take care of compressed PDB files
 * <pre>
 * public {@link Structure} loadStructureById() {
		String path = "/path/to/PDB/directory/";

		{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();
		pdbreader.setPath(path);
		{@link Structure} structure = null;
		try {
			structure = pdbreader.getStructureById("5pti");
		} catch (IOException e){
			e.printStackTrace();
		}
		return structure;

	}
	</pre>
 *
 *
 * @author Andreas Prlic
 *
 */
public class PDBFileReader extends LocalPDBDirectory {

	//private static final Logger logger = LoggerFactory.getLogger(PDBFileReader.class);

	// a list of big pdb files for testing
	//  "1htq",
	//  "1c2w",
	//  "1ffk",
	//  "1giy",
	//  "1j5a",
	//  "1jj2",
	//  "1jzx",
	//  "1jzy",
	//  "1jzz",
	//  "1k01",
	//  "1k73",
	//  "1k8a",
	//  "1k9m",
	//  "1kc8",
	//  "1kd1",
	//  "1kqs",
	//  "1m1k",
	//  "1m90",
	//  "1mkz",
	//  "1ml5",
	//  "1n8r",

	public static final String LOAD_CHEM_COMP_PROPERTY = "loadChemCompInfo";

	public static final String[] PDB_SPLIT_DIR    = new String[]{"data","structures","divided" ,"pdb"};
	public static final String[] PDB_OBSOLETE_DIR = new String[]{"data","structures","obsolete","pdb"};

	public static void main(String[] args){


		PDBFileReader pdbreader = new PDBFileReader();


		// set the path to cache files
		String property = "java.io.tmpdir";
		String tempdir = System.getProperty(property);
		// tempdir = "/path/to/local/PDB/installation/";
		pdbreader.setPath(tempdir);


		FileParsingParameters params = new FileParsingParameters();
		pdbreader.setFileParsingParameters(params);


		try{

			Structure struc = pdbreader.getStructureById("193D");
			System.out.println(struc);

			List<Compound>	compounds = struc.getCompounds();
			for (Compound comp : compounds  ){
				List<Chain> chains = comp.getChains();
				System.out.print(">Chains :" );
				for (Chain c : chains){
					System.out.print(c.getChainID() + " " );					
				}
				System.out.println();
				if ( chains.size() > 0)	{				
					System.out.println(chains.get(0).getAtomSequence());
					System.out.println(chains.get(0).getSeqResSequence());
					System.out.print("  Atom Ligands: ");

					for ( Group g: chains.get(0).getAtomLigands()){
						System.out.print( g.getPDBName() + " ");
					}

					System.out.println(" ");
				}
			}


			/*
			 GroupIterator gi = new GroupIterator(struc);
			while (gi.hasNext()){
				Group g = (Group) gi.next();
				Chain  c = g.getParent();
				if ( g instanceof AminoAcid ){

					AminoAcid aa = (AminoAcid)g;

					Map<String,String> sec = aa.getSecStruc();

					//System.out.println(c.getName() + " " + g + " " + sec);

					ChemComp cc = g.getChemComp();

					System.out.println(c.getName() + " " + g.getPDBCode() + " " + g.getPDBName() + " " + cc + " " +sec);
				}

			}
			 */
		} catch (IOException e) {
			e.printStackTrace();
		}
	}





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


	@Deprecated
	public void downloadPDB(String pdbId) throws IOException {
		downloadStructure(pdbId);
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchFileEvenIfObsolete the fetchFileEvenIfObsolete to set
	 * @deprecated Use {@link FileParsingParameters#setObsoleteBehavior(ObsoleteBehavior)}
	 */
	@Deprecated
	public void setFetchFileEvenIfObsolete(boolean fetchFileEvenIfObsolete) {
		if(fetchFileEvenIfObsolete) {
			setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		} else {
			if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE) {
				setObsoleteBehavior(ObsoleteBehavior.DEFAULT);
			}
		}
	}

	/**forces the reader to fetch the file if its status is OBSOLETE.
	 * This feature has a higher priority than {@link #setFetchCurrent(boolean)}. <br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchFileEvenIfObsolete
	 * @author Amr AL-Hossary
	 * @see #fetchCurrent
	 * @since 3.0.2
	 * @deprecated Use {@link FileParsingParameters#getObsoleteBehavior()}
	 */
	@Deprecated
	public boolean isFetchFileEvenIfObsolete() {
		return getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE;
	}


	/**if enabled, the reader searches for the newest possible PDB ID, if not present in he local installation.
	 * The {@link #setFetchFileEvenIfObsolete(boolean)} function has a higher priority than this function. <br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchCurrent the fetchCurrent to set
	 * @author Amr AL-Hossary
	 * @see #setFetchFileEvenIfObsolete(boolean)
	 * @since 3.0.2
	 * @deprecated Use {@link FileParsingParameters#setObsoleteBehavior(ObsoleteBehavior)}
	 */
	@Deprecated
	public void setFetchCurrent(boolean fetchNewestCurrent) {
		if(fetchNewestCurrent) {
			setObsoleteBehavior(ObsoleteBehavior.FETCH_CURRENT);
		} else {
			if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_CURRENT) {
				setObsoleteBehavior(ObsoleteBehavior.DEFAULT);
			}
		}
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchCurrent
	 * @deprecated Use {@link FileParsingParameters#getObsoleteBehavior()}
	 */
	@Deprecated
	public boolean isFetchCurrent() {
		return getObsoleteBehavior() == ObsoleteBehavior.FETCH_CURRENT;
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
