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
package org.biojava.bio.structure.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupIterator;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava.utils.io.InputStreamProvider;


/**
 * <p>
 *  The wrapper class for parsing a PDB file.
 *  </p>
 *
 *
 *  <p>
 *  Several flags can be set for this class
 *  <ul>
 *  <li> {@link #setParseCAOnly} - parse only the Atom records for C-alpha atoms (default:false)</li>
 * <li> {@link #setParseSecStruc} - a flag if the secondary structure information from the PDB file (author's assignment) should be parsed.
 *      If true the assignment can be accessed through {@link AminoAcid}.getSecStruc(); (default:false)</li>
 * <li> {@link #setAlignSeqRes(boolean)} - should the AminoAcid sequences from the SEQRES
 *      and ATOM records of a PDB file be aligned? (default:true)</li>
 * <li> {@link #setAutoFetch(boolean)} - if the PDB file can not be found locally, should it be fetched
 *  from the PDB ftp servers? (default:false)</li>
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
public class PDBFileReader implements StructureIOFile {

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


	String path                     ;
	List<String> extensions            ;
	boolean parseSecStruc;
	boolean autoFetch;
	boolean parseCAOnly;
	boolean alignSeqRes;
	boolean pdbDirectorySplit;
	boolean headerOnly;
	boolean loadChemCompInfo;

	public static final String lineSplit = System.getProperty("file.separator");

	public static final String DEFAULT_PDB_FILE_SERVER = "ftp.wwpdb.org";

	public static final String PDB_FILE_SERVER_PROPERTY = "PDB.FILE.SERVER";

	public static void main(String[] args){


		PDBFileReader pdbreader = new PDBFileReader();


		// set the path to cache files
		String property = "java.io.tmpdir";
		String tempdir = System.getProperty(property);
		// tempdir = "/path/to/local/PDB/installation/";
		pdbreader.setPath(tempdir);

		pdbreader.setParseSecStruc(true);
		pdbreader.setAlignSeqRes(true);
		pdbreader.setParseCAOnly(false);
		pdbreader.setAutoFetch(true);
		pdbreader.setLoadChemCompInfo(true);

		try{

			Structure struc = pdbreader.getStructureById("193D");
			System.out.println(struc);

			List<Compound>	compounds = struc.getCompounds();
			for (Compound comp : compounds  ){
				List<Chain> chains = comp.getChains();
				System.out.print(">Chains :" );
				for (Chain c : chains){
					System.out.print(c.getName() + " " );					
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
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	public PDBFileReader() {
		extensions    = new ArrayList<String>();
		path = "" ;
		extensions.add(".ent");
		extensions.add(".pdb");
		extensions.add(".ent.gz");
		extensions.add(".pdb.gz");
		extensions.add(".ent.Z");
		extensions.add(".pdb.Z");
		parseSecStruc = false;
		autoFetch     = false;
		parseCAOnly   = false;
		alignSeqRes   = true;
		pdbDirectorySplit = false;
		headerOnly    = false;
		loadChemCompInfo = false;
	}


	/** return the flag if only the CA atoms should be parsed
	 *
	 * @return flag if CA only should be read
	 */
	public boolean isParseCAOnly() {
		return parseCAOnly;
	}

	/** only the CA atoms should be parsed from the PDB file
	 *
	 * @param parseCAOnly
	 */
	public void setParseCAOnly(boolean parseCAOnly) {
		this.parseCAOnly = parseCAOnly;
	}

	/** get the flag if the SEQRES and ATOM amino acids are going to be aligned
	 *
	 * @return flag
	 */
	public boolean isAlignSeqRes() {
		return alignSeqRes;
	}


	/** set the flag if the SEQRES and ATOM amino acids should be aligned and linked
	 *
	 * @param alignSeqRes
	 */
	public void setAlignSeqRes(boolean alignSeqRes) {
		this.alignSeqRes = alignSeqRes;
	}


	/** should the parser to fetch missing PDB files from the PDB FTP server automatically?
	 *  default is false
	 * @return flag
	 */
	public boolean isAutoFetch() {
		return autoFetch;
	}

	/** Tell the parser to fetch missing PDB files from the PDB FTP servers automatically.
	 *
	 * default is false. If true, new PDB files will be automatically stored in the Path and gzip compressed.
	 *
	 * @param autoFetch
	 */
	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}

	/* A flag to tell the parser to parse the Author's secondary structure assignment from the file
	 * default is set to false, i.e. do NOT parse.
	 * @param parseSecStruc
	 */
	public boolean isParseSecStruc() {
		return parseSecStruc;
	}


	/*  A flag to tell the parser to parse the Author's secondary structure assignment from the file
	 *
	 */

	public void setParseSecStruc(boolean parseSecStruc) {
		this.parseSecStruc = parseSecStruc;
	}


	/** directory where to find PDB files */
	public void setPath(String p){
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,p);
		ChemCompGroupFactory.setPath(p);
		path = p ;
	}

	/**
	 * Returns the path value.
	 * @return a String representing the path value
	 * @see #setPath
	 *
	 */
	public String getPath() {
		return path ;
	}

	/** define supported file extensions
	 * compressed extensions .Z,.gz do not need to be specified
	 * they are dealt with automatically.

	 */
	public void addExtension(String s){
		//System.out.println("add Extension "+s);
		extensions.add(s);
	}

	/** clear the supported file extensions
	 *
	 */
	public void clearExtensions(){
		extensions.clear();
	}

	/** Flag that defines if the PDB directory is containing all PDB files or is split into sub dirs (like the FTP site).
	 *  
	 * @return boolean. default is false (all files in one directory)
	 */
	public boolean isPdbDirectorySplit() {
		return pdbDirectorySplit;
	}

	/** Flag that defines if the PDB directory is containing all PDB files or is split into sub dirs (like the FTP site).
	 *  
	 * @param pdbDirectorySplit boolean. If set to false all files are in one directory.
	 */
	public void setPdbDirectorySplit(boolean pdbDirectorySplit) {
		this.pdbDirectorySplit = pdbDirectorySplit;
	}


	/** try to find the file in the filesystem and return a filestream in order to parse it
	 * rules how to find file
	 * - first check: if file is in path specified by PDBpath
	 * - secnd check: if not found check in PDBpath/xy/ where xy is second and third char of PDBcode.
	 * if autoFetch is set it will try to download missing PDB files automatically.
	 */

	private InputStream getInputStream(String pdbId)
	throws IOException
	{

		if ( pdbId.length() < 4)
			throw new IOException("the provided ID does not look like a PDB ID : " + pdbId);


		InputStream inputStream =null;

		String pdbFile = null ;
		File f = null ;

		// this are the possible PDB file names...
		String fpath ;
		String ppath ;

		if ( pdbDirectorySplit){
			// pdb files are split into subdirectories based on their middle position...
			String middle = pdbId.substring(1,3).toLowerCase();
			fpath = path+lineSplit + middle + lineSplit + pdbId;
			ppath = path +lineSplit +  middle + lineSplit + "pdb"+pdbId;
		} else {
			fpath = path+lineSplit + pdbId;
			ppath = path +lineSplit + "pdb"+pdbId;
		}

		String[] paths = new String[]{fpath,ppath};

		for ( int p=0;p<paths.length;p++ ){
			String testpath = paths[p];
			//System.out.println(testpath);
			for (int i=0 ; i<extensions.size();i++){
				String ex = (String)extensions.get(i) ;
				//System.out.println("PDBFileReader testing: "+testpath+ex);
				f = new File(testpath+ex) ;

				if ( f.exists()) {
					//System.out.println("found!");
					pdbFile = testpath+ex ;

					InputStreamProvider isp = new InputStreamProvider();

					inputStream = isp.getInputStream(pdbFile);
					break;
				}

				if ( pdbFile != null) break;
			}
		}

		if ( pdbFile == null ) {
			if ( autoFetch)
				return downloadAndGetInputStream(pdbId);

			String message = "no structure with PDB code " + pdbId + " found!" ;
			throw new IOException (message);
		}

		return inputStream ;
	}


	private File downloadPDB(String pdbId){

		if ((path == null) || (path.equals(""))){
			// accessing temp. OS directory:         
			String property = "java.io.tmpdir";

			String tempdir = System.getProperty(property);

			if ( !(tempdir.endsWith(lineSplit) ) )
				tempdir = tempdir + lineSplit;
			System.err.println("you did not set the path in PDBFileReader, don;t know where to write the downloaded file to");
			System.err.println("assuming default location is temp directory: " + tempdir);
			path = tempdir;
		}


		File tempFile ;

		if ( pdbDirectorySplit) {
			String middle = pdbId.substring(1,3).toLowerCase();
			String dir = path+lineSplit+middle;
			File directoryCheck = new File (dir);
			if ( ! directoryCheck.exists()){
				directoryCheck.mkdir();
			}

			tempFile =new File(dir+lineSplit+"pdb"+ pdbId.toLowerCase()+".ent.gz");

		} else {


			tempFile = new File(path+lineSplit+"pdb"+pdbId.toLowerCase()+".ent.gz");
		}


		String serverName = System.getProperty(PDB_FILE_SERVER_PROPERTY);

		if ( serverName == null)
			serverName = DEFAULT_PDB_FILE_SERVER;

		String ftp = String.format("ftp://%s/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz", serverName,pdbId.toLowerCase());

		System.out.println("Fetching " + ftp);

		// prepare destination
		System.out.println("writing to " + tempFile);

		try {
			URL url = new URL(ftp);

			InputStream uStream = url.openStream();
			InputStream conn = new GZIPInputStream(uStream);


			FileOutputStream outPut = new FileOutputStream(tempFile);
			GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);

			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(conn));
			String line;
			while ((line = fileBuffer.readLine()) != null) {
				pw.println(line);
			}
			pw.flush();
			pw.close();

			outPut.close();
			conn.close();
			uStream.close();

		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return tempFile;
	}

	private InputStream downloadAndGetInputStream(String pdbId)
	throws IOException{
		//PDBURLReader reader = new PDBURLReader();
		//Structure s = reader.getStructureById(pdbId);
		File tmp = downloadPDB(pdbId);
		if ( tmp != null ) {
			InputStreamProvider prov = new InputStreamProvider();
			return prov.getInputStream(tmp);


		} else {
			throw new IOException("could not find PDB " + pdbId + " in file system and also could not download");
		}


	}


	/** load a structure from local file system and return a PDBStructure object

	 * @param pdbId  a String specifying the id value (PDB code)
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructureById(String pdbId)
	throws IOException
	{


		InputStream inStream = getInputStream(pdbId);

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setParseSecStruc(parseSecStruc);
		pdbpars.setAlignSeqRes(alignSeqRes);
		pdbpars.setHeaderOnly(headerOnly);
		pdbpars.setParseCAOnly(parseCAOnly);
		pdbpars.setLoadChemCompInfo(loadChemCompInfo);


		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;
	}

	/** opens filename, parses it and returns
	 * aStructure object .
	 * @param filename  a String
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(String filename)
	throws IOException
	{
		File f = new File(filename);
		return getStructure(f);

	}

	/** opens filename, parses it and returns a Structure object
	 *
	 * @param filename a File object
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(File filename) throws IOException {

		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(filename);

		return getStructure(inStream);
	}

	private Structure getStructure(InputStream inStream) throws IOException{

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setParseSecStruc(parseSecStruc);
		pdbpars.setAlignSeqRes(alignSeqRes);
		pdbpars.setParseCAOnly(parseCAOnly);
		pdbpars.setHeaderOnly(headerOnly);

		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;

	}

	public Structure getStructure(URL u) throws IOException{

		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(u);

		return getStructure(inStream);

	}


	public boolean isHeaderOnly() {
		return headerOnly;
	}


	public void setHeaderOnly(boolean flag) {
		headerOnly = flag;

	}

	public boolean isLoadChemCompInfo() {
		return loadChemCompInfo;
	}
	public void setLoadChemCompInfo(boolean loadChemCompInfo) {
		if ( loadChemCompInfo)
			System.setProperty(LOAD_CHEM_COMP_PROPERTY, "true");
		this.loadChemCompInfo = loadChemCompInfo;
	}
}
