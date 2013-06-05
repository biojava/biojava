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


import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;


import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Compound;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.PDBStatus;
import org.biojava.bio.structure.PDBStatus.Status;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.bio.structure.io.mmcif.ReducedChemCompProvider;
import org.biojava.bio.structure.io.util.FileDownloadUtils;
import org.biojava3.core.util.InputStreamProvider;

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


	String path;
	List<String> extensions;

	boolean autoFetch;

	private boolean fetchCurrent;

	private boolean fetchFileEvenIfObsolete;

	boolean pdbDirectorySplit;

	private int bioAssemblyId = 0; // number > 0 indicates the id of the biological assembly
	private boolean bioAssemblyFallback = true; // use regular PDB file as the biological assembly (i.e. for NMR structures)
	// in case no biological assembly file is available.
	private boolean loadedBioAssembly = false;

	public static final String lineSplit = System.getProperty("file.separator");

	public static final String DEFAULT_PDB_FILE_SERVER = "ftp.wwpdb.org";
	public static final String PDB_FILE_SERVER_PROPERTY = "PDB.FILE.SERVER";

	private static final String CURRENT_FILES_PATH  = "/pub/pdb/data/structures/divided/pdb/";
	private static final String OBSOLETE_FILES_PATH = "/pub/pdb/data/structures/obsolete/pdb/";
	private static final String BIO_ASSEMBLY_FILES_PATH  = "/pub/pdb/data/biounit/coordinates/divided/";
	private static final String LOCAL_BIO_ASSEMBLY_DIRECTORY = "bio_assembly";
	FileParsingParameters params ;

	public static final long lastRemediationDate ;

	static {

		SimpleDateFormat formatter = new SimpleDateFormat("yyyy/MM/dd", Locale.US);

		long t = 0;
		try {
			Date d = formatter.parse("2011/07/12");
			t = d.getTime();
		} catch (Exception e){
			e.printStackTrace();					
		}
		lastRemediationDate = t;
	}

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
		} catch (Exception e) {
			e.printStackTrace();
		}
	}





	public PDBFileReader() {
		extensions    = new ArrayList<String>();

		extensions.add(".ent");
		extensions.add(".pdb");
		extensions.add(".ent.gz");
		extensions.add(".pdb.gz");
		extensions.add(".ent.Z");
		extensions.add(".pdb.Z");

		autoFetch     = false;		
		pdbDirectorySplit = false;

		params = new FileParsingParameters();

		//checkPath();
	}

	/** Check the directory that contains the local cache of PDB and chemical component definiton files.
	 * By default, at startup we fall back to java.io.tmpdir
	 */
	private  void checkPath(){

		if ((path == null) || (path.equals("")) || path.equals("null")) {

			String syspath = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);

			if ((syspath != null) && (! syspath.equals("")) && (! syspath.equals("null"))){

				path = syspath;
				return;
			}

			// accessing temp. OS directory:         
			String property = "java.io.tmpdir";

			String tempdir = System.getProperty(property);

			if ( !(tempdir.endsWith(lineSplit) ) )
				tempdir = tempdir + lineSplit;

			// 
			System.err.println("you did not set the path in PDBFileReader, don't know where to write the downloaded file to");
			System.err.println("assuming default location is temp directory: " + tempdir);
			path = tempdir;
		}
	}






	/** directory where to find PDB files */
	public void setPath(String p){
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,p);

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

	/** 
	 * Sets the ID of the biological assembly to be loaded. By default, the bioAssemblyId=0, which corresponds to the
	 * original PDB file. If an ID > 0 is specified, the corresponding biological assembly file will be loaded. Note, the
	 * number of available biological unit files varies. Many entries don't have a biological unit specified (i.e. NMR structures),
	 * many entries have only one biological unit (bioAssemblyId=1), and a few structures have multiple biological assemblies.
	 *  
	 * @param bioAssemblyId ID of the biological assembly to be loaded
	 * @author Peter Rose
	 * @since 3.2
	 */
	public void setBioAssemblyId(int bioAssemblyId) {
		this.bioAssemblyId = bioAssemblyId;
	}

	/**
	 * Set bioAssemblyFallback to true, to download the original PDB file in cases that a biological assembly file is not available.
	 * 
	 * @param bioAssemblyFallback if true, tries reading original PDB file in case the biological assembly file is not available
	 * @author Peter Rose
	 * @since 3.2
	 */
	public void setBioAssemblyFallback(boolean bioAssemblyFallback) {
		this.bioAssemblyFallback = bioAssemblyFallback;
	}

	/** try to find the file in the filesystem and return a filestream in order to parse it
	 * rules how to find file
	 * - first check: if file is in path specified by PDBpath
	 * - second check: if not found check in PDBpath/xy/ where xy is second and third char of PDBcode.
	 * if autoFetch is set it will try to download missing PDB files automatically.
	 */

	public InputStream getInputStream(String pdbId)
			throws IOException {

		if ( pdbId.length() < 4)
			throw new IOException("the provided ID does not look like a PDB ID : " + pdbId);

		checkPath();

		InputStream inputStream =null;

		String pdbFile = getLocalPDBFilePath(pdbId);

		if ( pdbFile != null) {
			InputStreamProvider isp = new InputStreamProvider();

			try {
				inputStream = isp.getInputStream(pdbFile);
				return inputStream;
			} catch (Exception e){
				e.printStackTrace();
				// something is wrong with the file!
				// it probably should be downloaded again...
				pdbFile = null;
			}
		}
		if ( pdbFile == null ) {
			if (autoFetch){//from here we try our online search
				if(fetchCurrent && !fetchFileEvenIfObsolete) {
					String current = PDBStatus.getCurrent(pdbId);

					if(current == null) {
						// either an error or there is not current entry
						current = pdbId;
					}
					return downloadAndGetInputStream(current, CURRENT_FILES_PATH);
				} else if(fetchFileEvenIfObsolete && PDBStatus.getStatus(pdbId) == Status.OBSOLETE) {
					return downloadAndGetInputStream(pdbId, OBSOLETE_FILES_PATH);
				} else {
					return downloadAndGetInputStream(pdbId, CURRENT_FILES_PATH);
				}
			}else {
				String message = "no structure with PDB code " + pdbId + " found!" ;
				throw new IOException (message);
			}
		}
		return null ;
	}

	/** Returns null if local PDB file does not exist or should be downloaded again...
	 * 
	 * @param pdbId
	 * @return
	 */
	private String getLocalPDBFilePath(String pdbId) {

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

					pdbFile = testpath+ex ;

					// we have found the file locally

					if ( params.isUpdateRemediatedFiles()){
						long lastModified = f.lastModified();

						if (lastModified < lastRemediationDate) {
							// the file is too old, replace with newer version
							System.out.println("replacing file " + pdbFile +" with latest remediated file from PDB.");
							pdbFile = null;

							return null;
						}
					}

					return pdbFile;
				}

				if ( pdbFile != null) break;
			}
		}
		return null;
	}





	/**
	 * Returns an input stream for a biological assembly file based on the passed in pdbId and 
	 * the biological assembly id {@link #setBioAssemblyId(int)}. Files are cached in a local directory.
	 * @param pdbId
	 * @return InputStream
	 * @throws IOException
	 * @author Peter Rose
	 * @since 3.2
	 */
	public InputStream getInputStreamBioAssembly(String pdbId)
			throws IOException
			{
		loadedBioAssembly = true;
		InputStream inputStream = null;

		if ( pdbId.length() < 4)
			throw new IOException("the provided ID does not look like a PDB ID : " + pdbId);

		checkPath();

		// make sure path follows unix convention
		String uPath = FileDownloadUtils.toUnixPath(path);

		// create local subdirectory for biological assembly files if it doesn't exist
		String dir = uPath + LOCAL_BIO_ASSEMBLY_DIRECTORY;

		File tmp = new File(dir);
		if ( ! tmp.exists()){
			tmp.mkdir();
		}

		String fpath ;

		if (pdbDirectorySplit){
			// pdb files are split into subdirectories based on their middle position...
			String middle = pdbId.substring(1,3).toLowerCase();
			fpath = dir + "/" + middle + "/";
		} else {
			fpath = dir + "/";
		}  

		String fileName = fpath + getBiologicalAsssemblyFileName(pdbId.toLowerCase(), bioAssemblyId);
		File f = new File(fileName) ;

		// check if bio assembly file exists in local cache
		if ( f.exists()) {
			InputStreamProvider isp = new InputStreamProvider();

			try {
				inputStream = isp.getInputStream(fileName);
				return inputStream;
			} catch (Exception e){
				e.printStackTrace();
				// something is wrong with the file!
				// it probably should be downloaded again...
			}
		} else if (bioAssemblyFallback) {
			if (pdbDirectorySplit){
				// pdb files are split into subdirectories based on their middle position...
				String middle = pdbId.substring(1,3).toLowerCase();
				fpath = uPath + middle + "/";
			} else {
				fpath = uPath;
			}  
			fileName = fpath + "pdb" + pdbId + ".ent.gz";
			f = new File(fileName);

			if (f.exists()) {

				try {
					InputStreamProvider isp = new InputStreamProvider();
					inputStream = isp.getInputStream(fileName);
					System.out.println("Loaded original PDB file as a fallback." + fileName);
					loadedBioAssembly = false;
					return inputStream;
				} catch (Exception e) {
					// now try autofetch next
				}
			}
		}

		if (autoFetch){//from here we try our online search
			if(fetchCurrent && !fetchFileEvenIfObsolete) {
				String current = PDBStatus.getCurrent(pdbId);

				if(current == null) {
					// either an error or there is not current entry
					current = pdbId;
				}
				inputStream  = downloadAndGetInputStreamBioAssembly(current, BIO_ASSEMBLY_FILES_PATH);
				if (inputStream != null) {
					return inputStream;
				}
			} else if(fetchFileEvenIfObsolete && PDBStatus.getStatus(pdbId) == Status.OBSOLETE) {
				String message = "No biological assembly with PDB code " + pdbId + " found!" ;
				throw new IOException (message);
			} else {
				inputStream = downloadAndGetInputStreamBioAssembly(pdbId, BIO_ASSEMBLY_FILES_PATH);
				if (inputStream != null) {
					return inputStream;
				}
			}
		} 

		// if biological assembly file cannot be found, and bioAssemblyFallBack is true,
		// get the original PDB file as a fall back (i.e. for NMR structures, where the
		// PDB file represents the biological assembly).
		if (bioAssemblyFallback) {	
			inputStream = getInputStream(pdbId);		
			if (inputStream != null) {
				System.out.println("Biological assembly file for PDB ID: " + pdbId+  " is not available. " +
						"Loaded original PDB file as a fallback from local cache.");
				return getInputStream(pdbId);
			}
		}

		return null;
			}

	private  File downloadPDB(String pdbId, String pathOnServer){

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



		File realFile = null;

		pdbId = pdbId.toLowerCase();
		String middle = pdbId.substring(1,3);

		if ( pdbDirectorySplit) {
			String dir = path+lineSplit+middle;
			File directoryCheck = new File (dir);
			if ( ! directoryCheck.exists()){
				directoryCheck.mkdir();
			}

			realFile =new File(dir+lineSplit+"pdb"+ pdbId+".ent.gz");
		} else {

			realFile = new File(path+lineSplit+"pdb"+pdbId+".ent.gz");
		}


		String serverName = System.getProperty(PDB_FILE_SERVER_PROPERTY);

		if ( serverName == null)
			serverName = DEFAULT_PDB_FILE_SERVER;

		String ftp = String.format("ftp://%s%s%s/pdb%s.ent.gz", serverName,pathOnServer,middle, pdbId);

		//System.out.println("Fetching " + ftp);


		try {
			URL url = new URL(ftp);

			FileDownloadUtils.downloadGzipCompressedFile(url, realFile);

		} catch (Exception e){
			System.err.println("Problem while downloading PDB ID " + pdbId + " from FTP server." );
			System.err.println("Failed URL: " + ftp);
			e.printStackTrace();
			return null;
		}


		return realFile;
	}




	/**
	 * Downloads a biological assembly file. If this file cannot be found, it will download the original 
	 * PDB file (i.e. for NMR structures), if bioAssemblyFallback has been set to true;
	 * @param pdbId
	 * @param pathOnServer
	 * @return File of biological assembly
	 * @author Peter Rose
	 * @since 3.2
	 */
	private  File downloadPDBBiologicalAssembly(String pdbId, String pathOnServer){	
		loadedBioAssembly = true;

		checkPath();

		String fileName = getBiologicalAsssemblyFileName(pdbId, bioAssemblyId);

		pdbId = pdbId.toLowerCase();
		String middle = pdbId.substring(1,3);

		String serverName = System.getProperty(PDB_FILE_SERVER_PROPERTY);

		if ( serverName == null)
			serverName = DEFAULT_PDB_FILE_SERVER;

		String ftp = String.format("ftp://%s%s%s/%s", serverName,pathOnServer,middle,fileName);

		System.out.println("Fetching " + ftp);

		URL url = null;
		try {
			url = new URL(ftp);
		} catch (MalformedURLException e1) {
			System.err.println("Problem while downloading Biological Assembly " + pdbId + " from FTP server." );
			e1.printStackTrace();
			return null;
		}

		// check if the file exists on the FTP site. If biological assembly file does not exist
		// and the fallback has been set, get the original PDB file instead (i.e. for NMR structures).

		File f = null;

		try {
			f = downloadFileIfAvailable(url, pdbId,fileName);
		} catch (IOException ioe){
			// ignore here because it prob. means the file does not exist...
		}

		if ( f == null) {

			if (bioAssemblyFallback) {

				System.out.println("Biological unit file for PDB ID: " + pdbId+  " is not available. " +
						"Downloading original PDB file as a fallback.");

				loadedBioAssembly = false;

				String fallBackPDBF = getLocalPDBFilePath(pdbId);

				if ( fallBackPDBF != null)
					return new File(fallBackPDBF);

				return downloadPDB(pdbId, CURRENT_FILES_PATH);
			} 
			return null;
		}

		return f;


	}





	private File downloadFileIfAvailable(URL url, String pdbId, String fileName) throws IOException {

		String middle = pdbId.substring(1,3);
		String uPath = FileDownloadUtils.toUnixPath(path);
		File tempFile = null;
		if (pdbDirectorySplit) {
			String dir = uPath + LOCAL_BIO_ASSEMBLY_DIRECTORY + "/" + middle;
			File directoryCheck = new File (dir);
			if ( ! directoryCheck.exists()){
				directoryCheck.mkdir();
			}
			tempFile =new File(dir + "/" + fileName);
		} else {
			tempFile = new File(path + LOCAL_BIO_ASSEMBLY_DIRECTORY + "/" + fileName);
		}

		return FileDownloadUtils.downloadFileIfAvailable(url, tempFile);
	}





	private  InputStream downloadAndGetInputStream(String pdbId, String pathOnServer)
			throws IOException{

		File tmp = downloadPDB(pdbId, pathOnServer);

		if (tmp != null) {
			InputStreamProvider prov = new InputStreamProvider();
			InputStream is = prov.getInputStream(tmp);
			return is;
		} else {
			throw new IOException("could not find PDB " + pdbId + " in file system and also could not download");
		}


	}

	/**
	 * Downloads biological assembly file to local cache and provides input stream to cached file.
	 * @param pdbId
	 * @param pathOnServer
	 * @return inputStream to cached file
	 * @throws IOException
	 * @author Peter Rose
	 * @since 3.2
	 */
	private  InputStream downloadAndGetInputStreamBioAssembly(String pdbId, String pathOnServer)
			throws IOException{

		File tmp = downloadPDBBiologicalAssembly(pdbId, pathOnServer);

		if (tmp != null) {
			InputStreamProvider prov = new InputStreamProvider();
			InputStream is = prov.getInputStream(tmp);
			return is;
		} else {
			throw new IOException("Could not find Biological Assembly " + pdbId + " in file system and also could not download");
		}


	}

	public boolean checkFileExists(String pdbId){
		String path =  getLocalPDBFilePath(pdbId);
		if ( path != null)
			return true;
		return false;

	}

	public void downloadPDB(String pdbId){
		//don't overwrite existing files
		if ( checkFileExists(pdbId))
			return;

		if (autoFetch){//from here we try our online search
			if(fetchCurrent && !fetchFileEvenIfObsolete) {
				String current = PDBStatus.getCurrent(pdbId);

				if(current == null) {
					// either an error or there is not current entry
					current = pdbId;
				}
				downloadPDB(current, CURRENT_FILES_PATH);
			} else if(fetchFileEvenIfObsolete && PDBStatus.getStatus(pdbId) == Status.OBSOLETE) {
				downloadPDB(pdbId, OBSOLETE_FILES_PATH);
			} else {
				downloadPDB(pdbId, CURRENT_FILES_PATH);
			}
		}

	}


	/** load a structure from local file system and return a PDBStructure object

	 * @param pdbId  a String specifying the id value (PDB code)
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public  Structure getStructureById(String pdbId)throws IOException	{

		InputStream inStream = null;
		if (bioAssemblyId == 0) {
			inStream = getInputStream(pdbId);
		} else {
			//System.out.println("loading bioassembly " + bioAssemblyId);
			inStream = getInputStreamBioAssembly(pdbId);
		}

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setFileParsingParameters(params);

		Structure struc = pdbpars.parsePDBFile(inStream) ;
		struc.setBiologicalAssembly(loadedBioAssembly);
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
		pdbpars.setFileParsingParameters(params);

		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;

	}

	public Structure getStructure(URL u) throws IOException{
		InputStreamProvider isp = new InputStreamProvider();
		InputStream inStream = isp.getInputStream(u);
		return getStructure(inStream);
	}



	public void setFileParsingParameters(FileParsingParameters params){
		this.params= params;
		if ( ! params.isLoadChemCompInfo()) {
			ChemCompGroupFactory.setChemCompProvider(new ReducedChemCompProvider());
		}
	}

	public FileParsingParameters getFileParsingParameters(){
		return params;
	}


	public boolean isAutoFetch(){
		return autoFetch;
	}


	public void setAutoFetch(boolean autoFetch){
		this.autoFetch = autoFetch;
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchFileEvenIfObsolete the fetchFileEvenIfObsolete to set
	 */
	public void setFetchFileEvenIfObsolete(boolean fetchFileEvenIfObsolete) {
		this.fetchFileEvenIfObsolete = fetchFileEvenIfObsolete;
	}

	/**forces the reader to fetch the file if its status is OBSOLETE.
	 * This feature has a higher priority than {@link #setFetchCurrent(boolean)}. <br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchFileEvenIfObsolete
	 * @author Amr AL-Hossary
	 * @see #fetchCurrent
	 * @since 3.0.2
	 */
	public boolean isFetchFileEvenIfObsolete() {
		return fetchFileEvenIfObsolete;
	}


	/**if enabled, the reader searches for the newest possible PDB ID, if not present in he local installation.
	 * The {@link #setFetchFileEvenIfObsolete(boolean)} function has a higher priority than this function. <br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchCurrent the fetchCurrent to set
	 * @author Amr AL-Hossary
	 * @see #setFetchFileEvenIfObsolete(boolean)
	 * @since 3.0.2
	 */
	public void setFetchCurrent(boolean fetchNewestCurrent) {
		this.fetchCurrent = fetchNewestCurrent;
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchCurrent
	 */
	public boolean isFetchCurrent() {
		return fetchCurrent;
	}

	/**
	 * Returns the file name of a PDB biological unit file based on the pdbId and biologicalAssemblyID.
	 * 
	 * @param pdbId the protein data bank ID
	 * @param biologicalAssemblyId the ID of the biological assembly
	 * @return file name of PDB biological assembly file
	 * @author Peter Rose
	 * @since 3.2
	 */
	private String getBiologicalAsssemblyFileName(String pdbId, int biologicalAssemblyId) {
		return pdbId + ".pdb" + biologicalAssemblyId + ".gz";
	}




}
