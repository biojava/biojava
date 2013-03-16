package org.biojava.bio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.GZIPOutputStream;


import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.core.util.InputStreamProvider;



/** This provider of chemical components can download and cache chemical component definition files from the RCSB PDB web site.
 *  It is the default way to access these definitions.
 *  If this provider is called he first time, it will download and install all chemical 
 *  component definitions in a local directory. 
 *  Once the definition files have been installed, it has quick startup time and low memory requirements. 
 *
 *  An alternative provider, that keeps all definitions in memory is the {@link AllChemCompProvider}. Another provider, that
 *  does not require any network access, but only can support a limited set of chemical component definitions, is the {@link ReducedChemCompProvider}. 
 *  
 * 
 * @author Andreas Prlic
 *
 */
public class DownloadChemCompProvider implements ChemCompProvider {

	static String path; 
	private static final String FILE_SEPARATOR = System.getProperty("file.separator");
	static final String NEWLINE = System.getProperty("line.separator");

	public static String CHEM_COMP_CACHE_DIRECTORY = "chemcomp";

	private static String serverLocation = "http://www.rcsb.org/pdb/files/ligand/";

	// flags to make sure there is only one thread running that is loading the dictionary
	static AtomicBoolean loading = new AtomicBoolean(false);

	static final List<String> protectedIDs = new ArrayList<String> ();
	static {
		protectedIDs.add("CON");
		protectedIDs.add("PRN");
		protectedIDs.add("AUX");
		protectedIDs.add("NUL");
	}
	
	/** by default we will download only some of the files. User has to request that all files should be downloaded...
	 * 
	 */
	boolean downloadAll = false;

	public DownloadChemCompProvider(){
		//System.out.println("USING DOWNLOAD CHEM COMP PROVIDER");		
	}

	/** checks if the chemical components already have been installed into the PDB directory.
	 *  If not, will download the chemical components definitions file and split it up into small
	 *  subfiles.
	 */
	public void checkDoFirstInstall(){

		if ( ! downloadAll ) {
			return;
		}
		
		if ( path == null)
			path =  System.getProperty(AbstractUserArgumentProcessor.CACHE_DIR);
		if (path == null || path.equals(""))
			path = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);
		
		String filename = path + 	
		DownloadChemCompProvider.CHEM_COMP_CACHE_DIRECTORY +
		FILE_SEPARATOR + 
		"components.cif.gz";

		File f = new File(filename);

		if ( ! f.exists()) {

			downloadAllDefinitions();

		} else {
			// file exists.. did it get extracted?
			String directoryName = path + 	
			DownloadChemCompProvider.CHEM_COMP_CACHE_DIRECTORY +
			FILE_SEPARATOR;

			File dir = new File(directoryName);

			FilenameFilter filter =new FilenameFilter() {

				public boolean accept(File dir, String file) {
					return file.endsWith(".cif.gz");
				}
			};
			String[] files = dir.list(filter);
			if ( files.length < 500) {
				// not all did get unpacked
				split();
			}
		}
	}

	private void split(){

		System.out.println("Installing individual chem comp files ...");
		
		String filename = path + 
		DownloadChemCompProvider.CHEM_COMP_CACHE_DIRECTORY + FILE_SEPARATOR +
		"components.cif.gz";

		int counter = 0;
		InputStreamProvider prov = new InputStreamProvider();
		try {
			InputStream inStream = prov.getInputStream(filename);

			BufferedReader buf = new BufferedReader (new InputStreamReader (inStream));

			String line = null;
			line = buf.readLine ();
			StringWriter writer = new StringWriter();

			String currentID = null;
			while (line != null){

				if ( line.startsWith("data_")) {
					// a new record found!

					if ( currentID != null) {
						writeID(writer, currentID);
						counter++;
					}

					currentID = line.substring(5);
					writer = new StringWriter();
				}

				writer.append(line);
				writer.append(NEWLINE);

				line = buf.readLine ();
			}

			// write the last record...
			writeID(writer,currentID);
			counter++;
		} catch (Exception e){
			e.printStackTrace();
		}
		System.out.println("created " + counter + " chemical component files.");
	}

	private void writeID(StringWriter writer, String currentID) throws IOException{

		//System.out.println("writing ID " + currentID);

		
		
		String localName = DownloadChemCompProvider.getLocalFileName(currentID);

		FileOutputStream outPut = new FileOutputStream(localName);

		GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);

		PrintWriter pw = new PrintWriter(gzOutPut);

		pw.print(writer.toString());
		writer.close();
		pw.flush();
		pw.close();

		outPut.close();

	}

	/** Loads the definitions for this {@link ChemComp} from a local file and instantiates a new object.
	 * 
	 * @param recordName the ID of the {@link ChemComp}
	 * @return a new {@link ChemComp} definition.
	 */
	public  ChemComp getChemComp(String recordName) {

		checkPath();

		// make sure we work with upper case records		
		recordName = recordName.toUpperCase().trim();

		if ( recordName.equals("?")){
			return null;
		}
		try {
			if ( ! fileExists(recordName)) {
				// check if we should install all components				
				checkDoFirstInstall();
			}
			if ( ! fileExists(recordName)) {
				// we previously have installed already the definitions,
				// just do an incrememntal update
				downloadChemCompRecord(recordName);
			}

			String filename = getLocalFileName(recordName);

			System.out.println("reading " + filename);
			InputStreamProvider isp = new InputStreamProvider();

			InputStream inStream = isp.getInputStream(filename);

			MMcifParser parser = new SimpleMMcifParser();

			ChemCompConsumer consumer = new ChemCompConsumer();

			// The Consumer builds up the BioJava - structure object.
			// you could also hook in your own and build up you own data model.
			parser.addMMcifConsumer(consumer);

			parser.parse(new BufferedReader(new InputStreamReader(inStream)));

			ChemicalComponentDictionary dict = consumer.getDictionary();

			ChemComp chemComp = dict.getChemComp(recordName);

			return chemComp;

		} catch (Exception e) {

			e.printStackTrace();

		}
		return null;

	}

	private static void checkPath(){

		if ((path == null) || (path.equals("")) || path.equals("null")) {

			String syspath = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);

			if ((syspath != null) && (! syspath.equals("")) && (! syspath.equals("null"))){

				path = syspath;
				return;
			}

			// accessing temp. OS directory:         
			String property = "java.io.tmpdir";

			String tempdir = System.getProperty(property);

			if ( !(tempdir.endsWith(FILE_SEPARATOR) ) )
				tempdir = tempdir + FILE_SEPARATOR;

			System.err.println("you did not set the path in PDBFileReader, don't know where to write the downloaded file to");
			System.err.println("assuming default location is temp directory: " + tempdir);
			path = tempdir;
		}
	}

	/** Returns the file name that contains the definition for this {@link ChemComp}
	 *  
	 * @param recordName the ID of the {@link ChemComp}
	 * @return full path to the file
	 */
	public static String getLocalFileName(String recordName){

		if ( protectedIDs.contains(recordName)){
			recordName = "_" + recordName;
		}
		
		String dir = path + CHEM_COMP_CACHE_DIRECTORY + FILE_SEPARATOR;

		File f = new File(dir);
		if (! f.exists()){
			System.out.println("creating directory " + f);
			f.mkdir();
		}

		String fileName = path + CHEM_COMP_CACHE_DIRECTORY + FILE_SEPARATOR + recordName + ".cif.gz";

		return fileName;
	}

	private static  boolean fileExists(String recordName){

		String fileName = getLocalFileName(recordName);

		File f = new File(fileName);

		return f.exists();

	}

	private static void downloadChemCompRecord(String recordName) {
		String path = System.getProperty(AbstractUserArgumentProcessor.PDB_DIR);
		setPath(path);


		String localName = getLocalFileName(recordName);

		String u = serverLocation + recordName + ".cif";

		System.out.println("downloading " + u);

		try {

			URL url = new URL(u);

			HttpURLConnection uconn = HTTPConnectionTools.openHttpURLConnection(url);

			InputStream conn = uconn.getInputStream();

			FileOutputStream outPut = new FileOutputStream(localName);

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


		} catch (Exception e){
			e.printStackTrace();
		}


	}


	/** making sure we use the same path for the PDB installation as is used by the PdbFileReader
	 * 
	 * @param p path to PDB files.
	 */
	public static void setPath(String p) {
		if ( p == null)
			return;
		path = p;
		if ( ! path.endsWith(FILE_SEPARATOR))
			path += FILE_SEPARATOR;
		
		System.setProperty(AbstractUserArgumentProcessor.CACHE_DIR,path);
		
		

	}

	private void downloadAllDefinitions() {

		if ( loading.get()){
			System.out.println("Waiting for other thread to install chemical components...");
		}

		while ( loading.get() ) {

			// another thread is already downloading the components definitions
			// wait for the other thread to finish...

			try {
				// wait half a second

				Thread.sleep(500);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			System.out.println("Another thread installed the chemical components.");
			return;

		}

		loading.set(true);
		long timeS = System.currentTimeMillis();

		System.out.println("Performing first installation of chemical components.");
		System.out.println("Downloading components.cif.gz ...");

		AllChemCompProvider.checkPath();
		try {
			AllChemCompProvider.downloadFile();
		} catch (IOException e){
			e.printStackTrace();
		}
		
		split();
		long timeE = System.currentTimeMillis();		
		System.out.println("time to install chem comp dictionary: " + (timeE - timeS) / 1000 + " sec.");		
		loading.set(false);

	}

	/** By default this provider will download only some of the {@link ChemComp} files. 
	 * The user has to request that all files should be downloaded by setting this parameter to true.
	 * 
	 *  @return flag if the all components should be downloaded and installed at startup. (default: false)
	 */
	public boolean isDownloadAll() {
		return downloadAll;
	}

	/** By default this provider will download only some of the {@link ChemComp} files. 
	 * The user has to request that all files should be downloaded by setting this parameter to true.
	 * 
	 * @param  flag if the all components should be downloaded and installed at startup. (default: false)
	 */
	public void setDownloadAll(boolean downloadAll) {
		this.downloadAll = downloadAll;
	}





}
