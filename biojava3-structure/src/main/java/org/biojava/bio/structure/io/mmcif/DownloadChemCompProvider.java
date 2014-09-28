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

import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;



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

	private static final Logger logger = LoggerFactory.getLogger(DownloadChemCompProvider.class);
	
	private static File path; 
	//private static final String FILE_SEPARATOR = System.getProperty("file.separator");
	private static final String NEWLINE = System.getProperty("line.separator");

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
		logger.debug("Initialising DownloadChemCompProvider");	
		
		UserConfiguration config = new UserConfiguration();
		// TODO should this be getPdbFilePath() or getCacheFilePath()?, in AllChemCompProvider it's getCacheFilePath()
		path = new File(config.getPdbFilePath());
		
	}

	/** checks if the chemical components already have been installed into the PDB directory.
	 *  If not, will download the chemical components definitions file and split it up into small
	 *  subfiles.
	 */
	public void checkDoFirstInstall(){ 

		if ( ! downloadAll ) {
			return;
		}
		
		
		// this makes sure there is a file separator between every component,
		// if path has a trailing file separator or not, it will work for both cases
		File dir = new File(path, CHEM_COMP_CACHE_DIRECTORY); 
		File f = new File(dir, "components.cif.gz");
		
		if ( ! f.exists()) {

			downloadAllDefinitions();

		} else {
			// file exists.. did it get extracted?
			
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

		logger.info("Installing individual chem comp files ...");
		
		File dir = new File(path, CHEM_COMP_CACHE_DIRECTORY);
		File f = new File(dir, "components.cif.gz");
		

		int counter = 0;
		InputStreamProvider prov = new InputStreamProvider();
		try {
			InputStream inStream = prov.getInputStream(f);

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
		} catch (IOException e){
			e.printStackTrace();
		}
		logger.info("Created " + counter + " chemical component files.");
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

//			System.out.println("reading " + filename);
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

		} catch (IOException e) {

			e.printStackTrace();

		}
		return null;

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
		
		
		File f = new File(path, CHEM_COMP_CACHE_DIRECTORY);
		if (! f.exists()){
			logger.info("Creating directory " + f);
			
			boolean success = f.mkdir();
			// we've checked in initPath that path is writable, so there's no need to check if it succeeds
			// in the unlikely case that in the meantime it isn't writable at least we log an error 
			if (!success) logger.error("Directory {} could not be created",f);
			
		}

		File theFile = new File(f,recordName + ".cif.gz");
		
		return theFile.toString();
	}

	private static  boolean fileExists(String recordName){

		String fileName = getLocalFileName(recordName);

		File f = new File(fileName);

		return f.exists();

	}

	private static void downloadChemCompRecord(String recordName) {
		
		String localName = getLocalFileName(recordName);

		String u = serverLocation + recordName + ".cif";

//		System.out.println("downloading " + u);

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


		} catch (IOException e){
			e.printStackTrace();
		}


	}

	private void downloadAllDefinitions() {

		if ( loading.get()){
			logger.info("Waiting for other thread to install chemical components...");
		}

		while ( loading.get() ) {

			// another thread is already downloading the components definitions
			// wait for the other thread to finish...

			try {
				// wait half a second

				Thread.sleep(500);
			} catch (InterruptedException e) {
				//e.printStackTrace();
				logger.error("Thread interrupted "+e.getMessage());
			}

			logger.info("Another thread installed the chemical components.");
			return;

		}

		loading.set(true);
		long timeS = System.currentTimeMillis();

		logger.info("Performing first installation of chemical components.");
		logger.info("Downloading components.cif.gz ...");

		AllChemCompProvider.checkPath();
		try {
			AllChemCompProvider.downloadFile();
		} catch (IOException e){
			e.printStackTrace();
		}
		
		split();
		long timeE = System.currentTimeMillis();		
		logger.info("time to install chem comp dictionary: " + (timeE - timeS) / 1000 + " sec.");		
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
