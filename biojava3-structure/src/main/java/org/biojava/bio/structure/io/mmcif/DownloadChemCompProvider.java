package org.biojava.bio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.zip.GZIPOutputStream;


import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.util.HTTPConnectionTools;
import org.biojava.bio.structure.io.mmcif.model.ChemComp;
import org.biojava3.core.util.InputStreamProvider;


/** This class can download and caches chemical component files from the RCSB PDB web site.  
 * 
 * @author Andreas Prlic
 *
 */
public class DownloadChemCompProvider implements ChemCompProvider {

	static String path; 
	private static final String lineSplit = System.getProperty("file.separator");

	private static String dirName = "chemcomp";

	private static String serverLocation = "http://www.rcsb.org/pdb/files/ligand/";

	
	public DownloadChemCompProvider(){
		//System.out.println("USING DOWNLOAD CHEM COMP PROVIDER");
	}
	
	public  ChemComp getChemComp(String recordName) {

		checkPath();
		
		// make sure we work with upper case records		
		recordName = recordName.toUpperCase().trim();

		if ( recordName.equals("?")){
			return null;
		}
		try {
			if ( ! fileExists(recordName)) {

				downloadChemCompRecord(recordName);
			}

			String filename = getLocalFileName(recordName);

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

			if ( !(tempdir.endsWith(lineSplit) ) )
				tempdir = tempdir + lineSplit;

			System.err.println("you did not set the path in PDBFileReader, don't know where to write the downloaded file to");
			System.err.println("assuming default location is temp directory: " + tempdir);
			path = tempdir;
		}
	}


	private static String getLocalFileName(String recordName){

		String dir = path + dirName + lineSplit;

		File f = new File(dir);
		if (! f.exists()){
			System.out.println("creating directory " + f);
			f.mkdir();
		}

		String fileName = path + dirName + lineSplit + recordName + ".cif.gz";

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
		path = p;
		if ( ! path.endsWith(lineSplit))
			path += lineSplit;

	}



	

}
