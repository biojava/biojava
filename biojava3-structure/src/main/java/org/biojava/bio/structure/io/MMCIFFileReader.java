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
 * created at Oct 18, 2008
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
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.StructureIOFile;
import org.biojava.bio.structure.io.mmcif.MMcifParser;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava3.core.util.InputStreamProvider;


/** How to parse an mmCif file:
 * <pre>
  public static void main(String[] args){
        String filename =  "/path/to/something.cif.gz" ;

        StructureIOFile reader = new MMCIFFileReader();

        try{
            Structure struc = reader.getStructure(filename);
            System.out.println(struc);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    </pre>
 *
 * @author Andreas Prlic
 * @since 1.7
 *
 */
public class MMCIFFileReader implements StructureIOFile {

	String path;
	List<String> extensions;
	boolean autoFetch;
	boolean pdbDirectorySplit;
	
	public static final String lineSplit = System.getProperty("file.separator");
	
	FileParsingParameters params;
	SimpleMMcifConsumer consumer;
	
	public static void main(String[] args){
	
		MMCIFFileReader reader = new MMCIFFileReader();
		FileParsingParameters params = new FileParsingParameters();
		reader.setFileParsingParameters(params);
		
		try{
			Structure struc = reader.getStructureById("1m4x");
			System.out.println(struc);
			System.out.println(struc.toPDB());
			
		
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public MMCIFFileReader(){
		extensions    = new ArrayList<String>();
		extensions.add(".cif");
		extensions.add(".mmcif");
		extensions.add(".cif.gz");
		extensions.add(".mmcif.gz");

		UserConfiguration config = new UserConfiguration();
		path = config.getPdbFilePath() ;
		autoFetch     = config.getAutoFetch();
		pdbDirectorySplit = config.isSplit();
		params = new FileParsingParameters();

	}

	public void addExtension(String ext) {
		extensions.add(ext);

	}

	public void clearExtensions(){
		extensions.clear();
	}

	/** Opens filename, parses it and returns
	 * a Structure object .
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

	/** Opens filename, parses it and returns a Structure object.
	 *
	 * @param filename a File object
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(File filename) throws IOException {

		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(filename);

		return parseFromInputStream(inStream);
	}



	private Structure parseFromInputStream(InputStream inStream) throws IOException{

		MMcifParser parser = new SimpleMMcifParser();

		consumer = new SimpleMMcifConsumer();
		   
		consumer.setFileParsingParameters(params);
		
		
		// The Consumer builds up the BioJava - structure object.
		// you could also hook in your own and build up you own data model.
		parser.addMMcifConsumer(consumer);

		parser.parse(new BufferedReader(new InputStreamReader(inStream)));


		// now get the protein structure.
		Structure cifStructure = consumer.getStructure();

		return cifStructure;
	}

	public void setPath(String path) {
		this.path = path;

	}


	public String getPath() {
		return path;
	}

	/** Get a structure by PDB code. This works if a PATH has been set via setPath, or if setAutoFetch has been set to true.
	 *
	 * @param pdbId a 4 letter PDB code.
	 */
	public Structure getStructureById(String pdbId) throws IOException {
		InputStream inStream = getInputStream(pdbId);

		return parseFromInputStream(inStream);
	}

	private InputStream getInputStream(String pdbId) throws IOException{
		
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
					
					if ( params.isUpdateRemediatedFiles()){
						long lastModified = f.lastModified();

						if (lastModified < PDBFileReader.lastRemediationDate) {
							// the file is too old, replace with newer version
							System.out.println("replacing file " + pdbFile +" with latest remediated file from PDB.");
							pdbFile = null;

							return null;
						}
					}
					

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

	public File downloadPDB(String pdbId){

		if ((path == null) || (path.equals(""))){
			System.err.println("you did not set the path in PDBFileReader, don;t know where to write the downloaded file to");
			System.err.println("assuming default location is local directory.");
			path = ".";
		}
				
		File tempFile ;

		if ( pdbDirectorySplit) {
			String middle = pdbId.substring(1,3).toLowerCase();
			String dir = path+lineSplit+middle;
			File directoryCheck = new File (dir);
			if ( ! directoryCheck.exists()){
				directoryCheck.mkdir();
			}

			tempFile = new File(dir+lineSplit+ pdbId.toLowerCase()+".cif.gz");

		} else {

			tempFile = new File(path+lineSplit+pdbId.toLowerCase()+".cif.gz");
		}
		
		
		String ftp = String.format("ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/%s.cif.gz", pdbId.toLowerCase());

		System.out.println("Fetching " + ftp);
		try {
			URL url = new URL(ftp);
			InputStream conn = url.openStream();

			// prepare destination
			System.out.println("writing to " + tempFile);

			FileOutputStream outPut = new FileOutputStream(tempFile);
			GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);

			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(new GZIPInputStream(conn)));
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
			return null;
		}
		return tempFile;
	}

	public boolean isAutoFetch() {
		return autoFetch;
	}


	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;

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



   public FileParsingParameters getFileParsingParameters()
   {
      return params;
   }


   public void setFileParsingParameters(FileParsingParameters params)
   {
     this.params=params;
      
   }

   public SimpleMMcifConsumer getMMcifConsumer(){
	   return consumer;
   }
   
   public void setMMCifConsumer(SimpleMMcifConsumer consumer){
	   this.consumer = consumer;
   }

}
