package org.biojava.bio.structure.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;



/** Provides structures based on locally cached PDB files. The path where to cache this can be configured in FileParsingParameters
 * 
 * 
 * @author Andreas Prlic
 *
 */
public class LocalCacheStructureProvider implements StructureProvider{

	private static final Logger logger = LoggerFactory.getLogger(LocalCacheStructureProvider.class);
	
	FileParsingParameters params ;
	
	String path;
	public static final String lineSplit = System.getProperty("file.separator");
	
	private static final String BIOL_DIR  = "BIOL_UNITS";
	
	public LocalCacheStructureProvider(){
		params = new FileParsingParameters();
		
		// set the path to cache files

		UserConfiguration config = new UserConfiguration();
		
		setPath(config.getPdbFilePath());
		
	}
	
	/** directory where to find PDB files */
	public void setPath(String p){
		path = p ;
		
		if ( !(path.endsWith(lineSplit) ) )
			path = path + lineSplit;
				
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
	
	@Override
	public Structure getStructureById(String pdbId) throws IOException {
		
		PDBFileReader reader = getPdbFilereader();
		return reader.getStructureById(pdbId);
	}
	
	private PDBFileReader getPdbFilereader(){
		PDBFileReader reader = new PDBFileReader(path);
		reader.setFetchBehavior(FetchBehavior.DEFAULT);
		reader.setFileParsingParameters(params);
		return reader;
	}

	/** Returns the first available biological unit for the PDB ID.
	 * Note: the biological units are represented using multiple models.
	 * You need to work with all models to get the full set of atoms.
	 * 
	 * @throws IOException 
	 * 
	 */
	@Override
	public Structure getBiologicalUnit(String pdbId) throws StructureException, IOException {
		if (pdbId == null || pdbId.length()< 4)
			throw new StructureException("This does not look like a valid PDB ID! (" + pdbId + ")");
	
		String middle = pdbId.substring(1,3).toLowerCase();
		String dir = path + BIOL_DIR;
		File tmp1 = new File(dir);
		if ( ! tmp1.exists()){
			tmp1.mkdir();
		}
		dir = dir+lineSplit+middle;
		File directoryCheck = new File (dir);
		if ( ! directoryCheck.exists()){
			directoryCheck.mkdir();
		}
		
	
		File tempFile =new File(dir + lineSplit+"pdb"+ pdbId.toLowerCase()+".ent.gz");
	
		if (! tempFile.exists()){
			// we need to download the file
			downloadBiolUnit(pdbId.toLowerCase(),tempFile);
		}
		
		
		PDBFileReader reader = getPdbFilereader();
		
		return reader.getStructure(tempFile);
		
	}
	
	
	
	/** Download the biological unit file
	 * 
	 * @param pdbId the PDB ID
	 * @param localFile where to store the file locally
	 */
	private void downloadBiolUnit(String pdbId, File localFile) throws IOException{
		String u = "http://www.rcsb.org/pdb/files/%s.pdb1.gz";
		
		String ur = String.format(u,pdbId);
		
		logger.info("Fetching " + ur);

		// prepare destination
		logger.info("writing to " + localFile);


		URL url = new URL(ur);

		InputStream uStream = url.openStream();
		InputStream conn = new GZIPInputStream(uStream);


		FileOutputStream outPut = new FileOutputStream(localFile);
		GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
		PrintWriter pw = new PrintWriter(gzOutPut);

		BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(conn));
		String line;
		while ((line = fileBuffer.readLine()) != null) {
			pw.println(line);
		}
		pw.flush();
		pw.close();

		outPut.flush();
		outPut.close();
		conn.close();
		uStream.close();


	}
	
	@Override
	public void setFileParsingParameters(FileParsingParameters params) {
		this.params = params;
	}

	@Override
	public FileParsingParameters getFileParsingParameters() {
		return params;
	}

	
	
}
