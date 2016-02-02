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
package org.biojava.nbio.structure.io.mmcif;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;

/** This chemical component provider retrieves and caches chemical component definition files from a zip archive specified in its construction.
 * 	If the archive does not contain the record, an attempt is made to download it using DownloadChemCompProvider. The downloaded file is then added to the archive.
 *	It is up to the calling function to handle exceptions raised by constructing with a non-existent archive. 
 *
 * @author edlunde
 * @author larsonm
 * @since 12/05/12
 * 
 */
public class ZipChemCompProvider implements ChemCompProvider{

	private static final String FILE_SEPARATOR = System.getProperty("file.separator");
	private static ZipChemCompProvider s_provider = null;
	private String tempdir = "";  // Base path where $zipRootDir/ will be downloaded to.
	private String zipRootDir = "chemcomp/"; 
	private ZipFile zipDictionary;
	private String zipFile;
	private static final Logger s_logger = Logger.getLogger( ZipChemCompProvider.class.getPackage().getName() );
	private DownloadChemCompProvider dlProvider;
	
	// Missing IDs from library that cannot be download added here to prevent delays.
	private Set<String> unavailable = new HashSet<String>();

	/**
	 * 
	 * @param chemicalComponentDictionaryFile : zip filename
	 * @param tempDir : directory of the zip filename
	 * @throws IOException
	 */
	public ZipChemCompProvider(String chemicalComponentDictionaryFile, String tempDir) throws IOException{
		this.zipFile = chemicalComponentDictionaryFile;
		this.tempdir = tempDir;
		
		init();
	}
	
	// Create the zip file and setup a fallback DownloadChemCompProvider.
	private void init() throws IOException {
		try{
			s_logger.info("Using chemical component dictionary: " + zipFile.toString());
			
			// Check if file exists, if not create.
			File zipFile = new File(this.zipFile);
			if (!zipFile.exists()) {
				createNewZip(zipFile);
			}
			
			this.zipDictionary = new ZipFile(this.zipFile);
		}catch(ZipException e){
			s_logger.info(e.getMessage());
		}
		
		// Setup an instance of the download chemcomp provider.
		dlProvider = new DownloadChemCompProvider(this.tempdir);
	}

	/**
	 * Get a singleton instance of the provider.
	 * @return
	 */
	public static ZipChemCompProvider getInstance() throws IOException {
		ZipChemCompProvider provider = s_provider;
		if (provider == null){
			s_provider = provider = new ZipChemCompProvider();
		}

		return provider;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.io.mmcif.ChemCompProvider#getChemComp(java.lang.String)
	 * 
	 * @param recordName : three letter PDB name for a residue
	 * @return ChemComp from .zip or ChemComp from repository.  Will return empty ChemComp when unable to find a residue and will return null if not provided a valid recordName.
	 */
	@Override
	public ChemComp getChemComp(String recordName) {
		
		if (null == recordName) return null;
		
		// handle non-existent ChemComp codes and do not repeatedly attempt to add these.
		for (String str : unavailable) {
			if (recordName.equals(str)) return getEmptyChemComp(recordName);
		}
		
		try {
			ZipEntry ccEntry = zipDictionary.getEntry(zipRootDir+recordName+".cif.gz");
			
			//read single gzipped cif from archive and write to temporary file 
			InputStream zipStream = zipDictionary.getInputStream(ccEntry);
			InputStream inputStream = new GZIPInputStream(zipStream);
			s_logger.fine("reading "+recordName+" from " + zipFile); 
			MMcifParser parser = new SimpleMMcifParser();
			ChemCompConsumer consumer = new ChemCompConsumer();
			parser.addMMcifConsumer(consumer);
			parser.parse(inputStream);
			
			ChemicalComponentDictionary dict = consumer.getDictionary();
			ChemComp cc = dict.getChemComp(recordName);
			zipStream.close();
			if (cc != null){
				return cc;
			}
		}catch(NullPointerException npe){
			s_logger.info(npe.getMessage()+"\n"+"File "+recordName+" not found in archive. Attempting download from PDB.");
		}
		catch (Exception e) {
			s_logger.info(e.getMessage());
		}
		
		//if file isn't found in archive, download it/add to archive
		ChemComp cc = downloadAndAdd(recordName);
		
		if (cc == null) {
			// Could not find this ChemComp, add this to a list of unavailable. 
			unavailable.add(recordName);
			return getEmptyChemComp(recordName);
		}
		
		return cc;

	}
	/**	Use DownloadChemCompProvider to grab a gzipped cif record from the PDB. Zip all downloaded cif.gz files into the dictionary. 
	 * */ 
	public ChemComp downloadAndAdd(String recordName){
		ChemComp cc = dlProvider.getChemComp(recordName);
		
		if (cc != null && "?".equals(cc.getOne_letter_code())) {
			// Failed to find this chemical component and defaulted back to reduced chemical component.
			unavailable.add(recordName);
			return cc;  // Don't try to add this to the chemcomp.zip.
		}
		
		if(ZipAppendUtil.getInstance().isBusy()) {
			return cc;
		}
		
		String tmpdir = this.tempdir;

		// Fall-back to java temporary directory if not set.
		if (tmpdir == null || tmpdir.equals("")) {
			tmpdir = System.getProperty("java.io.tmpdir");	
		}
		if ( !(tmpdir.endsWith(FILE_SEPARATOR) ) ) {
			tmpdir += FILE_SEPARATOR;
		}
		tmpdir+="chemcomp";
		
		final File [] files;
		final File source;
		try{
			files = finder(tmpdir, "cif.gz");	//DownloadProvider places files in default temp dir\chemcomp
			source = new File(zipDictionary.getName());
		}catch(NullPointerException npe){
			return cc;
		}

		// Adding files could run on a separate thread, but more care is needed to prevent race
		// conditions such as having asssembled a list of files to add/adding files while still writing an individual
		// file to the temporary directory.  Keeping this on one thread to prevent race conditions.
		
		// TODO: assess race conditions and make reassembling the zip a background task.
		
		//Thread updateArchiveThread = new Thread(new Runnable(){
		//	public void run(){
				
		if(! ZipAppendUtil.getInstance().isBusy()) {
			ZipAppendUtil.getInstance().addToZipArchive(source, files, zipRootDir);
			for (File f : files){
				f.delete();
			}
			try {
				zipDictionary.close();
			} catch (IOException e) {
				s_logger.info(e.getMessage());
			}
		}
		//	}
		//});
		//updateArchiveThread.start();
		
		return cc;
	}

	/**
	 * Cleanup the temporary files that have been created within tmpdir.
	 */
	public void purgeAllTempFiles() {
		File[] ccOutFiles = finder(System.getProperty("java.io.tmpdir"),"cif");
		for(File f : ccOutFiles) f.delete();
		File[] chemcompTempFiles = finderPrefix(System.getProperty("java.io.tmpdir"), "chemcomp.zip");
		for(File f : chemcompTempFiles)f.delete();
		File chemcompDir = new File(System.getProperty("java.io.tmpdir") + FILE_SEPARATOR +"chemcomp");
		chemcompDir.delete();
	}

	/**
	 * Return an empty ChemComp group for a three-letter resName.
	 * @param resName
	 * @return
	 */
    private ChemComp getEmptyChemComp(String resName){
    	String pdbName = "";
    	if (null != resName && resName.length() >= 3) {
    		pdbName = resName.substring(0,3);
    	}
    	ChemComp comp = new ChemComp();   
    	comp.setOne_letter_code("?");
    	comp.setThree_letter_code(pdbName);
    	comp.setPolymerType(PolymerType.unknown);
    	comp.setResidueType(ResidueType.atomn);
    	return comp;
    }
    
	/**
	 * Return File(s) in dirName that match suffix.
	 * @param dirName
	 * @param suffix
	 * @return
	 */
	private File[] finder( String dirName, final String suffix){
		File dir = new File(dirName);
		return dir.listFiles(new FilenameFilter() { 
			public boolean accept(File dir, String filename)
			{ return filename.endsWith(suffix); }
		} );

	}
	
	/**
	 * Return File(s) in dirName that match prefix.
	 * @param dirName
	 * @param prefix
	 * @return
	 */
	private File[] finderPrefix( String dirName, final String prefix){
		File dir = new File(dirName);
		return dir.listFiles(new FilenameFilter() { 
			public boolean accept(File dir, String filename)
			{ return filename.startsWith(prefix); }
		} );
	}
	
	/**
	 * Construct a provider for a singleton instance.
	 */
	private ZipChemCompProvider() throws IOException {
		this.tempdir = System.getProperty("java.io.tmpdir");	
		init();
	}
	
	private void createNewZip(File f) throws IOException {
	      final ZipOutputStream out = new ZipOutputStream(new FileOutputStream(f));
	      out.close();
	}
}
