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

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
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
	private static final Logger s_logger = Logger.getLogger( ZipChemCompProvider.class.getPackage().getName() );

	private final Path m_tempDir;  // Base path where $m_zipRootDir/ will be downloaded to.
	private final Path m_zipRootDir;
	private final Path m_zipFile;
	private final DownloadChemCompProvider m_dlProvider;
	// private final ZipFile m_zipDictionary;
	
	// Missing IDs from library that cannot be download added here to prevent delays.
	private Set<String> unavailable = new HashSet<String>();

	/**
	 * Constructor for a ZipChemCompProvider, a Chemical Component provider that stores chemical component
	 * files in a single zip archive.  Missing chemical components will be downloaded and appended to this
	 * archive.
	 *
	 * @param chemicalComponentDictionaryFile : path to zip archive for chemical components. If non-existent,
	 *                                        an empty zip archive will be created and populated.
	 * @param tempDir : path for temporary directory, (null) will system default property "java.io.tmpdir".
	 * @throws IOException
	 */
	public ZipChemCompProvider(String chemicalComponentDictionaryFile, String tempDir) throws IOException {
		this.m_zipFile = Paths.get(chemicalComponentDictionaryFile);

		// Use a default temporary directory if not passed a value.
		if (tempDir == null || tempDir.equals("")) {
			this.m_tempDir = Paths.get(System.getProperty("java.io.tmpdir"));
		} else {
			this.m_tempDir = Paths.get(tempDir);
		}

		this.m_zipRootDir = Paths.get("chemcomp");

		// Setup an instance of the download chemcomp provider.
		this.m_dlProvider = new DownloadChemCompProvider(m_tempDir.toString());

		// Create a new zip dictionary file if one isn't already created.
		initializeZipFile();
	}
	
	// Create the zip file and setup a fallback DownloadChemCompProvider.
	private void initializeZipFile() throws IOException {
		s_logger.info("Using chemical component dictionary: " + m_zipFile.toString());

		// Check if file exists, if not create.
		final File zipFile = m_zipFile.toFile();
		if (!zipFile.exists()) {
			FileOutputStream f = new FileOutputStream(zipFile);
			ZipOutputStream zip = new ZipOutputStream(new BufferedOutputStream(f));
			try {
				zip.putNextEntry(new ZipEntry("chemcomp/"));
				zip.closeEntry();
			} finally {
				zip.close();
			}
		}
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

		ChemComp cc = null; // Will find in zip or download
		String filename = "chemcomp/" + recordName+".cif.gz";

		// Try to pull from zip.
		try {
			ZipFile m_zip = new ZipFile(m_zipFile.toFile());
			ZipEntry ccEntry = m_zip.getEntry(filename);
			if (ccEntry != null) {
				final InputStream zipStream = m_zip.getInputStream(ccEntry);
				final InputStream inputStream = new GZIPInputStream(zipStream);
				s_logger.fine("reading " + recordName + " from " + m_zipFile);
				final MMcifParser parser = new SimpleMMcifParser();
				final ChemCompConsumer consumer = new ChemCompConsumer();
				parser.addMMcifConsumer(consumer);
				parser.parse(inputStream);

				final ChemicalComponentDictionary dict = consumer.getDictionary();
				cc = dict.getChemComp(recordName);
				inputStream.close();
			}
			m_zip.close();
		} catch (IOException ex) {
			s_logger.severe("Error reading " + m_zipFile.toString() + " " + ex.getMessage());
		}

		// Try to download next.
		if (cc == null) {
			s_logger.info("File "+recordName+" not found in archive. Attempting download from PDB.");
			cc = downloadAndAdd(recordName);
		}

		// Still failed, mark unavailable and return empty chemcomp
		if (cc == null) {
			// Could not find this ChemComp, add this to a list of unavailable.
			s_logger.info("Unable to find or download " + recordName + " - excluding from future searches.");
			unavailable.add(recordName);
			return getEmptyChemComp(recordName);
		}
		return cc;
	}

	/** Use DownloadChemCompProvider to grab a gzipped cif record from the PDB.
	 *  Zip all downloaded cif.gz files into the dictionary.
	 *
	 * @param recordName is the three-letter chemical component code (i.e. residue name).
	 * @return ChemComp matching recordName
     */
	private ChemComp downloadAndAdd(String recordName){
		final ChemComp cc = m_dlProvider.getChemComp(recordName);
		
		if (cc != null && "?".equals(cc.getOne_letter_code())) {
			// Failed to find this chemical component and defaulted back to reduced chemical component.
			unavailable.add(recordName);
			return cc;  // Don't try to add this to the chemcomp.zip.
		}

		//DownloadProvider places files in default temp dir\chemcomp
		final File [] files = finder(m_tempDir.resolve("chemcomp").toString(), "cif.gz");
		if (files != null) {
			ZipAppendUtil.getInstance().addToZipFileSystem(m_zipFile, files, m_zipRootDir);
			for (File f : files) f.delete();
		}
		return cc;
	}

	/**
	 * Cleanup any temporary chemical component files that have been created within tmpdir.
	 * @param tempdir : path to temporary directory for chemical component downloads.
     */
	public static void purgeAllTempFiles(String tempdir) {
		if (null != tempdir) {
			File[] chemCompOutFiles = finder(tempdir, "cif");
			for (File f : chemCompOutFiles) f.delete();

			File[] chemCompTempFiles = finderPrefix(tempdir, "chemcomp.zip");
			for (File f : chemCompTempFiles) f.delete();
		}

		Path chemCompDirPath = Paths.get(System.getProperty(tempdir)).resolve("chemcomp");
		chemCompDirPath.toFile().delete();
	}

	/**
	 * Return an empty ChemComp group for a three-letter resName.
	 * @param resName
	 * @return
	 */
    private ChemComp getEmptyChemComp(String resName){
    	String pdbName = ""; // Empty string is default
    	if (null != resName && resName.length() >= 3) {
    		pdbName = resName.substring(0,3);
    	}
    	final ChemComp comp = new ChemComp();
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
	static private File[] finder( String dirName, final String suffix){
		if (null == dirName || null == suffix) {
			return null;
		}

		final File dir = new File(dirName);
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
	static private File[] finderPrefix( String dirName, final String prefix){
		final File dir = new File(dirName);
		return dir.listFiles(new FilenameFilter() { 
			public boolean accept(File dir, String filename)
			{ return filename.startsWith(prefix); }
		} );
	}
}