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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

import org.biojava.bio.structure.PDBStatus;
import org.biojava.bio.structure.PDBStatus.Status;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.io.util.FileDownloadUtils;
import org.biojava3.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Superclass for classes which download and interact with the PDB's FTP server,
 * specifically {@link PDBFileReader} and {@link MMCIFFileReader}. The basic
 * functionality of downloading structure files from the FTP site is gathered
 * here, making the child classes responsible for only the specific paths and
 * file formats needed.
 * 
 * @author Spencer Bliven
 *
 */
public abstract class LocalPDBDirectory implements StructureIOFile {

	private static final Logger logger = LoggerFactory.getLogger(LocalPDBDirectory.class);

	public static final String DEFAULT_PDB_FILE_SERVER = "ftp.wwpdb.org";
	public static final String PDB_FILE_SERVER_PROPERTY = "PDB.FILE.SERVER";

	/**
	 * Behaviors for when an obsolete structure is requested.
	 * @author Spencer Bliven
	 * @see StructureIOFile#setObsoleteBehavior(ObsoleteBehavior)
	 */
	public static enum ObsoleteBehavior {
		/** Fetch the most recent version of the PDB entry. */
		FETCH_CURRENT,
		/** Fetch the obsolete entry from the PDB archives. */
		FETCH_OBSOLETE,
		/** Throw a StructureException for obsolete entries.*/
		THROW_EXCEPTION;

		public static final ObsoleteBehavior DEFAULT=THROW_EXCEPTION;
	}

	/**
	 * Controls when the class should fetch files from the ftp server
	 * @author Spencer Bliven
	 *
	 */
	public static enum FetchBehavior {
		/** Never fetch from the server; Throw errors for missing files */
		LOCAL_ONLY,
		/** Fetch missing files from the server */
		FETCH_FILES,
		/**
		 * Fetch missing files from the server.
		 * Also force the download of files older than lastRemediationDate.
		 */
		FETCH_REMEDIATED,
		/** For every file, check the server and re-download if a newer version is available. */
		FORCE_DOWNLOAD;

		public static final FetchBehavior DEFAULT = FETCH_REMEDIATED;
	}


	/**
	 * Date of the latest PDB file remediation
	 */
	public static final long LAST_REMEDIATION_DATE ;

	static {

		SimpleDateFormat formatter = new SimpleDateFormat("yyyy/MM/dd", Locale.US);

		long t = 0;
		try {
			Date d = formatter.parse("2011/07/12");
			t = d.getTime();
		} catch (ParseException e){
			logger.warn("Could not parse date: "+e.getMessage());
		}
		LAST_REMEDIATION_DATE = t;
	}

	protected static final String lineSplit = System.getProperty("file.separator");

	private File path;
	private List<String> extensions;

	private String serverName;

	private FileParsingParameters params;

	private ObsoleteBehavior obsoleteBehavior;
	private FetchBehavior fetchBehavior;

	/**
	 * Subclasses should provide default and single-string constructors.
	 * They should use {@link #addExtension(String)} to add one or more extensions.
	 *
	 * <p>If path is null, initialize using the system property/environment variable
	 * {@link UserConfiguration#PDB_DIR}.
	 * @param path Path to the PDB file directory
	 */
	public LocalPDBDirectory(String path) {
		extensions    = new ArrayList<String>();

		params = new FileParsingParameters();

		if( path == null) {
			UserConfiguration config = new UserConfiguration();
			path = config.getPdbFilePath();
			logger.debug("Initialising from system property/environment variable to path: {}", path.toString());
		} else {
			path = FileDownloadUtils.expandUserHome(path);
			logger.debug("Initialising with path {}", path.toString());
		}
		this.path = new File(path);

		this.serverName = System.getProperty(PDBFileReader.PDB_FILE_SERVER_PROPERTY);

		if ( serverName == null || serverName.trim().isEmpty()) {
			serverName = PDBFileReader.DEFAULT_PDB_FILE_SERVER;
			logger.debug("Using default PDB file server {}",serverName);
		} else {
			logger.info("Using PDB file server {} read from system property {}",serverName,PDBFileReader.PDB_FILE_SERVER_PROPERTY);
		}

		fetchBehavior = FetchBehavior.DEFAULT;
		obsoleteBehavior = ObsoleteBehavior.DEFAULT;
	}

	public LocalPDBDirectory() {
		this(null);
	}

	/** 
	 * Sets the path for the directory where PDB files are read/written 
	 */
	@Override
	public void setPath(String p){
		path = new File(FileDownloadUtils.expandUserHome(p)) ;
	}

	/**
	 * Returns the path value.
	 * @return a String representing the path value
	 * @see #setPath
	 *
	 */
	@Override
	public String getPath() {
		return path.toString() ;
	}

	/** define supported file extensions
	 * compressed extensions .Z,.gz do not need to be specified
	 * they are dealt with automatically.
	 */
	@Override
	public void addExtension(String s){
		//System.out.println("add Extension "+s);
		extensions.add(s);
	}

	public List<String> getExtensions() {
		return Collections.unmodifiableList(extensions);
	}

	/** clear the supported file extensions
	 *
	 */
	@Override
	public void clearExtensions(){
		extensions.clear();
	}

	/**
	 * @deprecated Use {@link #getFetchBehavior()}
	 */
	@Deprecated
	@Override
	public boolean isAutoFetch() {
		return fetchBehavior != FetchBehavior.LOCAL_ONLY;
	}

	/**
	 * @deprecated Use {@link #getFetchBehavior()}
	 */
	@Deprecated
	@Override
	public void setAutoFetch(boolean autoFetch) {
		if(autoFetch) {
			setFetchBehavior(FetchBehavior.DEFAULT);
		} else {
			setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		}
	}

	@Override
	public void setFileParsingParameters(FileParsingParameters params){
		this.params= params;
	}

	@Override
	public FileParsingParameters getFileParsingParameters(){
		return params;
	}

	/**
	 * <b>[Optional]</b> This method changes the behavior when obsolete entries
	 * are requested. Current behaviors are:
	 * <ul>
	 * <li>{@link ObsoleteBehavior#THROW_EXCEPTION THROW_EXCEPTION}
	 *   Throw a {@link StructureException} (the default)
	 * <li>{@link ObsoleteBehavior#FETCH_OBSOLETE FETCH_OBSOLETE}
	 *   Load the requested ID from the PDB's obsolete repository
	 * <li>{@link ObsoleteBehavior#FETCH_CURRENT FETCH_CURRENT}
	 *   Load the most recent version of the requested structure
	 * 
	 * <p>This setting may be silently ignored by implementations which do not have
	 * access to the server to determine whether an entry is obsolete, such as
	 * if {@link #isAutoFetch()} is false. Note that an obsolete entry may still be
	 * returned even this is FETCH_CURRENT if the entry is found locally.
	 * 
	 * @param fetchFileEvenIfObsolete Whether to fetch obsolete records
	 * @see #setFetchCurrent(boolean)
	 * @since 4.0.0
	 */
	public void setObsoleteBehavior(ObsoleteBehavior behavior) {
		obsoleteBehavior = behavior;
	}

	/**
	 * Returns how this instance deals with obsolete entries. Note that this
	 * setting may be ignored by some implementations or in some situations,
	 * such as when {@link #isAutoFetch()} is false.
	 * 
	 * <p>For most implementations, the default value is
	 * {@link ObsoleteBehavior#THROW_EXCEPTION THROW_EXCEPTION}.
	 * 
	 * @return The ObsoleteBehavior
	 * @since 4.0.0
	 */
	public ObsoleteBehavior getObsoleteBehavior() {
		return obsoleteBehavior;
	}

	/**
	 * Get the behavior for fetching files from the server
	 * @return
	 */
	public FetchBehavior getFetchBehavior() {
		return fetchBehavior;
	}
	/**
	 * Set the behavior for fetching files from the server.
	 * This replaces the {@link #setAutoFetch(boolean)} method with a more
	 * extensive set of options.
	 * @param fetchBehavior
	 */
	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior;
	}


	/** opens filename, parses it and returns
	 * a Structure object .
	 * @param filename  a String
	 * @return the Structure object
	 * @throws IOException for errors reading or writing the file
	 */
	@Override
	public Structure getStructure(String filename) throws IOException
	{
		filename = FileDownloadUtils.expandUserHome(filename);
		File f = new File(filename);
		return getStructure(f);

	}

	public Structure getStructure(URL u) throws IOException{
		InputStreamProvider isp = new InputStreamProvider();
		InputStream inStream = isp.getInputStream(u);
		return getStructure(inStream);
	}

	/** opens filename, parses it and returns a Structure object
	 *
	 * @param filename a File object
	 * @return the Structure object
	 * @throws IOException ...
	 */
	@Override
	public Structure getStructure(File filename) throws IOException {
		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(filename);

		return getStructure(inStream);
	}


	/** Get a structure by PDB code. This works if a PATH has been set via setPath, or if setAutoFetch has been set to true.
	 *
	 * @param pdbId a 4 letter PDB code.
	 */
	@Override
	public Structure getStructureById(String pdbId) throws IOException {
		InputStream inStream = getInputStream(pdbId);

		return getStructure(inStream);
	}

	/**
	 * Handles the actual parsing of the file into a Structure object.
	 * @param inStream
	 * @return
	 * @throws IOException
	 */
	public abstract Structure getStructure(InputStream inStream) throws IOException;

	/**
	 * Load or download the specified structure and return it as an InputStream
	 * for direct parsing.
	 * @param pdbId
	 * @return
	 * @throws IOException
	 */
	protected InputStream getInputStream(String pdbId) throws IOException{

		if ( pdbId.length() < 4)
			throw new IOException("The provided ID does not look like a PDB ID : " + pdbId);

		// Check existing
		File file = getLocalFile(pdbId);

		// Redownload
		if(file == null ) {
			file = downloadStructure(pdbId);
		}

		if(!file.exists()) {
			throw new IOException("Structure "+pdbId+" not found and unable to download.");
		}

		InputStreamProvider isp = new InputStreamProvider();

		InputStream inputStream = isp.getInputStream(file);

		return inputStream;
	}

	/**
	 * Download a structure, but don't parse it yet or store it in memory.
	 * 
	 * Used to pre-fetch large numbers of structures.
	 * @param pdbId
	 * @throws IOException 
	 */
	public void prefetchStructure(String pdbId) throws IOException {
		if ( pdbId.length() < 4)
			throw new IOException("The provided ID does not look like a PDB ID : " + pdbId);

		// Check existing
		File file = getLocalFile(pdbId);

		// Redownload
		if(file == null ) {
			file = downloadStructure(pdbId);
		}

		if(!file.exists()) {
			throw new IOException("Structure "+pdbId+" not found and unable to download.");
		}
	}

	/**
	 * Downloads an MMCIF file from the PDB to the local path
	 * @param pdbId
	 * @return The file, or null if it was unavailable for download
	 * @throws IOException for errors downloading or writing, or if the
	 *  fetchBehavior is LOCAL_ONLY
	 */
	@SuppressWarnings("deprecation") //for isUpdateRemediatedFiles()
	protected File downloadStructure(String pdbId) throws IOException{
		// decide whether download is required
		File existing =  getLocalFile(pdbId);
		switch(fetchBehavior) {
		case LOCAL_ONLY:
			throw new IOException(String.format("Structure {} not found in {} "
					+ "and configured not to download.",pdbId,getPath()));
		case FETCH_FILES:
			// Use existing if present
			if( existing != null) {
				// Respect deprecated remediation parameter
				if(getFileParsingParameters().isUpdateRemediatedFiles()) {
					long lastModified = existing.lastModified();

					if (lastModified < LAST_REMEDIATION_DATE) {
						// the file is too old, replace with newer version
						logger.warn("Replacing file " + existing +" with latest remediated file from PDB.");
						break;
					}
				}
				return existing;
			}
			break;
		case FETCH_REMEDIATED:
			// Use existing if present and recent enough
			if( existing != null) {
				long lastModified = existing.lastModified();

				if (lastModified < LAST_REMEDIATION_DATE) {
					// the file is too old, replace with newer version
					logger.warn("Replacing file " + existing +" with latest remediated file from PDB.");
					break;
				} else {
					return existing;
				}
			}
		case FORCE_DOWNLOAD:
			//TODO check if file is up-to-date (SB 2015-01)
			break; //force re-download of everything
		}

		// Force the download now
		if(obsoleteBehavior == ObsoleteBehavior.FETCH_CURRENT) {
			String current = PDBStatus.getCurrent(pdbId);

			if(current == null) {
				// either an error or there is not current entry
				current = pdbId;
			}
			String path;
			path = String.join("/", getSplitDirPath());
			return downloadStructure(current, path,false);
		} else if(obsoleteBehavior == ObsoleteBehavior.FETCH_OBSOLETE
				&& PDBStatus.getStatus(pdbId) == Status.OBSOLETE) {
			String path = String.join("/", getObsoleteDirPath());
			return downloadStructure(pdbId, path,true);
		} else {
			String path;
			path = String.join("/", getSplitDirPath());
			return downloadStructure(pdbId, path, false);
		}
	}

	/**
	 * Download a file from the ftp server, replacing any existing files if needed
	 * @param pdbId PDB ID
	 * @param pathOnServer Path on the FTP server, e.g. data/structures/divided/pdb
	 * @param obsolete Whether or not file should be saved to the obsolete location locally
	 * @return
	 * @throws IOException
	 */
	private File downloadStructure(String pdbId, String pathOnServer, boolean obsolete) throws IOException{
		File dir = getDir(pdbId,obsolete);
		File realFile = new File(dir,getFilename(pdbId));

		String ftp = String.format("ftp://%s/pub/pdb/%s/%s/%s", 
				serverName, pathOnServer, pdbId.substring(1,3).toLowerCase(), getFilename(pdbId));

		logger.info("Fetching " + ftp);
		logger.info("Writing to "+ realFile);

		URL url = new URL(ftp);

		FileDownloadUtils.downloadGzipCompressedFile(url, realFile);

		return realFile;
	}

	/**
	 * Gets the directory in which the file for a given MMCIF file would live,
	 * creating it if necessary.
	 * 
	 * The obsolete parameter is necessary to avoid additional server queries.
	 * @param pdbId
	 * @param obsolete Whether the pdbId is obsolete or not
	 * @return File pointing to the directory, 
	 */
	protected File getDir(String pdbId, boolean obsolete) {

		File dir = null;

		if (obsolete) {
			// obsolete is always split
			String middle = pdbId.substring(1,3).toLowerCase();
			dir = new File(path, String.join(lineSplit, getObsoleteDirPath()) + lineSplit + middle);
		} else {
			String middle = pdbId.substring(1,3).toLowerCase();
			dir = new File(path, String.join(lineSplit, getSplitDirPath()) + lineSplit + middle);
		}


		if (!dir.exists()) {
			boolean success = dir.mkdirs();
			if (!success) logger.error("Could not create mmCIF dir {}",dir.toString());
		}

		return dir;
	}

	/**
	 * Searches for previously downloaded files
	 * @param pdbId
	 * @return A file pointing to the existing file, or null if not found
	 */
	protected File getLocalFile(String pdbId) {

		// Search for existing files

		// Search directories:
		// 1) LOCAL_MMCIF_SPLIT_DIR/<middle>/(pdb)?<pdbId>.<ext>
		// 2) LOCAL_MMCIF_ALL_DIR/<middle>/(pdb)?<pdbId>.<ext>
		LinkedList<File> searchdirs = new LinkedList<File>();
		String middle = pdbId.substring(1,3).toLowerCase();

		File splitdir = new File(new File(getPath(), String.join(lineSplit, getSplitDirPath()) ), middle);
		searchdirs.add(splitdir);
		// Search obsolete files if requested
		if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE) {
			File obsdir = new File(getPath(), String.join(lineSplit, getObsoleteDirPath()));
			searchdirs.add(obsdir);
		}

		// valid prefixes before the <pdbId> in the filename
		String[] prefixes = new String[] {"", "pdb"};

		for( File searchdir :searchdirs){
			for( String prefix : prefixes) {
				for(String ex : getExtensions() ){
					File f = new File(searchdir,prefix + pdbId + ex) ;
					if ( f.exists()) {
						return f;
					}
				}
			}
		}
		//Not found
		return null;
	}

	protected boolean checkFileExists(String pdbId){
		File path =  getLocalFile(pdbId);
		if ( path != null)
			return true;
		return false;
	}

	/**
	 * Converts a PDB ID into a filename with the proper extension
	 * @param pdbId
	 * @return The filename, e.g. "4hhb.pdb.gz"
	 */
	protected abstract String getFilename(String pdbId);

	protected abstract String[] getSplitDirPath();
	protected abstract String[] getObsoleteDirPath();
}
