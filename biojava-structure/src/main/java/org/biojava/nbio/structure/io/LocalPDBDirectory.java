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
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.PDBStatus;
import org.biojava.nbio.structure.PDBStatus.Status;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.rcsb.mmtf.utils.CodecUtils;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

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

	/**
	 * The default server name, prefixed by the protocol string (http:// or ftp://).
	 * Note that we don't support file stamp retrieving for ftp protocol, thus some of the
	 * fetch modes will not work properly with ftp protocol
	 */
	public static final String DEFAULT_PDB_FILE_SERVER = "http://ftp.wwpdb.org";
	public static final String PDB_FILE_SERVER_PROPERTY = "PDB.FILE.SERVER";

	/**
	 * Behaviors for when an obsolete structure is requested.
	 * @author Spencer Bliven
	 * @see LocalPDBDirectory#setObsoleteBehavior(ObsoleteBehavior)
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
		/** Fetch missing files from the server. Don't check for outdated files */
		FETCH_FILES,
		/**
		 * Fetch missing files from the server, also fetch if file present but older than the
		 * server file.
		 * This requires always querying the server for the last modified time of the file, thus
		 * it adds an overhead to getting files from cache.
		 */
		FETCH_IF_OUTDATED,
		/**
		 * Fetch missing files from the server.
		 * Also force the download of files older than {@value #LAST_REMEDIATION_DATE_STRING}.
		 */
		FETCH_REMEDIATED,
		/** For every file, force downloading from the server */
		FORCE_DOWNLOAD;

		public static final FetchBehavior DEFAULT = FETCH_REMEDIATED;
	}

	/**
	 * Date of the latest PDB file remediation
	 */
	public static final long LAST_REMEDIATION_DATE ;
	private static final String LAST_REMEDIATION_DATE_STRING = "2011/07/12";

	static {

		SimpleDateFormat formatter = new SimpleDateFormat("yyyy/MM/dd", Locale.US);

		long t = 0;
		try {
			Date d = formatter.parse(LAST_REMEDIATION_DATE_STRING);
			t = d.getTime();
		} catch (ParseException e){
			logger.error("Unexpected error! could not parse LAST_REMEDIATION_DATE: "+e.getMessage());
		}
		LAST_REMEDIATION_DATE = t;
	}

	protected static final String lineSplit = System.getProperty("file.separator");

	private File path;
	private List<String> extensions;

	/**
	 * The server name, prefixed by the protocol string (http:// or ftp://).
	 * Note that we don't support file stamp retrieving for ftp protocol, thus some of the
	 * fetch modes will not work properly with ftp protocol
	 */
	private String serverName;

	private FileParsingParameters params;

	private ObsoleteBehavior obsoleteBehavior;
	private FetchBehavior fetchBehavior;


	// Cache results of get*DirPath()
	private String splitDirURL; // path on the server, starting with a slash and ending before the 2-char split directories
	private String obsoleteDirURL;
	private File splitDirPath; // path to the directory before the 2-char split
	private File obsoleteDirPath;

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
			logger.debug("Initialising from system property/environment variable to path: {}", path);
		} else {
			path = FileDownloadUtils.expandUserHome(path);
			logger.debug("Initialising with path {}", path);
		}
		this.path = new File(path);

		this.serverName = getServerName();

		// Initialize splitDirURL,obsoleteDirURL,splitDirPath,obsoleteDirPath
		initPaths();

		fetchBehavior = FetchBehavior.DEFAULT;
		obsoleteBehavior = ObsoleteBehavior.DEFAULT;
	}

	public LocalPDBDirectory() {
		this(null);
	}

	/**
	 * Sets the path for the directory where PDB files are read/written
	 */
	public void setPath(String p){
		path = new File(FileDownloadUtils.expandUserHome(p)) ;
		initPaths();
	}

	/**
	 * Returns the path value.
	 * @return a String representing the path value
	 * @see #setPath
	 *
	 */
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

	@Override
	public List<String> getExtensions() {
		return Collections.unmodifiableList(extensions);
	}

	/** clear the supported file extensions
	 *
	 */
	public void clearExtensions(){
		extensions.clear();
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

	@Override
	public Structure getStructure(File filename) throws IOException {
		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(filename);

		return getStructure(inStream);
	}


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

		if ( pdbId.length() != 4)
			throw new IOException("The provided ID does not look like a PDB ID : " + pdbId);

		// Check existing
		File file = downloadStructure(pdbId);

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
		if ( pdbId.length() != 4)
			throw new IOException("The provided ID does not look like a PDB ID : " + pdbId);

		// Check existing
		File file = downloadStructure(pdbId);

		if(!file.exists()) {
			throw new IOException("Structure "+pdbId+" not found and unable to download.");
		}
	}

	/**
	 * Attempts to delete all versions of a structure from the local directory.
	 * @param pdbId
	 * @return True if one or more files were deleted
	 */
	public boolean deleteStructure(String pdbId){
		boolean deleted = false;
		// Force getLocalFile to check in obsolete locations
		ObsoleteBehavior obsolete = getObsoleteBehavior();
		setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);

		try {
			File existing = getLocalFile(pdbId);
			while(existing != null) {
				assert(existing.exists()); // should exist unless concurrency problems

				if( getFetchBehavior() == FetchBehavior.LOCAL_ONLY) {
					throw new RuntimeException("Refusing to delete from LOCAL_ONLY directory");
				}

				// delete file
				boolean success = existing.delete();
				if(success) {
					logger.info("Deleting "+existing.getAbsolutePath());
				}
				deleted = deleted || success;

				// delete parent if empty
				File parent = existing.getParentFile();
				if(parent != null) {
					success = parent.delete();
					if(success) {
						logger.info("Deleting "+parent.getAbsolutePath());
					}
				}

				existing = getLocalFile(pdbId);
			}
			return deleted;

		} finally {
			setObsoleteBehavior(obsolete);
		}
	}

	/**
	 * Downloads an MMCIF file from the PDB to the local path
	 * @param pdbId
	 * @return The file, or null if it was unavailable for download
	 * @throws IOException for errors downloading or writing, or if the
	 *  fetchBehavior is {@link FetchBehavior#LOCAL_ONLY}
	 */
	protected File downloadStructure(String pdbId) throws IOException{
		if ( pdbId.length() != 4)
			throw new IOException("The provided ID does not look like a PDB ID : " + pdbId);

		// decide whether download is required
		File existing =  getLocalFile(pdbId);
		switch(fetchBehavior) {
		case LOCAL_ONLY:
			if( existing == null ) {
				throw new IOException(String.format("Structure %s not found in %s "
						+ "and configured not to download.",pdbId,getPath()));
			} else {
				return existing;
			}
		case FETCH_FILES:
			// Use existing if present
			if( existing != null) {
				return existing;
			}
			// existing is null, downloadStructure(String,String,boolean,File) will download it
			break;
		case FETCH_IF_OUTDATED:
			// here existing can be null or not:
			// existing == null : downloadStructure(String,String,boolean,File) will download it
			// existing != null : downloadStructure(String,String,boolean,File) will check its date and download if older
			break;
		case FETCH_REMEDIATED:
			// Use existing if present and recent enough
			if( existing != null) {
				long lastModified = existing.lastModified();

				if (lastModified < LAST_REMEDIATION_DATE) {
					// the file is too old, replace with newer version
					logger.warn("Replacing file {} with latest remediated (remediation of {}) file from PDB.",
							existing, LAST_REMEDIATION_DATE_STRING);
					existing = null;
					break;
				} else {
					return existing;
				}
			}
		case FORCE_DOWNLOAD:
			// discard the existing file to force redownload
			existing = null; // downloadStructure(String,String,boolean,File) will download it
			break;
		}

		// Force the download now
		if(obsoleteBehavior == ObsoleteBehavior.FETCH_CURRENT) {
			String current = PDBStatus.getCurrent(pdbId);

			if(current == null) {
				// either an error or there is not current entry
				current = pdbId;
			}
			return downloadStructure(current, splitDirURL,false, existing);
		} else if(obsoleteBehavior == ObsoleteBehavior.FETCH_OBSOLETE
				&& PDBStatus.getStatus(pdbId) == Status.OBSOLETE) {
			return downloadStructure(pdbId, obsoleteDirURL, true, existing);
		} else {
			return downloadStructure(pdbId, splitDirURL, false, existing);
		}
	}

	/**
	 * Download a file from the ftp server, replacing any existing files if needed
	 * @param pdbId PDB ID
	 * @param pathOnServer Path on the FTP server, e.g. data/structures/divided/pdb
	 * @param obsolete Whether or not file should be saved to the obsolete location locally
	 * @param existingFile if not null and checkServerFileDate is true, the last modified date of the
	 * server file and this file will be compared to decide whether to download or not
	 * @return
	 * @throws IOException
	 */
	private File downloadStructure(String pdbId, String pathOnServer, boolean obsolete, File existingFile)
			throws IOException{

		File dir = getDir(pdbId,obsolete);
		File realFile = new File(dir,getFilename(pdbId));

		String ftp;
		
		if (getFilename(pdbId).endsWith(".mmtf.gz")){			
			ftp = CodecUtils.getMmtfEntryUrl(pdbId, true, false);
		} else {
			ftp = String.format("%s%s/%s/%s",
			serverName, pathOnServer, pdbId.substring(1,3).toLowerCase(), getFilename(pdbId));
		}

		URL url = new URL(ftp);

		Date serverFileDate = null;
		if (existingFile!=null) {

			serverFileDate = getLastModifiedTime(url);

			if (serverFileDate!=null) {
				if (existingFile.lastModified()>=serverFileDate.getTime()) {
					return existingFile;
				} else {
					// otherwise we go ahead and download, warning about it first
					logger.warn("File {} is outdated, will download new one from PDB (updated on {})",
							existingFile, serverFileDate.toString());
				}
			} else {
				logger.warn("Could not determine if file {} is outdated (could not get timestamp from server). Will force redownload", existingFile);
			}
		}

		logger.info("Fetching " + ftp);
		logger.info("Writing to "+ realFile);

		FileDownloadUtils.downloadFile(url, realFile);

		// Commented out following code used for setting the modified date to the downloaded file - JD 2015-01-15
		// The only reason to have it was in order to get an rsync-like behavior, respecting the timestamps
		// but the issue is that it would make the FETCH_REMEDIATED mode redownload files with timestamps before
		// the remediation.
		//if (serverFileDate==null)
		//	serverFileDate = getLastModifiedTime(url);
		//
		//if (serverFileDate!=null) {
		//	logger.debug("Setting modified time of downloaded file {} to {}",realFile,serverFileDate.toString());
		//	realFile.setLastModified(serverFileDate.getTime());
		//} else {
		//	logger.warn("Unknown modified time of file {}, will set its modified time to now.", ftp);
		//}


		return realFile;
	}

	/**
	 * Get the last modified time of the file in given url by retrieveing the "Last-Modified" header.
	 * Note that this only works for http URLs
	 * @param url
	 * @return the last modified date or null if it couldn't be retrieved (in that case a warning will be logged)
	 */
	private Date getLastModifiedTime(URL url) {

		// see http://stackoverflow.com/questions/2416872/how-do-you-obtain-modified-date-from-a-remote-file-java
		Date date = null;
		try {
			String lastModified = url.openConnection().getHeaderField("Last-Modified");
			logger.debug("Last modified date of server file ({}) is {}",url.toString(),lastModified);


			if (lastModified!=null) {

				try {
					date = new SimpleDateFormat("E, d MMM yyyy HH:mm:ss Z", Locale.ENGLISH).parse(lastModified);
				} catch (ParseException e) {
					logger.warn("Could not parse last modified time from string '{}', no last modified time available for file {}",
							lastModified, url.toString());
					// this will return null
				}

			}
		} catch (IOException e) {
			logger.warn("Problems while retrieving last modified time for file {}", url.toString());
		}
		return date;

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
			dir = new File(obsoleteDirPath, middle);
		} else {
			String middle = pdbId.substring(1,3).toLowerCase();
			dir = new File(splitDirPath, middle);
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
	public File getLocalFile(String pdbId) {

		// Search for existing files

		// Search directories:
		// 1) LOCAL_MMCIF_SPLIT_DIR/<middle>/(pdb)?<pdbId>.<ext>
		// 2) LOCAL_MMCIF_ALL_DIR/<middle>/(pdb)?<pdbId>.<ext>
		LinkedList<File> searchdirs = new LinkedList<File>();
		String middle = pdbId.substring(1,3).toLowerCase();

		File splitdir = new File(splitDirPath, middle);
		searchdirs.add(splitdir);
		// Search obsolete files if requested
		if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE) {
			File obsdir = new File(obsoleteDirPath,middle);
			searchdirs.add(obsdir);
		}

		// valid prefixes before the <pdbId> in the filename
		String[] prefixes = new String[] {"", "pdb"};

		for( File searchdir :searchdirs){
			for( String prefix : prefixes) {
				for(String ex : getExtensions() ){
					File f = new File(searchdir,prefix + pdbId.toLowerCase() + ex) ;
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
	 * Return the String with the PDB server name, including the leading protocol
	 * String (http:// or ftp://).
	 * The server name will be by default the value {@value #DEFAULT_PDB_FILE_SERVER} or the one
	 * read from system property {@value #PDB_FILE_SERVER_PROPERTY}
	 *
	 * @return the server name including the leading protocol string
	 */
	public static String getServerName() {
		String name = System.getProperty(PDB_FILE_SERVER_PROPERTY);

		if ( name == null || name.trim().isEmpty()) {
			name = DEFAULT_PDB_FILE_SERVER;
			logger.debug("Using default PDB file server {}",name);
		} else {
			if (!name.startsWith("http://") && !name.startsWith("ftp://") && !name.startsWith("https://")) {
				logger.warn("Server name {} read from system property {} does not have a leading protocol string. Adding http:// to it", name, PDB_FILE_SERVER_PROPERTY);
				name = "http://"+name;
			}
			logger.info("Using PDB file server {} read from system property {}", name, PDB_FILE_SERVER_PROPERTY);
		}
		return name;
	}

	/**
	 * Should be called whenever any of the path variables change.
	 * Thus, if {@link getSplitDirPath()} or {@link getObsoleteDirPath()}
	 * depend on anything, they should call this function when that thing
	 * changes (possibly including at the end of the constructor).
	 */
	protected void initPaths() {
		// Hand-rolled String.join(), for java 6
		String[] split = getSplitDirPath();
		String[] obsolete = getObsoleteDirPath();

		//URLs are joined with '/'
		StringBuilder splitURL = new StringBuilder("/pub/pdb");
		for(int i=0;i<split.length;i++) {
			splitURL.append("/");
			splitURL.append(split[i]);
		}
		StringBuilder obsoleteURL = new StringBuilder("/pub/pdb");
		for(int i=0;i<obsolete.length;i++) {
			obsoleteURL.append("/");
			obsoleteURL.append(obsolete[i]);
		}

		splitDirURL = splitURL.toString();
		obsoleteDirURL = obsoleteURL.toString();


		//Files join themselves iteratively
		splitDirPath = path;
		for(int i=0;i<split.length;i++) {
			splitDirPath = new File(splitDirPath,split[i]);
		}
		obsoleteDirPath = path;
		for(int i=0;i<obsolete.length;i++) {
			obsoleteDirPath = new File(obsoleteDirPath,obsolete[i]);
		}
	}

	/**
	 * Converts a PDB ID into a filename with the proper extension
	 * @param pdbId
	 * @return The filename, e.g. "4hhb.pdb.gz"
	 */
	protected abstract String getFilename(String pdbId);

	/**
	 * Location of split files within the directory, as an array of paths.
	 * These will be joined with either slashes (for the URL) or the file
	 * separator (for directories). The returned results should be constant,
	 * to allow for caching.
	 * @return A list of directories, relative to the /pub/pdb directory on the server
	 */
	protected abstract String[] getSplitDirPath();
	/**
	 * Location of obsolete files within the directory, as an array of paths.
	 * These will be joined with either slashes (for the URL) or the file
	 * separator (for directories). The returned results should be constant,
	 * to allow for caching.
	 * @return A list of directories, relative to the /pub/pdb directory on the server
	 */
	protected abstract String[] getObsoleteDirPath();
}
