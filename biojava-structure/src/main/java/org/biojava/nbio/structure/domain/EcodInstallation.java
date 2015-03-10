package org.biojava.nbio.structure.domain;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EcodInstallation {
	private static final Logger logger = LoggerFactory.getLogger(EcodInstallation.class);

	public static final String DEFAULT_VERSION = "latest";
	private static final String DOMAINS_FILENAME_FORMAT = "ecod.%s.domains.txt";

	public static final String ECOD_URL = "http://prodata.swmed.edu";
	public static final String DOMAINS_PATH = "/ecod/distributions/";

	public static final Pattern ECOD_RE = Pattern.compile("e(....)(.)(.)");


	private String cacheLocation;
	private String version;

	// lock to prevent multiple threads from downloading simultaneously
	// Should hold the lock when reading/writing allDomains or domainMap
	private ReadWriteLock domainsFileLock;
	private List<EcodDomain> allDomains;
	private Map<String,List<EcodDomain>> domainMap;//PDB ID -> domains

	private String url;

	public EcodInstallation(String cacheLocation) {
		domainsFileLock = new ReentrantReadWriteLock();

		this.cacheLocation = cacheLocation;

		this.version = DEFAULT_VERSION;
		this.url = ECOD_URL;

		allDomains = null; // null signals it needs to be parsed
		domainMap = null; // null signals it needs to be constructed from allDomains
	}

	public EcodInstallation() {
		this( new UserConfiguration().getCacheFilePath() );
	}

	/**
	 * Get a list of all ECOD domains for a particular PDB ID
	 * @param pdbId
	 * @return the list of domains, or null if no matching domains were found
	 * @throws IOException
	 */
	public List<EcodDomain> getDomainsForPDB(String pdbId) throws IOException {
		domainsFileLock.readLock().lock();
		try {
			logger.trace("LOCK readlock");
			while( domainMap == null ) {
				// unlock to allow ensureDomainsFileInstalled to get the write lock
				logger.trace("UNLOCK readlock");
				domainsFileLock.readLock().unlock();
				indexDomains();
				domainsFileLock.readLock().lock();
				logger.trace("LOCK readlock");
			}

			if(pdbId != null)
				pdbId = pdbId.toLowerCase();
			List<EcodDomain> doms = domainMap.get(pdbId);
			if(doms == null) {
				return null;
			}
			// Deep clone
			List<EcodDomain> clonedDoms = new ArrayList<EcodDomain>(doms.size());
			for(EcodDomain d : doms) {
				clonedDoms.add( new EcodDomain(d) );
			}
			return clonedDoms;
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}
	}

	public EcodDomain getDomainsById(String ecodId) throws IOException {
		if(ecodId == null || ecodId.isEmpty()) {
			return null;
		}

		Matcher match = ECOD_RE.matcher(ecodId);
		String pdbId = null;
		if( match.matches() )
			pdbId = match.group(1);
		List<EcodDomain> doms = getDomainsForPDB(pdbId);
		if(doms == null) {
			logger.debug("Null domains for {} from {}",pdbId,ecodId);
			return null;
		}
		logger.debug("Got {} domains from {}",doms.size(),pdbId);
		for(EcodDomain d: doms) {
			if(ecodId.equals(d.getDomainId())) {
				return d;
			}
		}
		return null;
	}

	public List<EcodDomain> getAllDomains() throws IOException {
		domainsFileLock.readLock().lock();
		logger.trace("LOCK readlock");
		try {
			while( allDomains == null) {
				// unlock to allow ensureDomainsFileInstalled to get the write lock
				logger.trace("UNLOCK readlock");
				domainsFileLock.readLock().unlock();
				ensureDomainsFileInstalled();
				domainsFileLock.readLock().lock();
				logger.trace("LOCK readlock");
			}
			return allDomains;
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}

	}

	/**
	 * Clears all domains, requiring the file to be reparsed for subsequent accesses
	 */
	public void clear() {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		allDomains = null;
		domainMap = null;
		logger.trace("UNLOCK writelock");
		domainsFileLock.writeLock().unlock();
	}
	public String getVersion() {
		return version;
	}
	public void setVersion(String version) {
		domainsFileLock.readLock().lock();
		logger.trace("LOCK readlock");
		try {
			if(version.equals(this.version)) {
				return; //no change
			}
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}

		// update version and force reparsing
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			this.version = version;
			this.clear();
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
	}

	/**
	 * Get the top-level ECOD server URL. Defaults to "http://prodata.swmed.edu"
	 * @return the url to the ecod server
	 */
	public String getUrl() {
		return url;
	}

	/**
	 * Specify a different mirror for the ECOD server.
	 * @param urlFormat the urlFormat to set
	 */
	public void setUrl(String url) {
		this.url = url;
	}

	/**
	 * Get the location of the cache directory (usually set to the PDB_CACHE_DIR
	 * property). ECOD files will be downloaded to this directory
	 * @return
	 */
	public String getCacheLocation() {
		return cacheLocation;
	}
	/**
	 * Set an alternate download location for files
	 * @param cacheLocation
	 */
	public void setCacheLocation(String cacheLocation) {
		if(cacheLocation.equals(this.cacheLocation)) {
			return; //no change
		}
		// update location
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		this.cacheLocation = cacheLocation;
		logger.trace("UNLOCK writelock");
		domainsFileLock.writeLock().unlock();
	}

	/**
	 * Blocks until ECOD domains file has been downloaded and parsed.
	 * 
	 * This may be useful in multithreaded environments
	 * @throws IOException
	 */
	public void ensureDomainsFileInstalled() throws IOException{
		// Quick check for availability
		domainsFileLock.readLock().lock();
		logger.trace("LOCK readlock");
		try {
			if( allDomains != null ) {
				return;
			}
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}

		// Download domains
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			if( !domainsAvailable() ) {
				downloadDomains();
			}
			parseDomains();
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
	}

	/**
	 * Checks that the domains file has been downloaded
	 * @return
	 */
	private boolean domainsAvailable() {
		domainsFileLock.readLock().lock();
		logger.trace("LOCK readlock");
		try {
			File f = getDomainFile();

			return f.exists() && f.length()>0;
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}
	}

	/**
	 * Downloads the domains file, overwriting any existing file
	 * @throws IOException
	 */
	private void downloadDomains() throws IOException {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			URL domainsURL = new URL( url + DOMAINS_PATH + getDomainFilename());
			File localFile = getDomainFile();

			logger.info("Downloading {} to: {}",domainsURL, localFile);
			FileDownloadUtils.downloadFile(domainsURL, localFile);
		} catch (MalformedURLException e) {
			logger.error("Malformed url: "+ url + DOMAINS_PATH + getDomainFilename(),e);
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
	}

	/**
	 * Basename for the domains file with the current version.
	 * @return
	 */
	private String getDomainFilename() {
		return  String.format(DOMAINS_FILENAME_FORMAT,getVersion());
	}

	/**
	 * Local location for the domain file
	 * @return
	 */
	private File getDomainFile() {
		return new File(getCacheLocation(),getDomainFilename());
	}

	/**
	 * Parses the domains from the local file
	 * @throws IOException
	 */
	private void parseDomains() throws IOException {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			EcodParser parser = new EcodParser(getDomainFile());
			allDomains = parser.getDomains();
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
	}

	/**
	 * Populates domainMap
	 * @throws IOException 
	 */
	private void indexDomains() throws IOException {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			if( allDomains == null) {
				ensureDomainsFileInstalled();
			}

			// Leave enough space for all PDBs as of 2015
			domainMap = new HashMap<String, List<EcodDomain>>((int) (150000/.85),.85f);

			// Index with domainMap
			for(EcodDomain d : allDomains) {
				// Get the PDB ID, either directly or from the domain ID
				String pdbId = d.getPdbId();
				if( pdbId == null ) {
					String ecodId = d.getDomainId();
					if( ecodId != null && !ecodId.isEmpty() ) {
						Matcher match = ECOD_RE.matcher(ecodId);
						pdbId = match.group(1);
					}
				}

				// Add current domain to the map
				List<EcodDomain> currDomains;
				if( domainMap.containsKey(pdbId) ) {
					currDomains = domainMap.get(pdbId);
				} else {
					currDomains = new LinkedList<EcodDomain>();
					domainMap.put(pdbId,currDomains);
				}
				currDomains.add(d);
			}
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}

	}


	public static class EcodParser {

		private final List<EcodDomain> domains;

		public EcodParser(String filename) throws IOException {
			this(new File(filename));
		}
		public EcodParser(File file) throws IOException {
			this(new FileReader(file));
		}
		public EcodParser(Reader reader) throws IOException {
			this(new BufferedReader(reader));
		}
		public EcodParser(BufferedReader reader) throws IOException {
			domains =Collections.unmodifiableList( parse(reader) );
		}

		private static List<EcodDomain> parse(BufferedReader in) throws IOException {
			ArrayList<EcodDomain> domains = null;
			try {
				// Allocate plenty of space for ECOD as of 2015 
				domains = new ArrayList<EcodDomain>(500000);

				String line = in.readLine();
				int lineNum = 0;
				while( line != null ) {
					// Ignore comments
					if( line.charAt(0) != '#' ) {
						String[] fields = line.split("\t");
						if( fields.length == 13 ) {
							try {
								int i = 0; // field number, to allow future insertion of fields

								Long uid = Long.parseLong(fields[i++]);
								String domainId = fields[i++];
								Boolean manual = null;

								// heirarchical field, e.g. "1.1.4"
								String[] xhtGroup = fields[i++].split("\\.");
								if(xhtGroup.length != 3) {
									logger.warn("Unexpected format for heirarchical field \"{}\" in line {}",fields[i-1],lineNum);
								}
								Integer xGroup = xhtGroup.length>0 ? Integer.parseInt(xhtGroup[0]) : null;
								Integer hGroup = xhtGroup.length>1 ? Integer.parseInt(xhtGroup[1]) : null;
								Integer tGroup = xhtGroup.length>2 ? Integer.parseInt(xhtGroup[2]) : null;

								String pdbId = fields[i++];
								String chainId = fields[i++];
								String range = fields[i++];

								// Intern strings likely to be shared by many domains
								String architectureName = fields[i++].intern();
								String xGroupName = fields[i++].intern();
								String hGroupName = fields[i++].intern();
								String tGroupName = fields[i++].intern();
								String fGroupName = fields[i++].intern();

								Boolean isAssembly = null;
								String assemblyStr = fields[i++];
								if(assemblyStr.equals("NOT_DOMAIN_ASSEMBLY")) {
									isAssembly = false;
								} else if(assemblyStr.equals("IS_DOMAIN_ASSEMBLY")) {
									isAssembly = true;
								} else {
									logger.warn("Unexpected value for assembly field \"{}\" in line {}",assemblyStr,lineNum);
								}

								String ligandStr = fields[i++];
								Set<String> ligands = null;
								if( ligandStr.equals("NO_LIGANDS_4A") || ligandStr.isEmpty() ) {
									ligands = Collections.emptySet();
								} else {
									String[] ligSplit = ligandStr.split(",");
									ligands = new LinkedHashSet<String>(ligSplit.length);
									for(String s : ligSplit) {
										ligands.add(s.intern());
									}
								}


								EcodDomain domain = new EcodDomain(uid, domainId, manual, xGroup, hGroup, tGroup, pdbId, chainId, range, architectureName, xGroupName, hGroupName, tGroupName, fGroupName, isAssembly, ligands);
								domains.add(domain);
							} catch(NumberFormatException e) {
								logger.warn("Error in ECOD parsing at line "+lineNum,e);
							}
						} else {
							logger.warn("Unexpected number of fields in line {}",lineNum);
						}
					}

					line = in.readLine();
					lineNum++;
				}
				logger.info("Parsed {} ECOD domains",domains.size());
			} finally {
				if(in != null) {
					in.close();
				}
			}

			return domains;
		}

		/**
		 * @return a list of all EcodDomains
		 */
		public List<EcodDomain> getDomains() {
			return domains;
		}
	}


	public static void main(String[] args) {
		if( args.length!= 1) {
			System.out.println("usage: ecod_domains.txt");
			System.exit(1); return;
		}

		String filename = args[0];

		try {
			EcodParser parser = new EcodParser(filename);

			List<EcodDomain> domains = parser.getDomains();

			System.out.format("Found %d ECOD domains.%n",domains.size());

			System.out.println("First 10 domains:");
			int i = 0;
			for(EcodDomain d: domains) {
				if( i>10) break;

				System.out.println(d.getDomainId());
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
