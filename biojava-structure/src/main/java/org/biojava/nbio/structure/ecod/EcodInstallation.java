/*
 * BioJava development code
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
 */

package org.biojava.nbio.structure.ecod;

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

/**
 * Provides access to the Evolutionary Classification of Protein Domains (ECOD).
 * 
 * The preferred mechanism for obtaining instances of this class is through the
 * {@link EcodFactory} class.
 * 
 * Reference:
 * H. Cheng, R. D. Schaeffer, Y. Liao, L. N. Kinch, J. Pei, S. Shi, B. H.\
 *   Kim, N. V. Grishin. (2014) ECOD: An evolutionary classification of protein
 *   domains. PLoS Comput Biol 10(12): e1003926.
 * http://prodata.swmed.edu/ecod/
 * 
 * @author Spencer Bliven
 *
 */
public class EcodInstallation implements EcodDatabase {
	private static final Logger logger = LoggerFactory.getLogger(EcodInstallation.class);

	public static final String DEFAULT_VERSION = "latest";
	private static final String DOMAINS_FILENAME_FORMAT = "ecod.%s.domains.txt";

	public static final String ECOD_URL = "http://prodata.swmed.edu";
	public static final String DOMAINS_PATH = "/ecod/distributions/";

	// ECOD identifiers are e<pdbID><chain><domain>, where chain and domain
	// Chain and domain can both be multi-letter (e.g. e2q7zA10)
	public static final Pattern ECOD_RE = Pattern.compile("^e(....).+\\d+$");


	private String cacheLocation;
	private String requestedVersion; // version requested, e.g. "latest". Used for the paths
	private String parsedVersion; // actual version parsed

	// lock to prevent multiple threads from downloading simultaneously
	// Should hold the lock when reading/writing allDomains or domainMap
	private ReadWriteLock domainsFileLock;
	private List<EcodDomain> allDomains;
	private Map<String,List<EcodDomain>> domainMap;//PDB ID -> domains, lazily constructed from allDomains

	private String url;

	/**
	 * Use EcodFactory to create instances. The instantiation of multiple
	 * installations at the same path can lead to race conditions when downloading
	 * files.
	 * @param cacheLocation Location to save files, typically from the PDB_CACHE_DIR parameter
	 * @param requestedVersion ECOD requestedVersion to fetch
	 */
	public EcodInstallation(String cacheLocation, String version) {
		domainsFileLock = new ReentrantReadWriteLock();

		this.cacheLocation = cacheLocation;

		this.requestedVersion = version;
		this.url = ECOD_URL;

		allDomains = null; // null signals it needs to be parsed
		domainMap = null; // null signals it needs to be constructed from allDomains
	}

	/**
	 * @see EcodFactory#getEcodDatabase()
	 */
	EcodInstallation() {
		this( new UserConfiguration().getCacheFilePath(), DEFAULT_VERSION );
	}
	/**
	public EcodInstallation(String cacheLocation) {
		this( cacheLocation, DEFAULT_VERSION );
	}

	/**
	 * Get a list of all ECOD domains for a particular PDB ID
	 * @param pdbId
	 * @return the list of domains, or null if no matching domains were found
	 * @throws IOException
	 */
	@Override
	public List<EcodDomain> getDomainsForPdb(String pdbId) throws IOException {
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

	/**
	 * Get a particular ECOD domain by the domain ID (e.g. "e4hhbA1")
	 * @param ecodId
	 * @return
	 * @throws IOException
	 */
	@Override
	public EcodDomain getDomainsById(String ecodId) throws IOException {
		if(ecodId == null || ecodId.isEmpty()) {
			return null;
		}

		Matcher match = ECOD_RE.matcher(ecodId);
		String pdbId = null;
		if( match.matches() )
			pdbId = match.group(1);
		List<EcodDomain> doms = getDomainsForPdb(pdbId);
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

	/**
	 * Get all ECOD domains
	 * @return
	 * @throws IOException
	 */
	@Override
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
	/**
	 * Return the ECOD version, as parsed from the file.
	 * 
	 * Note that this may differ from the version requested in the constructor
	 * for the special case of "latest"
	 * @return the ECOD version
	 */
	@Override
	public String getVersion() {
		if( parsedVersion == null) {
			return requestedVersion;
		}
		return parsedVersion;
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
	 * This may be useful in multithreaded environments.
	 * @throws IOException
	 */
	// Populates allDomains
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
	 * Basename for the domains file with the current requestedVersion.
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
			parsedVersion = parser.getVersion();
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
	}

	/**
	 * Populates domainMap from allDomains
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

		private List<EcodDomain> domains;
		private String version;

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
			version = null;
			parse(reader);
		}

		private void parse(BufferedReader in) throws IOException {
			try {
				// Allocate plenty of space for ECOD as of 2015 
				ArrayList<EcodDomain> domainsList = new ArrayList<EcodDomain>(500000);

				Pattern versionRE = Pattern.compile("^\\s*#.*ECOD\\s*requestedVersion\\s+(\\w+)");
				Pattern commentRE = Pattern.compile("^\\s*#");

				String line = in.readLine();
				int lineNum = 0;
				while( line != null ) {
					// Check for requestedVersion string
					Matcher match = versionRE.matcher(line);
					if(match.matches()) {
						// special requestedVersion comment
						this.version = match.group(1);
					} else {
						match = commentRE.matcher(line);
						if(match.matches()) {
							// ignore comments
						} else {
							// data line
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
									domainsList.add(domain);
								} catch(NumberFormatException e) {
									logger.warn("Error in ECOD parsing at line "+lineNum,e);
								}
							} else {
								logger.warn("Unexpected number of fields in line {}",lineNum);
							}
						}
					}

					line = in.readLine();
					lineNum++;
				}
				if(this.version == null)
					logger.info("Parsed {} ECOD domains",domainsList.size());
				else
					logger.info("Parsed {} ECOD domains from requestedVersion {}",domainsList.size(),this.version);


				this.domains = Collections.unmodifiableList( domainsList );

			} finally {
				if(in != null) {
					in.close();
				}
			}
		}

		/**
		 * @return a list of all EcodDomains
		 */
		public List<EcodDomain> getDomains() {
			return domains;
		}
		
		/**
		 * @return the requestedVersion for this file, or null if none was parsed
		 */
		public String getVersion() {
			return version;
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
