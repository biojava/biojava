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
import java.util.Calendar;
import java.util.Collections;
import java.util.Date;
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

import org.biojava.nbio.structure.PdbId;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.core.util.FileDownloadUtils;
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
	private Map<PdbId,List<EcodDomain>> domainMap;//PDB ID -> domains, lazily constructed from allDomains

	private String url;

	// Frequency of ECOD updates, in days. If non-null, redownloads "latest" if older than this.
	private Integer updateFrequency = 14;

	/**
	 * Use EcodFactory to create instances. The instantiation of multiple
	 * installations at the same path can lead to race conditions when downloading
	 * files.
	 * @param cacheLocation Location to save files, typically from the PDB_CACHE_DIR parameter
	 * @param version ECOD requestedVersion to fetch
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
	public EcodInstallation() {
		this( new UserConfiguration().getCacheFilePath(), DEFAULT_VERSION );
	}
	/**
	public EcodInstallation(String cacheLocation) {
		this( cacheLocation, DEFAULT_VERSION );
	}

	/**
	 * Get a list of all ECOD domains for a particular PDB ID
	 * @param id
	 * @return the list of domains, or null if no matching domains were found
	 * @throws IOException
	 */
	@Override
	public List<EcodDomain> getDomainsForPdb(String id) throws IOException {
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

			PdbId pdbId = null;
			try {
				pdbId = new PdbId(id);
			} catch (IllegalArgumentException e) {
				return null;
			}
			List<EcodDomain> doms = domainMap.get(pdbId);
			if(doms == null) {
				return null;
			}
			// Deep clone
			List<EcodDomain> clonedDoms = new ArrayList<>(doms.size());
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
	 * Get a list of domains within a particular level of the hierarchy
	 * @param hierarchy A dot-separated list giving the X-group, H-group, and/or
	 *  T-group (e.g. "1.1" for all members of the RIFT-related H-group)
	 * @return
	 * @throws IOException
	 */
	@Override
	public List<EcodDomain> filterByHierarchy(String hierarchy) throws IOException {
		String[] xhtGroup = hierarchy.split("\\.");
		Integer xGroup = xhtGroup.length>0 ? Integer.parseInt(xhtGroup[0]) : null;
		Integer hGroup = xhtGroup.length>1 ? Integer.parseInt(xhtGroup[1]) : null;
		Integer tGroup = xhtGroup.length>2 ? Integer.parseInt(xhtGroup[2]) : null;

		List<EcodDomain> filtered = new ArrayList<>();
		for(EcodDomain d: getAllDomains()) {
			boolean match = true;
			if(xhtGroup.length>0) {
				match = match && xGroup.equals(d.getXGroup());
			}
			if(xhtGroup.length>1) {
				match = match && hGroup.equals(d.getHGroup());
			}
			if(xhtGroup.length>2) {
				match = match && tGroup.equals(d.getTGroup());
			}
			if(xhtGroup.length>3) {
				logger.warn("Ignoring unexpected additional parts of ECOD {}",hierarchy);
			}
			if(match) {
				filtered.add(d);
			}
		}
		return filtered;
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
	 * <p>
	 * Since 7.3.0 this reads only the file's header rather than parsing the whole
	 * file, so it no longer has the side effect of loading every domain.
	 * @return the ECOD version
	 * @throws IOException If an error occurs while downloading or parsing the file
	 */
	@Override
	public String getVersion() throws IOException {
		domainsFileLock.readLock().lock();
		logger.trace("LOCK readlock");
		try {
			if( parsedVersion != null ) {
				return parsedVersion;
			}
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}

		// The version is declared in the first few lines of the file, so read those rather
		// than the millions of domain records behind them. The current release is 657 MB and
		// holds nearly three million records; parsing it in full to answer this question
		// costs over a gigabyte of heap and several seconds.
		ensureDomainsFileDownloaded();

		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			if( parsedVersion == null ) {
				parsedVersion = parseVersionOnly();
			}
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}

		if( parsedVersion == null) {
			return requestedVersion;
		}
		return parsedVersion;
	}

	/**
	 * Reads the version from the header of the local domains file without parsing the
	 * domains themselves.
	 * @return the version, or null if the header does not declare one
	 * @throws IOException if the file cannot be read
	 * @since 7.3.0
	 */
	private String parseVersionOnly() throws IOException {
		try( BufferedReader in = new BufferedReader(new FileReader(getDomainFile())) ) {
			String line;
			while( (line = in.readLine()) != null ) {
				Matcher match = EcodParser.VERSION_RE.matcher(line);
				if( match.matches() ) {
					return match.group(1);
				}
				if( !line.startsWith("#") ) {
					// past the header block; from v294.1 the column names are not commented
					return null;
				}
			}
		}
		return null;
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
	 * @param url the urlFormat to set
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
	 * Ensures the domains file is present and current locally, without parsing it.
	 * @throws IOException in cases of file I/O, including failure to download a healthy file
	 * @since 7.3.0
	 */
	private void ensureDomainsFileDownloaded() throws IOException {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			if( !domainsAvailable() ) {
				downloadDomains();
			}
		} finally {
			logger.trace("UNLOCK writelock");
			domainsFileLock.writeLock().unlock();
		}
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

			if (! (f.exists() && FileDownloadUtils.validateFile(f)))
				return false;

			// Re-download old copies of "latest"
			if(updateFrequency != null && requestedVersion.equals(DEFAULT_VERSION)) {
				long mod = f.lastModified();
				// Time of last update
				Date lastUpdate = new Date();
				Calendar cal = Calendar.getInstance();
				cal.setTime(lastUpdate);
				cal.add(Calendar.DAY_OF_WEEK, -updateFrequency);
				long updateTime = cal.getTimeInMillis();
				// Check if file predates last update
				if( mod < updateTime ) {
					logger.info("{} is out of date.",f);
					return false;
				}
			}
			return true;
		} finally {
			logger.trace("UNLOCK readlock");
			domainsFileLock.readLock().unlock();
		}
	}

	/**
	 * Downloads the domains file +/- its validation metadata, overwriting any existing file
	 * @throws IOException in cases of file I/O, including failure to download a healthy (non-corrupted) file.
	 */
	private void downloadDomains() throws IOException {
		domainsFileLock.writeLock().lock();
		logger.trace("LOCK writelock");
		try {
			URL domainsURL = new URL( url + DOMAINS_PATH + getDomainFilename());
			File localFile = getDomainFile();

			logger.info("Downloading {} to: {}",domainsURL, localFile);
			FileDownloadUtils.createValidationFiles(domainsURL, localFile, null, FileDownloadUtils.Hash.UNKNOWN);
			FileDownloadUtils.downloadFile(domainsURL, localFile);
			if(! FileDownloadUtils.validateFile(localFile))
				throw new IOException("Downloaded file invalid: "+ localFile);
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
		return  String.format(DOMAINS_FILENAME_FORMAT,requestedVersion);
	}

	/**
	 * Local location for the domain file
	 * @return
	 */
	private File getDomainFile() {
		return new File(getCacheLocation(),getDomainFilename());
	}

	/**
	 * The expected ECOD update frequency determines whether the version
	 * "latest" should be re-downloaded
	 * @return the expected ECOD update frequency, in days
	 */
	public Integer getUpdateFrequency() {
		return updateFrequency;
	}

	/**
	 * The "latest" version will be re-downloaded if it is older than
	 * {@link #getUpdateFrequency()} days. Setting this to null disables
	 * re-downloading (delete $PDB_CACHE_DIR/ecod.latest.domains.txt manually
	 * to force updating). Setting to 0 will force downloading for every
	 * program execution.
	 * @param updateFrequency the updateFrequency to set
	 */
	public void setUpdateFrequency(Integer updateFrequency) {
		this.updateFrequency = updateFrequency;
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
			domainMap = new HashMap<>((int) (150000/.85),.85f);

			// Index with domainMap
			for(EcodDomain d : allDomains) {
				// Get the PDB ID, either directly or from the domain ID
				PdbId pdbId = d.getPdbId();
				if( pdbId == null ) {
					String ecodId = d.getDomainId();
					if( ecodId != null && !ecodId.isEmpty() ) {
						Matcher match = ECOD_RE.matcher(ecodId);
						pdbId = new PdbId(match.group(1));
					}
				}

				// Add current domain to the map
				List<EcodDomain> currDomains;
				if( domainMap.containsKey(pdbId) ) {
					currDomains = domainMap.get(pdbId);
				} else {
					currDomains = new LinkedList<>();
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
		/*
Version Notes

Current version (1.4) contains the following columns:

Column 1: ECOD uid - internal domain unique identifier
Column 2: ECOD domain id - domain identifier
Column 3: ECOD representative status - manual (curated) or automated nonrep
Column 4: ECOD hierachy identifier - [X-group].[H-group].[T-group].[F-group]
	* In develop45-66 these also include single numbers in the range 1-265
Column 5: PDB identifier
Column 6: Chain identifier (note: case-sensitive)
Column 7: PDB residue number range
	* These are sometimes incorrect up to at least develop124. Examples are:
	  e4lxaA2 (should be A:184-385), e4lxmC3 (should be C:46P-183)
Column 8: seq_id number range (based on internal PDB indices)
Column 9: Architecture name
Column 10: X-group name
Column 11: H-group name
Column 12: T-group name
Column 13: F-group name (F_UNCLASSIFIED denotes that domain has not been assigned to an F-group)
Column 14: Domain assembly status (if domain is member of assembly, partners' ecod domain ids listed)
Column 15: Comma-separated value list of non-polymer entities within 4 A of at least one residue of domain

Notes older versions:
changelog:
v1.0 - original version (8/04/2014)
v1.1 - added rep/nonrep data (1/15/2015)
v1.2 - added f-group identifiers to fasta file, domain description file. ECODf identifiers now used when available for F-group name.
	Domain assemblies now represented by assembly uid in domain assembly status.
v1.4 - added seqid_range and headers (develop101)
v1.6 - renamed column 4 from f_id to t_id and inserted unp_acc (UniProt accession) as
	column 9, giving 16 columns (seen in develop291)

From v294.1 the distribution was redesigned. The header comment changed from
"#ECOD version develop291" to "# Version: v294.1", the column header row is no longer
commented out, and the columns became:

	uid ecod_domain_id manual_rep f_id pdb chain pdb_range seqid_range architecture_name
	x_name h_name t_name f_name assembly_id domain_id_short range_count arch_manual
	x_manual h_manual t_manual f_manual valid_structure ligand_binding

v295 appends ligand_comp_ids and ligand_pdbnum, for 25 columns. Also note that
manual_rep now holds True/False rather than MANUAL_REP/AUTO_NONREP, that assembly_id
and domain_id_short are empty on every row, that f_name is empty rather than
F_UNCLASSIFIED for unclassified domains, and that uid restarts from 0.

Because the columns have been renamed, reordered and added to repeatedly, files that
declare a column header are read by column name rather than by position.
		 */

		/** String for unclassified F-groups */
		public static final String F_UNCLASSIFIED = "F_UNCLASSIFIED";
		/** String for single-domain assemblies */
		public static final String NOT_DOMAIN_ASSEMBLY = "NOT_DOMAIN_ASSEMBLY";
		/** Deprecated way of indicating there is an assembly. replaced by the assembly id */
		public static final String IS_DOMAIN_ASSEMBLY = "IS_DOMAIN_ASSEMBLY";
		/** Indicates a manual representative */
		public static final String IS_REPRESENTATIVE = "MANUAL_REP";
		/** Indicates not a manual representative */
		public static final String NOT_REPRESENTATIVE = "AUTO_NONREP";
		/**
		 * Matches the comment declaring the version, which has taken two forms:
		 * {@code #ECOD version develop291} up to develop292, and {@code # Version: v295}
		 * from v294.1 onwards.
		 * @since 7.3.0
		 */
		static final Pattern VERSION_RE = Pattern.compile(
				"^\\s*#\\s*(?:ECOD\\s+)?version\\s*:?\\s*(\\S+).*", Pattern.CASE_INSENSITIVE);

		private List<EcodDomain> domains;
		private String version;

		// prevent too many warnings; negative numbers print all warnings
		private int warnIsDomainAssembly = 1;
		private int warnHierarchicalFormat = 5;
		private int warnNumberOfFields = 10;
		private int warnNumberFormat = 10;
		/** Data lines that could not be turned into a domain, for the summary at the end */
		private int skippedLines = 0;
		/** Data lines describing a domain in a computed model rather than a PDB entry */
		private int modelLines = 0;

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
				ArrayList<EcodDomain> domainsList = new ArrayList<>(500000);

				Pattern commentRE = Pattern.compile("^\\s*#.*");

				ColumnLayout layout = null;

				String line = in.readLine();
				int lineNum = 1;
				while( line != null ) {
					// Check for requestedVersion string
					Matcher match = VERSION_RE.matcher(line);
					if(match.matches()) {
						// special requestedVersion comment
						this.version = match.group(1);
					} else if( ColumnLayout.isColumnHeader(line) ) {
						// The column names. Since the columns have been renamed, reordered and
						// added to several times, later lines are read by name rather than by
						// position wherever this header is present (develop101 onwards).
						layout = ColumnLayout.fromHeader(line);
						logger.debug("Read ECOD column header at line {}: {} columns",lineNum,layout.size());
					} else {
						match = commentRE.matcher(line);
						if(match.matches()) {
							// ignore comments
						} else {
							// data line. The last column is frequently empty, so keep trailing
							// empty fields rather than letting split() discard them.
							String[] fields = line.split("\t", -1);
							if( layout != null ) {
								String pdb = layout.get(fields, "pdb");
								if( pdb != null && pdb.isEmpty() ) {
									// From v294.1 the distribution also classifies domains
									// found in computed (AlphaFold) models, which have no PDB
									// entry and so cannot be represented by an EcodDomain.
									modelLines++;
								} else {
									try {
										EcodDomain domain = parseDomain(fields, layout, lineNum);
										if(domain != null) {
											domainsList.add(domain);
										} else {
											skippedLines++;
											warnMissingColumns(lineNum);
										}
									} catch(IllegalArgumentException e) {
										// includes NumberFormatException and an unusable PDB id
										skippedLines++;
										warnUnparseableLine(lineNum, e);
									}
								}
							} else if( fields.length == 13 || fields.length == 14 || fields.length == 15) {
								try {
									int i = 0; // field number, to allow future insertion of fields

									//Column 1: ECOD uid - internal domain unique identifier
									Long uid = Long.parseLong(fields[i++]);
									//Column 2: ECOD domain id - domain identifier
									String domainId = fields[i++];

									//Column 3: ECOD representative status - manual (curated) or automated nonrep
									// Manual column may be missing in version 1.0 files
									Boolean manual = null;
									if( fields.length >= 14) {
										manual = parseManualRep(fields[i++], lineNum);
									}

									//Column 4: ECOD hierachy identifier - [X-group].[H-group].[T-group].[F-group]
									// hierarchical field, e.g. "1.1.4.1"
									Integer[] xhtfGroup = parseHierarchy(fields[i++], lineNum);
									Integer xGroup = xhtfGroup[0];
									Integer hGroup = xhtfGroup[1];
									Integer tGroup = xhtfGroup[2];
									Integer fGroup = xhtfGroup[3];

									//Column 5: PDB identifier
									String pdbId = fields[i++];
									//Column 6: Chain identifier (note: case-sensitive)
									String chainId = fields[i++];
									//Column 7: PDB residue number range
									String range = fields[i++];

									//Column 8: seq_id number range (based on internal PDB indices)
									//Added in version 1.4
									String seqId = null;
									if( fields.length >= 15) {
										seqId = fields[i++];
									}

									//Column 9: Architecture name
									// Intern strings likely to be shared by many domains
									String architectureName = fields[i++].intern();
									//Column 10: X-group name
									String xGroupName = fields[i++].intern();
									//Column 11: H-group name
									String hGroupName = fields[i++].intern();
									//Column 12: T-group name
									String tGroupName = fields[i++].intern();
									//Column 13: F-group name (F_UNCLASSIFIED denotes that domain has not been assigned to an F-group)
									//Contents changed in version 1.3
									String fGroupName = fields[i++].intern();


									hGroupName = clearStringQuotes(hGroupName);
									tGroupName = clearStringQuotes(tGroupName);
									fGroupName = clearStringQuotes(fGroupName);
									xGroupName = clearStringQuotes(xGroupName);

									//Column 14: Domain assembly status (if domain is member of assembly, partners' ecod domain ids listed)
									//Column 15: Comma-separated value list of non-polymer entities within 4 A of at least one residue of domain
									Long assemblyId = null;
									String assemblyStr = fields[i++];
									if(assemblyStr.equals(NOT_DOMAIN_ASSEMBLY)) {
										assemblyId = uid;
									} else if("IS_DOMAIN_ASSEMBLY".equals(assemblyStr) ) {
										if(warnIsDomainAssembly > 1) {
											logger.info("Deprecated 'IS_DOMAIN_ASSEMBLY' value ignored in line {}.",lineNum);
											warnIsDomainAssembly--;
										} else if(warnIsDomainAssembly == 0) {
											logger.info("Deprecated 'IS_DOMAIN_ASSEMBLY' value ignored in line {}. Not printing future similar warnings.",lineNum);
											warnIsDomainAssembly--;
										}
										//assemblyId = null;
									} else {
										assemblyId = Long.parseLong(assemblyStr);
									}

									Set<String> ligands = parseLigands(fields[i++]);


									EcodDomain domain = new EcodDomain(uid, domainId, manual, xGroup, hGroup, tGroup, fGroup,pdbId, chainId, range, seqId, architectureName, xGroupName, hGroupName, tGroupName, fGroupName, assemblyId, ligands);
									domainsList.add(domain);
								} catch(NumberFormatException e) {
									skippedLines++;
									warnUnparseableLine(lineNum, e);
								}
							} else {
								skippedLines++;
								warnMissingColumns(lineNum);
							}
						}
					}

					line = in.readLine();
					lineNum++;
				}
				if(this.version == null)
					logger.info("Parsed {} ECOD domains",domainsList.size());
				else
					logger.info("Parsed {} ECOD domains from version {}",domainsList.size(),this.version);

				if(modelLines > 0) {
					logger.info("Ignored {} ECOD domains classified from computed models, "
							+ "which have no PDB entry", modelLines);
				}

				if(domainsList.isEmpty() && skippedLines > 0) {
					// Returning an empty list quietly is how an upstream format change went
					// unnoticed for eight months. Say so instead.
					logger.error("Parsed no ECOD domains from {} data lines of version {}. "
							+ "The file format has probably changed; please report this at "
							+ "https://github.com/biojava/biojava/issues", skippedLines,
							this.version == null ? "unknown" : this.version);
				} else if(skippedLines > 0) {
					logger.warn("Skipped {} of {} ECOD data lines that could not be parsed",
							skippedLines, skippedLines + domainsList.size());
				}


				this.domains = Collections.unmodifiableList( domainsList );

			} finally {
				if(in != null) {
					in.close();
				}
			}
		}

		/**
		 * Builds a domain from a data line using the column names the file declares in its
		 * header, rather than fixed offsets. This is what allows one parser to read the
		 * 15-column develop101 layout, the 16-column develop291 layout (which inserts
		 * {@code unp_acc}) and the 23- and 25-column v294.1 and v295 layouts.
		 * @param fields the tab-separated values of one data line
		 * @param layout the column names read from the file's header
		 * @param lineNum the line number, for warnings
		 * @return the domain, or null if the line does not carry every required column
		 * @throws NumberFormatException if a numeric column does not hold a number
		 * @since 7.3.0
		 */
		private EcodDomain parseDomain(String[] fields, ColumnLayout layout, int lineNum) {
			String uidStr    = layout.get(fields, "uid");
			String domainId  = layout.get(fields, "ecod_domain_id");
			// renamed from t_id to f_id when the hierarchy gained a fourth level
			String hierarchy = layout.get(fields, "f_id", "t_id");
			String pdbId     = layout.get(fields, "pdb");
			String chainId   = layout.get(fields, "chain");
			String range     = layout.get(fields, "pdb_range");
			if( uidStr == null || domainId == null || hierarchy == null
					|| pdbId == null || chainId == null || range == null ) {
				return null;
			}

			Long uid = Long.parseLong(uidStr);
			Boolean manual = parseManualRep(layout.get(fields, "manual_rep"), lineNum);
			Integer[] xhtfGroup = parseHierarchy(hierarchy, lineNum);
			// absent before version 1.4
			String seqId = layout.get(fields, "seqid_range");

			String architectureName = internName(layout.get(fields, "architecture_name", "arch_name"));
			String xGroupName = internName(layout.get(fields, "x_name"));
			String hGroupName = internName(layout.get(fields, "h_name"));
			String tGroupName = internName(layout.get(fields, "t_name"));
			// Up to develop292 an unclassified domain carried F_UNCLASSIFIED here. From
			// v294.1 the name is simply empty, while f_id still classifies the domain to
			// four levels, so the two are no longer equivalent and the empty value is
			// deliberately left as it is rather than translated.
			String fGroupName = internName(layout.get(fields, "f_name"));

			// v294.1 and later declare assembly_id but leave it empty on every row, which
			// means the same as the NOT_DOMAIN_ASSEMBLY of earlier versions.
			Long assemblyId = null;
			String assemblyStr = layout.get(fields, "assembly_id", "asm_status");
			if( assemblyStr == null || assemblyStr.isEmpty() || NOT_DOMAIN_ASSEMBLY.equals(assemblyStr) ) {
				assemblyId = uid;
			} else if( IS_DOMAIN_ASSEMBLY.equals(assemblyStr) ) {
				warnDomainAssembly(lineNum);
			} else {
				assemblyId = Long.parseLong(assemblyStr);
			}

			// the ligand list moved from the last column to ligand_comp_ids in v295
			Set<String> ligands = parseLigands(layout.get(fields, "ligand_comp_ids", "ligand"));

			return new EcodDomain(uid, domainId, manual, xhtfGroup[0], xhtfGroup[1], xhtfGroup[2],
					xhtfGroup[3], pdbId, chainId, range, seqId, architectureName, xGroupName,
					hGroupName, tGroupName, fGroupName, assemblyId, ligands);
		}

		/**
		 * Reads the representative-status column, which held MANUAL_REP or AUTO_NONREP up to
		 * develop292 and True or False from v294.1 onwards.
		 * @return true, false, or null if the column is absent or unrecognised
		 * @since 7.3.0
		 */
		private Boolean parseManualRep(String manualString, int lineNum) {
			if(manualString == null) {
				return null;
			}
			if(manualString.equalsIgnoreCase(IS_REPRESENTATIVE) || manualString.equalsIgnoreCase("true")) {
				return true;
			}
			if(manualString.equalsIgnoreCase(NOT_REPRESENTATIVE) || manualString.equalsIgnoreCase("false")) {
				return false;
			}
			logger.warn("Unexpected value for manual field: {} in line {}",manualString,lineNum);
			return null;
		}

		/**
		 * Splits the hierarchical identifier, e.g. "1.1.4.1".
		 * @return the X, H, T and F group numbers, any of which may be null if absent
		 * @since 7.3.0
		 */
		private Integer[] parseHierarchy(String hierarchy, int lineNum) {
			String[] xhtGroup = hierarchy.split("\\.");
			if(xhtGroup.length < 3 || 4 < xhtGroup.length) {
				if(warnHierarchicalFormat > 1) {
					logger.warn("Unexpected format for hierarchical field \"{}\" in line {}",hierarchy,lineNum);
					warnHierarchicalFormat--;
				} else if(warnHierarchicalFormat != 0) {
					logger.warn("Unexpected format for hierarchical field \"{}\" in line {}. Not printing future similar warnings.",hierarchy,lineNum);
					warnHierarchicalFormat--;
				}
			}
			Integer[] groups = new Integer[4];
			for(int j = 0; j < groups.length && j < xhtGroup.length; j++) {
				groups[j] = Integer.parseInt(xhtGroup[j]);
			}
			return groups;
		}

		/**
		 * Reads a comma-separated list of non-polymer entities close to the domain.
		 * @return the ligands, or an empty set for NO_LIGANDS_4A, an empty value or no column
		 * @since 7.3.0
		 */
		private Set<String> parseLigands(String ligandStr) {
			if( ligandStr == null || ligandStr.isEmpty() || "NO_LIGANDS_4A".equals(ligandStr) ) {
				return Collections.emptySet();
			}
			String[] ligSplit = ligandStr.split(",");
			Set<String> ligands = new LinkedHashSet<>(ligSplit.length);
			for(String s : ligSplit) {
				ligands.add(s.intern());
			}
			return ligands;
		}

		/**
		 * Interns a name likely to be shared by many domains, stripping the quotes that
		 * versions up to develop292 wrapped some of them in.
		 * @since 7.3.0
		 */
		private String internName(String name) {
			if(name == null) {
				return null;
			}
			return clearStringQuotes(name).intern();
		}

		private void warnDomainAssembly(int lineNum) {
			if(warnIsDomainAssembly > 1) {
				logger.info("Deprecated 'IS_DOMAIN_ASSEMBLY' value ignored in line {}.",lineNum);
				warnIsDomainAssembly--;
			} else if(warnIsDomainAssembly == 0) {
				logger.info("Deprecated 'IS_DOMAIN_ASSEMBLY' value ignored in line {}. Not printing future similar warnings.",lineNum);
				warnIsDomainAssembly--;
			}
		}

		private void warnMissingColumns(int lineNum) {
			if(warnNumberOfFields > 1) {
				logger.warn("Unexpected number of fields in line {}.",lineNum);
				warnNumberOfFields--;
			} else if(warnNumberOfFields == 1) {
				logger.warn("Unexpected number of fields in line {}. Not printing future similar warnings",lineNum);
				warnNumberOfFields--;
			}
		}

		private void warnUnparseableLine(int lineNum, IllegalArgumentException e) {
			if(warnNumberFormat > 1) {
				logger.warn("Error in ECOD parsing at line {}: {}", lineNum, e.getMessage());
				warnNumberFormat--;
			} else if(warnNumberFormat == 1) {
				logger.warn("Error in ECOD parsing at line {}: {}. Not printing future similar warnings", lineNum, e.getMessage());
				warnNumberFormat--;
			}
		}

		private String clearStringQuotes(String name) {
			if ( name.startsWith("\""))
				name = name.substring(1);

			if ( name.endsWith("\""))
				name = name.substring(0,name.length()-1);

			return name;
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

		/**
		 * Maps the column names an ECOD domain file declares in its header onto their
		 * positions, so a data line can be read by name rather than by offset.
		 * <p>
		 * Every distribution since develop101 carries such a header. It is commented
		 * (<code>#uid&lt;tab&gt;ecod_domain_id&lt;tab&gt;...</code>) up to develop292 and
		 * uncommented (<code>uid&lt;tab&gt;ecod_domain_id&lt;tab&gt;...</code>) from v294.1
		 * onwards. Because names have also been changed between versions, lookups accept
		 * aliases and any name the file does not declare simply reads as absent.
		 *
		 * @author Amr ALHOSSARY
		 * @since 7.3.0
		 */
		private static class ColumnLayout {
			private final Map<String,Integer> columns;

			private ColumnLayout(Map<String,Integer> columns) {
				this.columns = columns;
			}

			/**
			 * @return true if this line names the columns rather than holding domain data
			 */
			public static boolean isColumnHeader(String line) {
				int tab = line.indexOf('\t');
				if(tab < 0) {
					return false;
				}
				String first = line.substring(0, tab).trim();
				if(first.startsWith("#")) {
					first = first.substring(1).trim();
				}
				return first.equalsIgnoreCase("uid");
			}

			public static ColumnLayout fromHeader(String line) {
				String[] names = line.split("\t", -1);
				Map<String,Integer> columns = new HashMap<>(names.length * 2);
				for(int i = 0; i < names.length; i++) {
					String name = names[i].trim();
					if(i == 0 && name.startsWith("#")) {
						name = name.substring(1).trim();
					}
					if(!name.isEmpty()) {
						columns.put(name.toLowerCase(), i);
					}
				}
				return new ColumnLayout(columns);
			}

			/**
			 * @param fields the values of one data line
			 * @param aliases the names this column has gone by, most recent first
			 * @return the value of the first alias the file declares, or null if it declares
			 *  none of them or this line is too short to reach it
			 */
			public String get(String[] fields, String... aliases) {
				for(String alias : aliases) {
					Integer i = columns.get(alias);
					if(i != null) {
						return i < fields.length ? fields[i] : null;
					}
				}
				return null;
			}

			public int size() {
				return columns.size();
			}
		}
	}


	@Override
	public String toString() {
		String version = null;
		try {
			version = getVersion();
		} catch (IOException e) {
			// For parsing errors, use the requested version
			version = requestedVersion;
		}

		return "EcodInstallation [cacheLocation=" + cacheLocation
				+ ", version=" + version + "]";
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
