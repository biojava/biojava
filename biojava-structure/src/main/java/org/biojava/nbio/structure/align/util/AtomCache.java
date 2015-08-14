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
package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.cath.CathSegment;
import org.biojava.nbio.structure.domain.PDPProvider;
import org.biojava.nbio.structure.domain.RemotePDPProvider;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.io.util.FileDownloadUtils;
import org.biojava.nbio.structure.quaternary.io.BioUnitDataProviderFactory;
import org.biojava.nbio.structure.scop.*;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A utility class that provides easy access to Structure objects. If you are running a script that is frequently
 * re-using the same PDB structures, the AtomCache keeps an in-memory cache of the files for quicker access. The cache
 * is a soft-cache, this means it won't cause out of memory exceptions, but garbage collects the data if the Java
 * virtual machine needs to free up space. The AtomCache is thread-safe.
 * 
 * @author Andreas Prlic
 * @author Spencer Bliven
 * @author Peter Rose
 * @since 3.0
 */
public class AtomCache {
	
	private static final Logger logger = LoggerFactory.getLogger(AtomCache.class);

	public static final String BIOL_ASSEMBLY_IDENTIFIER = "BIO:";
	public static final String CHAIN_NR_SYMBOL = ":";
	public static final String CHAIN_SPLIT_SYMBOL = ".";

	public static final String PDP_DOMAIN_IDENTIFIER = "PDP:";

	public static final Pattern scopIDregex = Pattern.compile("d(....)(.)(.)");

	public static final String UNDERSCORE = "_";

	private static final String FILE_SEPARATOR = System.getProperty("file.separator");

	protected FileParsingParameters params;
	protected PDPProvider pdpprovider;

	private FetchBehavior fetchBehavior;
	private ObsoleteBehavior obsoleteBehavior;

	private String cachePath;

	// make sure IDs are loaded uniquely
	private Collection<String> currentlyLoading = Collections.synchronizedCollection(new TreeSet<String>());

	private String path;

	private boolean strictSCOP;
	private boolean useMmCif;

	/**
	 * Default AtomCache constructor.
	 * 
	 * Usually stores files in a temp directory, but this can be overriden by setting the PDB_DIR variable at runtime.
	 * 
	 * @see UserConfiguration#UserConfiguration()
	 */
	public AtomCache() {
		this(new UserConfiguration());
	}

	/**
	 * Creates an instance of an AtomCache that is pointed to the a particular path in the file system. It will use the same value for pdbFilePath and cachePath.
	 * 
	 * @param pdbFilePath
	 *            a directory in the file system to use as a location to cache files.
	 */
	public AtomCache(String pdbFilePath) {
		this(pdbFilePath,pdbFilePath);
	}

	/**
	 * Creates an instance of an AtomCache that is pointed to the a particular path in the file system.
	 * 
	 * @param pdbFilePath
	 *            a directory in the file system to use as a location to cache files.
	 * @param cachePath
	 */
	public AtomCache(String pdbFilePath, String cachePath) {
		
		logger.debug("Initialising AtomCache with pdbFilePath={}, cachePath={}",pdbFilePath, cachePath);

		if (!pdbFilePath.endsWith(FILE_SEPARATOR)) {
			pdbFilePath += FILE_SEPARATOR;
		}

		// we are caching the binary files that contain the PDBs gzipped
		// that is the most memory efficient way of caching...
		// set the input stream provider to caching mode
		System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");

		setPath(pdbFilePath);

		this.cachePath = cachePath;

		fetchBehavior = FetchBehavior.DEFAULT;
		obsoleteBehavior = ObsoleteBehavior.DEFAULT;

		currentlyLoading.clear();
		params = new FileParsingParameters();

		// we don't need this here
		params.setAlignSeqRes(false);
		// no secstruc either
		params.setParseSecStruc(false);
		//

		strictSCOP = true;

		setUseMmCif(true);

	}
	
	/**
	 * @param isSplit Ignored
	 * @deprecated isSplit parameter is ignored (4.0.0)
	 */
	@Deprecated
	public AtomCache(String pdbFilePath,boolean isSplit) {
		this(pdbFilePath);
	}
	/**
	 * @param isSplit Ignored
	 * @deprecated isSplit parameter is ignored (4.0.0)
	 */
	@Deprecated
	public AtomCache(String pdbFilePath, String cachePath,boolean isSplit) {
		this(pdbFilePath,cachePath);
	}

	/**
	 * Creates a new AtomCache object based on the provided UserConfiguration.
	 * 
	 * @param config
	 *            the UserConfiguration to use for this cache.
	 */
	public AtomCache(UserConfiguration config) {
		this(config.getPdbFilePath(), config.getCacheFilePath());
		fetchBehavior = config.getFetchBehavior();
		obsoleteBehavior = config.getObsoleteBehavior();
	}

	/**
	 * Returns the CA atoms for the provided name. See {@link #getStructure(String)} for supported naming conventions.
	 * <p>
	 * This method only works with protein chains. Use {@link #getRepresentativeAtoms(String)}
	 * for a more general solution.
	 * @param name
	 * @return an array of Atoms.
	 * @throws IOException
	 * @throws StructureException
	 * @see 
	 */
	public Atom[] getAtoms(String name) throws IOException, StructureException {

		Atom[] atoms = null;

		// System.out.println("loading " + name);
		Structure s = getStructure(name);

		atoms = StructureTools.getAtomCAArray(s);

		/*
		 * synchronized (cache){ cache.put(name, atoms); }
		 */

		return atoms;
	}
	/**
	 * Returns the representative atoms for the provided name.
	 * See {@link #getStructure(String)} for supported naming conventions.
	 * 
	 * @param name
	 * @return an array of Atoms.
	 * @throws IOException
	 * @throws StructureException
	 * @see 
	 */
	public Atom[] getRepresentativeAtoms(String name) throws IOException, StructureException {

		Atom[] atoms = null;

		// System.out.println("loading " + name);
		Structure s = getStructure(name);

		atoms = StructureTools.getRepresentativeAtomArray(s);

		/*
		 * synchronized (cache){ cache.put(name, atoms); }
		 */

		return atoms;
	}
	/**
	 * Loads the biological assembly for a given PDB ID and bioAssemblyId. If a bioAssemblyId > 0 is specified, the
	 * corresponding biological assembly file will be loaded. Note, the number of available biological unit files
	 * varies. Many entries don't have a biological assembly specified (i.e. NMR structures), many entries have only one
	 * biological assembly (bioAssemblyId=1), and a few structures have multiple biological assemblies. Set
	 * bioAssemblyFallback to true, to download the original PDB file in cases that a biological assembly file is not
	 * available.
	 * 
	 * @param pdbId
	 *            the PDB ID
	 * @param bioAssemblyId
	 *            the ID of the biological assembly
	 * @param bioAssemblyFallback
	 *            if true, try reading original PDB file in case the biological assembly file is not available
	 * @return a structure object
	 * @throws IOException
	 * @throws StructureException
	 * @author Peter Rose
	 * @since 3.2
	 */
	public Structure getBiologicalAssembly(String pdbId, int bioAssemblyId, boolean bioAssemblyFallback)
			throws StructureException, IOException {

		if (bioAssemblyId < 1) {
			throw new StructureException("bioAssemblyID must be greater than zero: " + pdbId + " bioAssemblyId "
					+ bioAssemblyId);
		}	
		Structure s = StructureIO.getBiologicalAssembly(pdbId, bioAssemblyId);
		
		if ( s == null && bioAssemblyFallback)
			return StructureIO.getBiologicalAssembly(pdbId, 0);
		
		return s;
	}

	/**
	 * Loads the default biological unit (*.pdb1.gz) file. If it is not available, the original PDB file will be loaded,
	 * i.e., for NMR structures, where the original files is also the biological assembly.
	 * 
	 * @param pdbId
	 *            the PDB ID
	 * @return a structure object
	 * @throws IOException
	 * @throws StructureException
	 * @since 3.2
	 */
	public Structure getBiologicalUnit(String pdbId) throws StructureException, IOException {
		int bioAssemblyId = 1;
		boolean bioAssemblyFallback = true;
		return getBiologicalAssembly(pdbId, bioAssemblyId, bioAssemblyFallback);
	}

	/**
	 * Returns the path that contains the caching file for utility data, such as domain definitions.
	 * 
	 * @return
	 */
	public String getCachePath() {
		return cachePath;
	}

	public FileParsingParameters getFileParsingParams() {
		return params;
	}

	/**
	 * Get the path that is used to cache PDB files.
	 * 
	 * @return path to a directory
	 */
	public String getPath() {
		return path;
	}

	public PDPProvider getPdpprovider() {
		return pdpprovider;
	}

	/**
	 * Request a Structure based on a <i>name</i>.
	 * 
	 * <pre>
	 * 		Formal specification for how to specify the <i>name</i>:
	 * 
	 * 		name     := pdbID
	 * 		               | pdbID '.' chainID
	 * 		               | pdbID '.' range
	 * 		               | scopID
	 * 		range         := '('? range (',' range)? ')'?
	 * 		               | chainID
	 * 		               | chainID '_' resNum '-' resNum
	 * 		pdbID         := [0-9][a-zA-Z0-9]{3}
	 * 		chainID       := [a-zA-Z0-9]
	 * 		scopID        := 'd' pdbID [a-z_][0-9_]
	 * 		resNum        := [-+]?[0-9]+[A-Za-z]?
	 * 
	 * 
	 * 		Example structures:
	 * 		1TIM     #whole structure
	 * 		4HHB.C     #single chain
	 * 		4GCR.A_1-83     #one domain, by residue number
	 * 		3AA0.A,B     #two chains treated as one structure
	 * 		d2bq6a1     #scop domain
	 * </pre>
	 * 
	 * With the additional set of rules:
	 * 
	 * <ul>
	 * <li>If only a PDB code is provided, the whole structure will be return including ligands, but the first model
	 * only (for NMR).
	 * <li>Chain IDs are case sensitive, PDB ids are not. To specify a particular chain write as: 4hhb.A or 4HHB.A</li>
	 * <li>To specify a SCOP domain write a scopId e.g. d2bq6a1. Some flexibility can be allowed in SCOP domain names,
	 * see {@link #setStrictSCOP(boolean)}</li>
	 * <li>URLs are accepted as well</li>
	 * </ul>
	 * 
	 * @param name
	 * @return a Structure object, or null if name appears improperly formated (eg too short, etc)
	 * @throws IOException
	 *             The PDB file cannot be cached due to IO errors
	 * @throws StructureException
	 *             The name appeared valid but did not correspond to a structure. Also thrown by some submethods upon
	 *             errors, eg for poorly formatted subranges.
	 */
	public Structure getStructure(String name) throws IOException, StructureException {

		if (name.length() < 4) {
			throw new IllegalArgumentException("Can't interpret IDs that are shorter than 4 characters!");
		}

		Structure n = null;

		boolean useChainNr = false;
		boolean useDomainInfo = false;
		String range = null;
		int chainNr = -1;


		StructureName structureName = new StructureName(name);

		String pdbId = null;
		String chainId = null;

		if (name.length() == 4) {

			pdbId = name;
			Structure s;
			if (useMmCif) {
				s = loadStructureFromCifByPdbId(pdbId);
			} else {
				s = loadStructureFromPdbByPdbId(pdbId);
			}
			return s;
		} else if (structureName.isScopName()) {

			// return based on SCOP domain ID
			return getStructureFromSCOPDomain(name);
		} else if (structureName.isCathID()) {
			return getStructureForCathDomain(structureName, CathFactory.getCathDatabase());
		} else if (name.length() == 6) {
			// name is PDB.CHAINID style (e.g. 4hhb.A)

			pdbId = name.substring(0, 4);
			if (name.substring(4, 5).equals(CHAIN_SPLIT_SYMBOL)) {
				chainId = name.substring(5, 6);
			} else if (name.substring(4, 5).equals(CHAIN_NR_SYMBOL)) {

				useChainNr = true;
				chainNr = Integer.parseInt(name.substring(5, 6));
			}

		} else if (name.startsWith("file:/") || name.startsWith("http:/")) {
			// this is a URL
			
			URL url = new URL(name);
			return getStructureFromURL(url);
			

		} else if (structureName.isPDPDomain()) {

			// this is a PDP domain definition

			return getPDPStructure(name);
			
		} else if (name.startsWith(BIOL_ASSEMBLY_IDENTIFIER)) {

			return getBioAssembly(name);

		} else if (name.length() > 6 && !name.startsWith(PDP_DOMAIN_IDENTIFIER)
				&& (name.contains(CHAIN_NR_SYMBOL) || name.contains(UNDERSCORE))
				&& !(name.startsWith("file:/") || name.startsWith("http:/"))

				) {

			// this is a name + range

			pdbId = name.substring(0, 4);
			// this ID has domain split information...
			useDomainInfo = true;
			range = name.substring(5);

		}

		// System.out.println("got: >" + name + "< " + pdbId + " " + chainId + " useChainNr:" + useChainNr + " "
		// +chainNr + " useDomainInfo:" + useDomainInfo + " " + range);

		if (pdbId == null) {

			return null;
		}

		while (checkLoading(pdbId)) {
			// waiting for loading to be finished...

			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				logger.error(e.getMessage());
			}

		}

		// long start = System.currentTimeMillis();

		Structure s;
		if (useMmCif) {
			s = loadStructureFromCifByPdbId(pdbId);
		} else {
			s = loadStructureFromPdbByPdbId(pdbId);
		}

		// long end = System.currentTimeMillis();
		// System.out.println("time to load " + pdbId + " " + (end-start) + "\t  size :" +
		// StructureTools.getNrAtoms(s) + "\t cached: " + cache.size());

		if (chainId == null && chainNr < 0 && range == null) {
			// we only want the 1st model in this case
			n = StructureTools.getReducedStructure(s, -1);
		} else {

			if (useChainNr) {
				// System.out.println("using ChainNr");
				n = StructureTools.getReducedStructure(s, chainNr);
			} else if (useDomainInfo) {
				// System.out.println("calling getSubRanges");
				n = StructureTools.getSubRanges(s, range);
			} else {
				// System.out.println("reducing Chain Id " + chainId);
				n = StructureTools.getReducedStructure(s, chainId);
			}
		}



		n.setName(name);
		return n;

	}

	/**
	 * Returns the representation of a {@link ScopDomain} as a BioJava {@link Structure} object.
	 * 
	 * @param domain
	 *            a SCOP domain
	 * @return a Structure object
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForDomain(ScopDomain domain) throws IOException, StructureException {
		return getStructureForDomain(domain, ScopFactory.getSCOP());
	}

	/**
	 * Returns the representation of a {@link ScopDomain} as a BioJava {@link Structure} object.
	 * 
	 * @param domain
	 *            a SCOP domain
	 * @param scopDatabase
	 *            A {@link ScopDatabase} to use
	 * @return a Structure object
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForDomain(ScopDomain domain, ScopDatabase scopDatabase) throws IOException,
			StructureException {
		return getStructureForDomain(domain, scopDatabase, false);
	}

	/**
	 * Returns the representation of a {@link ScopDomain} as a BioJava {@link Structure} object.
	 * 
	 * @param domain
	 *            a SCOP domain
	 * @param scopDatabase
	 *            A {@link ScopDatabase} to use
	 * @param strictLigandHandling
	 *            If set to false, hetero-atoms are included if and only if they belong to a chain to which the SCOP
	 *            domain belongs; if set to true, hetero-atoms are included if and only if they are strictly within the
	 *            definition (residue numbers) of the SCOP domain
	 * @return a Structure object
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForDomain(ScopDomain domain, ScopDatabase scopDatabase, boolean strictLigandHandling)
			throws IOException, StructureException {

		String pdbId = domain.getPdbId();
		Structure fullStructure = getStructure(pdbId);

		// build the substructure
		StringBuilder rangeString = new StringBuilder();
		Iterator<String> iter = domain.getRanges().iterator();
		while (iter.hasNext()) {
			rangeString.append(iter.next());
			if (iter.hasNext()) {
				rangeString.append(",");
			}
		}
		Structure structure = StructureTools.getSubRanges(fullStructure, rangeString.toString());
		structure.setName(domain.getScopId());
		structure.setPDBCode(domain.getScopId());

		// because ligands sometimes occur after TER records in PDB files, we may need to add some ligands back in
		// specifically, we add a ligand if and only if it occurs within the domain
		AtomPositionMap map = null;
		List<ResidueRangeAndLength> rrs = null;
		if (strictLigandHandling) {
			map = new AtomPositionMap(StructureTools.getAllAtomArray(fullStructure), AtomPositionMap.ANYTHING_MATCHER);
			rrs = ResidueRangeAndLength.parseMultiple(domain.getRanges(), map);
		}
		for (Chain chain : fullStructure.getChains()) {
			if (!structure.hasChain(chain.getChainID())) {
				continue; // we can't do anything with a chain our domain
			}
			// doesn't contain
			Chain newChain = structure.getChainByPDB(chain.getChainID());
			List<Group> ligands = StructureTools.filterLigands(chain.getAtomGroups());
			for (Group group : ligands) {
				boolean shouldContain = true;
				if (strictLigandHandling) {
					shouldContain = false; // whether the ligand occurs within the domain
					for (ResidueRange rr : rrs) {
						if (rr.contains(group.getResidueNumber(), map)) {
							shouldContain = true;
						}
					}
				}
				boolean alreadyContains = newChain.getAtomGroups().contains(group); // we don't want to add duplicate
																					// ligands
				if (shouldContain && !alreadyContains) {
					newChain.addGroup(group);
				}
			}
		}

		// build a more meaningful description for the new structure
		StringBuilder header = new StringBuilder();
		header.append(domain.getClassificationId());
		if (scopDatabase != null) {
			int sf = domain.getSuperfamilyId();
			ScopDescription description = scopDatabase.getScopDescriptionBySunid(sf);
			if (description != null) {
				header.append(" | ");
				header.append(description.getDescription());
			}
		}
		structure.getPDBHeader().setDescription(header.toString());

		return structure;

	}

	/**
	 * Returns the representation of a {@link ScopDomain} as a BioJava {@link Structure} object.
	 * 
	 * @param scopId
	 *            a SCOP Id
	 * @return a Structure object
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForDomain(String scopId) throws IOException, StructureException {
		return getStructureForDomain(scopId, ScopFactory.getSCOP());
	}

	/**
	 * Returns the representation of a {@link ScopDomain} as a BioJava {@link Structure} object.
	 * 
	 * @param scopId
	 *            a SCOP Id
	 * @param scopDatabase
	 *            A {@link ScopDatabase} to use
	 * @return a Structure object
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForDomain(String scopId, ScopDatabase scopDatabase) throws IOException,
			StructureException {
		ScopDomain domain = scopDatabase.getDomainByScopID(scopId);
		return getStructureForDomain(domain, scopDatabase);
	}

	/**
	 * Does the cache automatically download files that are missing from the local installation from the PDB FTP site?
	 *
	 * @return flag
	 * @deprecated Use {@link #getFetchBehavior()}
	 */
	@Deprecated
	public boolean isAutoFetch() {
		return fetchBehavior != FetchBehavior.LOCAL_ONLY;
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 *
	 * @return the fetchCurrent
	 * @deprecated Use {@link FileParsingParameters#getObsoleteBehavior()} instead (4.0.0)
	 */
	@Deprecated
	public boolean isFetchCurrent() {
		return getObsoleteBehavior() == ObsoleteBehavior.FETCH_CURRENT;
	}

	/**
	 * forces the cache to fetch the file if its status is OBSOLETE. This feature has a higher priority than
	 * {@link #setFetchCurrent(boolean)}.<br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 *
	 * @return the fetchFileEvenIfObsolete
	 * @author Amr AL-Hossary
	 * @see #fetchCurrent
	 * @since 3.0.2
	 * @deprecated Use {@link FileParsingParameters#getObsoleteBehavior()} instead (4.0.0)
	 */
	@Deprecated
	public boolean isFetchFileEvenIfObsolete() {
		return getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE;
	}


	/**
	 * Reports whether strict scop naming will be enforced, or whether this AtomCache should try to guess some simple
	 * variants on scop domains.
	 * 
	 * @return true if scop names should be used strictly with no guessing
	 */
	public boolean isStrictSCOP() {
		return strictSCOP;
	}

	/**
	 * Send a signal to the cache that the system is shutting down. Notifies underlying SerializableCache instances to
	 * flush themselves...
	 */
	public void notifyShutdown() {
		// System.out.println(" AtomCache got notify shutdown..");
		if (pdpprovider != null) {
			if (pdpprovider instanceof RemotePDPProvider) {
				RemotePDPProvider remotePDP = (RemotePDPProvider) pdpprovider;
				remotePDP.flushCache();
			}
		}

		// todo: use a SCOP implementation that is backed by SerializableCache
		ScopDatabase scopInstallation = ScopFactory.getSCOP();
		if (scopInstallation != null) {
			if (scopInstallation instanceof CachedRemoteScopInstallation) {
				CachedRemoteScopInstallation cacheScop = (CachedRemoteScopInstallation) scopInstallation;
				cacheScop.flushCache();
			}
		}

	}

	/**
	 * Does the cache automatically download files that are missing from the local installation from the PDB FTP site?
	 * 
	 * @param autoFetch
	 *            flag
	 * @deprecated Use {@link #getFetchBehavior()}
	 */
	@Deprecated
	public void setAutoFetch(boolean autoFetch) {
		if(autoFetch) {
			setFetchBehavior(FetchBehavior.DEFAULT);
		} else {
			setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		}
	}

	/**
	 * set the location at which utility data should be cached.
	 * 
	 * @param cachePath
	 */
	public void setCachePath(String cachePath) {
		this.cachePath = cachePath;
	}

	/**
	 * if enabled, the reader searches for the newest possible PDB ID, if not present in he local installation. The
	 * {@link #setFetchFileEvenIfObsolete(boolean)} function has a higher priority than this function.<br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 *
	 * @param fetchCurrent
	 *            the fetchCurrent to set
	 * @author Amr AL-Hossary
	 * @see #setFetchFileEvenIfObsolete(boolean)
	 * @since 3.0.2
	 * @deprecated Use {@link FileParsingParameters#setObsoleteBehavior()} instead (4.0.0)
	 */
	@Deprecated
	public void setFetchCurrent(boolean fetchNewestCurrent) {
		if(fetchNewestCurrent) {
			setObsoleteBehavior(ObsoleteBehavior.FETCH_CURRENT);
		} else {
			if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_CURRENT) {
				setObsoleteBehavior(ObsoleteBehavior.DEFAULT);
			}
		}
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * 
	 * @param fetchFileEvenIfObsolete
	 *            the fetchFileEvenIfObsolete to set
	 * @deprecated Use {@link FileParsingParameters#setObsoleteBehavior()} instead (4.0.0)
	 */
	@Deprecated
	public void setFetchFileEvenIfObsolete(boolean fetchFileEvenIfObsolete) {
		if(fetchFileEvenIfObsolete) {
			setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		} else {
			if(getObsoleteBehavior() == ObsoleteBehavior.FETCH_OBSOLETE) {
				setObsoleteBehavior(ObsoleteBehavior.DEFAULT);
			}
		}
	}

	public void setFileParsingParams(FileParsingParameters params) {
		this.params = params;
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
	 * Set the behavior for fetching files from the server
	 * @param fetchBehavior
	 */
	public void setFetchBehavior(FetchBehavior fetchBehavior) {
		this.fetchBehavior = fetchBehavior;
	}

	/**
	 * Set the path that is used to cache PDB files.
	 * 
	 * @param path
	 *            to a directory
	 */
	public void setPath(String path) {
		this.path = FileDownloadUtils.expandUserHome(path);
	}

	public void setPdpprovider(PDPProvider pdpprovider) {
		this.pdpprovider = pdpprovider;
	}


	/**
	 * When strictSCOP is enabled, SCOP domain identifiers (eg 'd1gbga_') are matched literally to the SCOP database.
	 * 
	 * When disabled, some simple mistakes are corrected automatically. For instance, the invalid identifier 'd1gbg__'
	 * would be corrected to 'd1gbga_' automatically.
	 * 
	 * @param strictSCOP
	 *            Indicates whether strict scop names should be used.
	 */
	public void setStrictSCOP(boolean strictSCOP) {
		this.strictSCOP = strictSCOP;
	}

	/**
	 * @return the useMmCif
	 */
	public boolean isUseMmCif() {
		return useMmCif;
	}

	/**
	 * @param useMmCif
	 *            the useMmCif to set
	 */
	public void setUseMmCif(boolean useMmCif) {
		this.useMmCif = useMmCif;
		
		if ( useMmCif) {
			// get bio assembly from mmcif file

			BioUnitDataProviderFactory.setBioUnitDataProvider(BioUnitDataProviderFactory.mmcifProviderClassName);

		} else {
		
			BioUnitDataProviderFactory.setBioUnitDataProvider(BioUnitDataProviderFactory.pdbProviderClassName);
			
		}
	}

	private boolean checkLoading(String name) {
		return currentlyLoading.contains(name);

	}

	private Structure getBioAssembly(String name) throws IOException, StructureException {

	
		// can be specified as:
		// BIO:1fah - first one
		// BIO:1fah:0 - asym unit
		// BIO:1fah:1 - first one
		// BIO:1fah:2 - second one

		String pdbId = name.substring(4, 8);
		int biolNr = 1;
		if (name.length() > 8) {
			biolNr = Integer.parseInt(name.substring(9, name.length()));
		}

		Structure s= StructureIO.getBiologicalAssembly(pdbId, biolNr);

		return s;
	}

	private Structure getPDPStructure(String pdpDomainName) {

		// System.out.println("loading PDP domain from server " + pdpDomainName);
		if (pdpprovider == null) {
			pdpprovider = new RemotePDPProvider(true);
		}

		return pdpprovider.getDomain(pdpDomainName, this);

	}

	private ScopDomain getScopDomain(String scopId) {
		return ScopFactory.getSCOP().getDomainByScopID(scopId);
	}

	/**
	 * Returns a {@link Structure} corresponding to the CATH identifier supplied in {@code structureName}, using the the {@link CathDatabase}
	 * at {@link CathFactory#getCathDatabase()}.
	 */
	public Structure getStructureForCathDomain(StructureName structureName) throws IOException, StructureException {
		return getStructureForCathDomain(structureName, CathFactory.getCathDatabase());
	}

	/**
	 * Returns a {@link Structure} corresponding to the CATH identifier supplied in {@code structureName}, using the specified {@link CathDatabase}.
	 */
	public Structure getStructureForCathDomain(StructureName structureName, CathDatabase cathInstall) throws IOException, StructureException {

		CathDomain cathDomain = cathInstall.getDomainByCathId(structureName.getName());

		List<CathSegment> segments = cathDomain.getSegments();

		StringWriter range = new StringWriter();

		int rangePos = 0;
		String chainId = structureName.getChainId();
		for (CathSegment segment : segments) {
			rangePos++;

			range.append(chainId);
			range.append("_");

			range.append(segment.getStart());
			range.append("-");
			range.append(segment.getStop());
			if (segments.size() > 1 && rangePos < segments.size()) {
				range.append(",");
			}
		}

		String pdbId = structureName.getPdbId();

		Structure s = getStructure(pdbId);

		String rangeS = range.toString();
		
		Structure n = StructureTools.getSubRanges(s, rangeS);
		
		// add the ligands of the chain...

		Chain newChain = n.getChainByPDB(structureName.getChainId());
		Chain origChain = s.getChainByPDB(structureName.getChainId());
		List<Group> ligands = origChain.getAtomLigands();

		for (Group g : ligands) {
			if (!newChain.getAtomGroups().contains(g)) {
				newChain.addGroup(g);
			}
		}

		// set new Header..
		n.setName(structureName.getName());
		n.setPDBCode(structureName.getPdbId());

		n.getPDBHeader().setDescription(cathDomain.getDomainName());

		return n;
	}

	private Structure getStructureFromSCOPDomain(String name) throws IOException, StructureException {
		// looks like a SCOP domain!
		ScopDomain domain;
		if (strictSCOP) {
			domain = getScopDomain(name);
		} else {
			domain = guessScopDomain(name);
		}

		//System.out.println(domain);
		if (domain != null) {
			Structure s = getStructureForDomain(domain);
			return s;
		}

		// Guessing didn't work, so just use the PDBID and Chain from name
		if (!strictSCOP) {
			Matcher scopMatch = scopIDregex.matcher(name);
			if (scopMatch.matches()) {
				String pdbID = scopMatch.group(1);
				String chainID = scopMatch.group(2);

				// None of the actual SCOP domains match. Guess that '_' means 'whole chain'
				if (!chainID.equals("_")) {
					// Add chain identifier
					pdbID += "." + scopMatch.group(2);
				}
				// Fetch the structure by pdb id
				Structure struct = getStructure(pdbID);
				if (struct != null) {
					System.err.println("Trying chain " + pdbID);
				}

				return struct;
			}
		}

		throw new StructureException("Unable to get structure for SCOP domain: " + name);
	}

	private Structure getStructureFromURL(URL url) throws IOException, StructureException {
		// looks like a URL for a file was provided:
		System.out.println("fetching structure from URL:" + url);

		String queryS = url.getQuery();
		String chainId = null;
		
		String fullu = url.toString();
		if (fullu.startsWith("file:") && fullu.endsWith("?" + queryS)) {
			String newu = fullu.substring(0, fullu.length() - ("?" + queryS).length());
			url = new URL(newu);
		}
		
		if (queryS != null && queryS.startsWith("chainId=")) {
			chainId = queryS.substring(8);
		}
		
		if ( fullu.contains(".cif")) {

			// need to do mmcif parsing!

			InputStreamProvider prov = new InputStreamProvider();
			InputStream inStream =  prov.getInputStream(url);

			MMcifParser parser = new SimpleMMcifParser();

			SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
			consumer.setFileParsingParameters(params);


			parser.addMMcifConsumer(consumer);

			try {
				parser.parse(new BufferedReader(new InputStreamReader(inStream)));
			} catch (IOException e){
				e.printStackTrace();
			}

			// now get the protein structure.
			Structure cifStructure = consumer.getStructure();
			if (chainId == null) {
				return StructureTools.getReducedStructure(cifStructure, -1);
			} else {
				return StructureTools.getReducedStructure(cifStructure, chainId);
			}

		} else {

			// pdb file based parsing

			PDBFileReader reader = new PDBFileReader(path);
			reader.setFetchBehavior(fetchBehavior);
			reader.setObsoleteBehavior(obsoleteBehavior);

			reader.setFileParsingParameters(params);
			Structure s = reader.getStructure(url);
			if (chainId == null) {
				return StructureTools.getReducedStructure(s, -1);
			} else {
				return StructureTools.getReducedStructure(s, chainId);
			}
		}
	}

	/**
	 * <p>
	 * Guess a scop domain. If an exact match is found, return that.
	 * 
	 * <p>
	 * Otherwise, return the first scop domain found for the specified protein such that
	 * <ul>
	 * <li>The chains match, or one of the chains is '_' or '.'.
	 * <li>The domains match, or one of the domains is '_'.
	 * </ul>
	 * 
	 * 
	 * @param name
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	private ScopDomain guessScopDomain(String name) throws IOException, StructureException {
		List<ScopDomain> matches = new LinkedList<ScopDomain>();

		// Try exact match first
		ScopDomain domain = getScopDomain(name);
		if (domain != null) {
			return domain;
		}

		// Didn't work. Guess it!
		logger.warn("Warning, could not find SCOP domain: " + name);

		Matcher scopMatch = scopIDregex.matcher(name);
		if (scopMatch.matches()) {
			String pdbID = scopMatch.group(1);
			String chainID = scopMatch.group(2);
			String domainID = scopMatch.group(3);

			for (ScopDomain potentialSCOP : ScopFactory.getSCOP().getDomainsForPDB(pdbID)) {
				Matcher potMatch = scopIDregex.matcher(potentialSCOP.getScopId());
				if (potMatch.matches()) {
					if (chainID.equals(potMatch.group(2)) || chainID.equals("_") || chainID.equals(".")
							|| potMatch.group(2).equals("_") || potMatch.group(2).equals(".")) {
						if (domainID.equals(potMatch.group(3)) || domainID.equals("_") || potMatch.group(3).equals("_")) {
							// Match, or near match
							matches.add(potentialSCOP);
						}
					}
				}
			}
		}

		Iterator<ScopDomain> match = matches.iterator();
		if (match.hasNext()) {
			ScopDomain bestMatch = match.next();
			StringBuilder warnMsg = new StringBuilder();
			warnMsg.append("Trying domain " + bestMatch.getScopId() + ".");
			if (match.hasNext()) {
				warnMsg.append(" Other possibilities: ");
				while (match.hasNext()) {
					warnMsg.append(match.next().getScopId() + " ");
				}
			}
			warnMsg.append(System.getProperty("line.separator"));
			logger.warn(warnMsg.toString());
			return bestMatch;
		} else {
			return null;
		}
	}

	protected void flagLoading(String name) {
		if (!currentlyLoading.contains(name)) {
			
			currentlyLoading.add(name);
		}
	}

	protected void flagLoadingFinished(String name) {
	
		currentlyLoading.remove(name);
	}

	protected Structure loadStructureFromCifByPdbId(String pdbId) throws IOException, StructureException {

		Structure s;
		flagLoading(pdbId);
		try {
			MMCIFFileReader reader = new MMCIFFileReader(path);
			reader.setFetchBehavior(fetchBehavior);
			reader.setObsoleteBehavior(obsoleteBehavior);

			reader.setFileParsingParameters(params);

			s = reader.getStructureById(pdbId.toLowerCase());

		} finally {
			flagLoadingFinished(pdbId);			
		}

		return s;
	}

	protected Structure loadStructureFromPdbByPdbId(String pdbId) throws IOException, StructureException {

		Structure s;
		flagLoading(pdbId);
		try {
			PDBFileReader reader = new PDBFileReader(path);
			reader.setFetchBehavior(fetchBehavior);
			reader.setObsoleteBehavior(obsoleteBehavior);

			reader.setFileParsingParameters(params);

			s = reader.getStructureById(pdbId.toLowerCase());

		} finally {
			flagLoadingFinished(pdbId);
		}
		
		return s;
	}

}