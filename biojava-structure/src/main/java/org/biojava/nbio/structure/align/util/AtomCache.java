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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;
import org.biojava.nbio.structure.domain.PDPProvider;
import org.biojava.nbio.structure.domain.RemotePDPProvider;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.MMTFFileReader;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.biojava.nbio.structure.scop.CachedRemoteScopInstallation;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopDescription;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
	
	/**
	 * The default output bioassembly style: if true the bioassemblies are multimodel,
	 * if false the bioassemblies are flat with renamed chains for symmetry-partners.
	 */
	public static final boolean DEFAULT_BIOASSEMBLY_STYLE = false;

	public static final String BIOL_ASSEMBLY_IDENTIFIER = "BIO:";
	public static final String CHAIN_NR_SYMBOL = ":";
	public static final String CHAIN_SPLIT_SYMBOL = ".";

	public static final String PDP_DOMAIN_IDENTIFIER = "PDP:";

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

	private boolean useMmCif;
	private boolean useMmtf;

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

		setUseMmCif(false);
		setUseMmtf(true);

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
		useMmCif = config.getFileFormat().equals( UserConfiguration.MMCIF_FORMAT );

		if ( useMmCif)
			useMmtf = false;

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
		return getAtoms(new StructureName(name));
	}
	public Atom[] getAtoms(StructureIdentifier name) throws IOException, StructureException {

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
		return getRepresentativeAtoms(new StructureName(name));
	}
	
	public Atom[] getRepresentativeAtoms(StructureIdentifier name) throws IOException, StructureException {

		Atom[] atoms = null;

		Structure s = getStructure(name);

		atoms = StructureTools.getRepresentativeAtomArray(s);

		/*
		 * synchronized (cache){ cache.put(name, atoms); }
		 */

		return atoms;
	}
	
	/**
	 * Returns the biological assembly for a given PDB ID and bioAssemblyId, by building the 
	 * assembly from the biounit annotations found in {@link Structure#getPDBHeader()}
	 * <p>
	 * Note, the number of available biological unit files
	 * varies. Many entries don't have a biological assembly specified (e.g. NMR structures), many entries have only one
	 * biological assembly (bioAssemblyId=1), and some structures have multiple biological assemblies.
	 *
	 * @param pdbId
	 *            the PDB ID
	 * @param bioAssemblyId
	 *            the 1-based index of the biological assembly (0 gets the asymmetric unit)
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).             
	 * @return a structure object
	 * @throws IOException
	 * @throws StructureException if biassemblyId < 0 or other problems while loading structure
	 * @author Peter Rose
	 * @since 3.2
	 */
	public Structure getBiologicalAssembly(String pdbId, int bioAssemblyId, boolean multiModel)
			throws StructureException, IOException {

		if (bioAssemblyId < 0) {
			throw new StructureException("bioAssemblyID must be nonnegative: " + pdbId + " bioAssemblyId "
					+ bioAssemblyId);
		}
		
		boolean prevIsParseBioAssembly = getFileParsingParams().isParseBioAssembly();
		
		if (!getFileParsingParams().isParseBioAssembly()) {
			getFileParsingParams().setParseBioAssembly(true);
		}
		
		Structure asymUnit = getStructureForPdbId(pdbId);
		
		getFileParsingParams().setParseBioAssembly(prevIsParseBioAssembly);
		
		if (asymUnit.getPDBHeader() == null || asymUnit.getPDBHeader().getBioAssemblies()==null) {
			logger.info("No bioassembly information found for {}, returning asymmetric unit as biological assembly", pdbId);
			return asymUnit; 
		}

		// 0 ... asym unit
		if ( bioAssemblyId == 0) {
			logger.info("Requested biological assembly 0 for PDB id "+pdbId+", returning asymmetric unit");
			return asymUnit;
		}
		// does it exist?
		if (!asymUnit.getPDBHeader().getBioAssemblies().containsKey(bioAssemblyId)) {
			throw new StructureException("No biological assembly available for biological assembly id " + bioAssemblyId + " of " + pdbId);
		}

		List<BiologicalAssemblyTransformation> transformations =
				asymUnit.getPDBHeader().getBioAssemblies().get(bioAssemblyId).getTransforms();


		if ( transformations == null || transformations.size() == 0){

			throw new StructureException("Could not load transformations to recreate biological assembly id " + bioAssemblyId + " of " + pdbId);
			
		}
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		// if we use mmcif or mmtf, then we need to pass useAsymIds=true
		boolean useAsymIds = false;
		if (useMmCif) useAsymIds = true;
		if (useMmtf) useAsymIds = true;
		return builder.rebuildQuaternaryStructure(asymUnit, transformations, useAsymIds, multiModel);
		
	}

	/**
	 * Returns the default biological unit (bioassemblyId=1, known in PDB as pdb1.gz). If it is not available,
	 * the asymmetric unit will be returned, e.g. for NMR structures.
	 *
	 * <p>Biological assemblies can also be accessed using
	 * <tt>getStructure("BIO:<i>[pdbId]</i>")</tt>
	 * @param pdbId the PDB id
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).  
	 * @return a structure object
	 * @throws IOException
	 * @throws StructureException
	 * @since 4.2
	 */
	public Structure getBiologicalAssembly(String pdbId, boolean multiModel) throws StructureException, IOException {
		
		boolean prevIsParseBioAssembly = getFileParsingParams().isParseBioAssembly();
		
		if (!getFileParsingParams().isParseBioAssembly()) {
			getFileParsingParams().setParseBioAssembly(true);
		}
		
		Structure asymUnit = getStructureForPdbId(pdbId);
		
		getFileParsingParams().setParseBioAssembly(prevIsParseBioAssembly);

		
		if (asymUnit.getPDBHeader() == null || asymUnit.getPDBHeader().getBioAssemblies()==null) {
			logger.info("No bioassembly information found for {}, returning asymmetric unit as biological assembly", pdbId);
			return asymUnit; 
		}

		int bioAssemblyId = 1;
		
		// does it exist?
		if (!asymUnit.getPDBHeader().getBioAssemblies().containsKey(bioAssemblyId)) {
			return asymUnit;
		}

		List<BiologicalAssemblyTransformation> transformations =
				asymUnit.getPDBHeader().getBioAssemblies().get(bioAssemblyId).getTransforms();


		if ( transformations == null || transformations.size() == 0){

			throw new StructureException("Could not load transformations to recreate biological assembly id " + bioAssemblyId + " of " + pdbId);
			
		}
		
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		// if we use mmcif or mmtf, then we need to pass useAsymIds=true
		boolean useAsymIds = false;
		if (useMmCif) useAsymIds = true;
		if (useMmtf) useAsymIds = true;
		return builder.rebuildQuaternaryStructure(asymUnit, transformations, useAsymIds, multiModel);
		
	}

	/**
	 * Returns all biological assemblies for given PDB id.
	 * @param pdbId
	 * @param multiModel if true the output Structure will be a multi-model one with one transformId per model, 
	 * if false the outputStructure will be as the original with added chains with renamed asymIds (in the form originalAsymId_transformId and originalAuthId_transformId).  
	 * @return
	 * @throws StructureException
	 * @throws IOException
	 * @since 5.0
	 */
	public List<Structure> getBiologicalAssemblies(String pdbId, boolean multiModel) throws StructureException, IOException {
		
		List<Structure> assemblies = new ArrayList<>();
		
		boolean prevIsParseBioAssembly = getFileParsingParams().isParseBioAssembly();
		
		if (!getFileParsingParams().isParseBioAssembly()) {
			getFileParsingParams().setParseBioAssembly(true);
		}
		
		Structure asymUnit = getStructureForPdbId(pdbId);
		
		getFileParsingParams().setParseBioAssembly(prevIsParseBioAssembly);
		

		if (asymUnit.getPDBHeader() == null || asymUnit.getPDBHeader().getBioAssemblies()==null) {
			logger.info("No bioassembly information found for {}, returning asymmetric unit as the only biological assembly", pdbId);
			assemblies.add(asymUnit);
			return assemblies; 
		}


		for (int bioAssemblyId : asymUnit.getPDBHeader().getBioAssemblies().keySet()) {	
			List<BiologicalAssemblyTransformation> transformations =
					asymUnit.getPDBHeader().getBioAssemblies().get(bioAssemblyId).getTransforms();


			if ( transformations == null || transformations.size() == 0){

				logger.info("Could not load transformations to recreate biological assembly id " + bioAssemblyId + " of " + pdbId+". Assembly id will be missing in biological assemblies.");
				continue;
			}

			BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

			// if we use mmcif or mmtf, then we need to pass useAsymIds=true
			boolean useAsymIds = false;
			if (useMmCif) useAsymIds = true;
			if (useMmtf) useAsymIds = true;
			Structure s = builder.rebuildQuaternaryStructure(asymUnit, transformations, useAsymIds, multiModel);
			assemblies.add(s);
		}
		return assemblies;
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
	 * <p>Note that this method should not be used in StructureIdentifier
	 * implementations to avoid circular calls.
	 * @param name
	 * @return a Structure object, or null if name appears improperly formated (eg too short, etc)
	 * @throws IOException
	 *             The PDB file cannot be cached due to IO errors
	 * @throws StructureException
	 *             The name appeared valid but did not correspond to a structure. Also thrown by some submethods upon
	 *             errors, eg for poorly formatted subranges.
	 */
	public Structure getStructure(String name) throws IOException, StructureException {
		StructureName structureName = new StructureName(name);

		return getStructure(structureName);
	}

	/**
	 * Get the structure corresponding to the given {@link StructureIdentifier}.
	 * Equivalent to calling {@link StructureIdentifier#loadStructure(AtomCache)}
	 * followed by {@link StructureIdentifier#reduce(Structure)}.
	 *
	 * <p>Note that this method should not be used in StructureIdentifier
	 * implementations to avoid circular calls.
	 * @param strucId
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructure(StructureIdentifier strucId) throws IOException, StructureException {
		Structure s = strucId.loadStructure(this);
		Structure r = strucId.reduce(s);
		r.setStructureIdentifier(strucId);
		return r;
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
		Structure fullStructure = getStructureForPdbId(pdbId);
		Structure structure = domain.reduce(fullStructure);

		// TODO It would be better to move all of this into the reduce method,
		// but that would require ligand handling properties in StructureIdentifiers

		// because ligands sometimes occur after TER records in PDB files, we may need to add some ligands back in
		// specifically, we add a ligand if and only if it occurs within the domain
		AtomPositionMap map = null;
		List<ResidueRangeAndLength> rrs = null;
		if (strictLigandHandling) {
			map = new AtomPositionMap(StructureTools.getAllAtomArray(fullStructure), AtomPositionMap.ANYTHING_MATCHER);
			rrs = ResidueRangeAndLength.parseMultiple(domain.getRanges(), map);
		}
		for (Chain chain : fullStructure.getNonPolyChains()) {

			if (!structure.hasPdbChain(chain.getName())) {
				continue; // we can't do anything with a chain our domain
			}

			Chain newChain;
			if (! structure.hasNonPolyChain(chain.getId())) {
				newChain = new ChainImpl();
				newChain.setId(chain.getId());
				newChain.setName(chain.getName());
				newChain.setEntityInfo(chain.getEntityInfo());
				structure.addChain(newChain);
			} else {
				newChain = structure.getNonPolyChain(chain.getId());
			}
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
	 * set the location at which utility data should be cached.
	 *
	 * @param cachePath
	 */
	public void setCachePath(String cachePath) {
		this.cachePath = cachePath;
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
		// Either way the user wants to use PDB or MMCIF
		this.useMmtf = false;
	}
	
	/**
	 * Set whether to use mmtf.
	 * @param bool the input boolean to set
	 */
	public void setUseMmtf(boolean useMmtf) {
		this.useMmtf = useMmtf;
		if(useMmtf){
			useMmCif=false;
		}
		
	}

	/** Returns useMmtf flag
	 *
	 * @return true if will load data via mmtf file format
     */
	public boolean isUseMmtf(){
		return this.useMmtf;
	}

	private boolean checkLoading(String name) {
		return currentlyLoading.contains(name);

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

		CathDomain cathDomain = cathInstall.getDomainByCathId(structureName.getIdentifier());

		Structure s = getStructureForPdbId(cathDomain.getIdentifier());
		Structure n = cathDomain.reduce(s);

		// add the ligands of the chain...

		Chain newChain = n.getPolyChainByPDB(structureName.getChainId());
		List<Chain> origChains = s.getNonPolyChainsByPDB(structureName.getChainId());
		for ( Chain origChain : origChains) {
			List<Group> ligands = origChain.getAtomGroups();

			for (Group g : ligands) {
				if (!newChain.getAtomGroups().contains(g)) {
					newChain.addGroup(g);
				}
			}
		}

		return n;
	}

	protected void flagLoading(String name) {
		if (!currentlyLoading.contains(name)) {

			currentlyLoading.add(name);
		}
	}

	protected void flagLoadingFinished(String name) {

		currentlyLoading.remove(name);
	}

	/**
	 * Loads a structure directly by PDB ID
	 * @param pdbId
	 * @return
	 * @throws IOException
	 * @throws StructureException
	 */
	public Structure getStructureForPdbId(String pdbId) throws IOException, StructureException {
		if(pdbId == null)
			return null;
		if(pdbId.length() != 4) {
			throw new StructureException("Unrecognized PDB ID: "+pdbId);
		}
		while (checkLoading(pdbId)) {
			// waiting for loading to be finished...

			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				logger.error(e.getMessage());
			}

		}

		Structure s;
		if (useMmtf) {
			logger.debug("loading from mmtf");
			s = loadStructureFromMmtfByPdbId(pdbId);
		}
		else if (useMmCif) {
			logger.debug("loading from mmcif");
			s = loadStructureFromCifByPdbId(pdbId);
		} else {
			logger.debug("loading from pdb");
			s = loadStructureFromPdbByPdbId(pdbId);
		}
		return s;
	}

	/**
	 * Load a {@link Structure} from MMTF either from the local file system.
	 * @param pdbId the input PDB id
	 * @return the {@link Structure} object of the parsed structure
	 * @throws IOException error reading from Web or file system
	 */
	private Structure loadStructureFromMmtfByPdbId(String pdbId) throws IOException {
		logger.debug("Loading structure {} from mmtf file.", pdbId);
		MMTFFileReader reader = new MMTFFileReader();
		reader.setFetchBehavior(fetchBehavior);
		reader.setObsoleteBehavior(obsoleteBehavior);
		Structure structure = reader.getStructureById(pdbId.toLowerCase());
		return structure;
	}

	protected Structure loadStructureFromCifByPdbId(String pdbId) throws IOException, StructureException {

		logger.debug("Loading structure {} from mmCIF file {}.", pdbId, path);
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

		logger.debug("Loading structure {} from PDB file {}.", pdbId, path);
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
