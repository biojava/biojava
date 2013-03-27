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
 * 
 * 
 *
 */
package org.biojava.bio.structure.align.util;

import java.io.IOException;
import java.io.StringWriter;
import java.net.URL;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.Iterator;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.cath.CathDomain;
import org.biojava.bio.structure.cath.CathInstallation;
import org.biojava.bio.structure.cath.CathSegment;
import org.biojava.bio.structure.domain.PDPProvider;
import org.biojava.bio.structure.domain.RemotePDPProvider;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.scop.CachedRemoteScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.core.util.InputStreamProvider;
import org.biojava3.structure.StructureIO;



/** A utility class that provides easy access to Structure objects. If you are running a
 *  script that is frequently re-using the same PDB structures, the AtomCache keeps an 
 *  in-memory cache of the files for quicker access. The cache is a soft-cache, this 
 *  means it won't cause out of memory exceptions, but garbage collects the data if the 
 *  Java virtual machine needs to free up space. The AtomCache is thread-safe.
 * 
 * @author Andreas Prlic
 * @author Spencer Bliven
 * @author Peter Rose
 * @since 3.0
 */
public class AtomCache {

	public static final String CHAIN_NR_SYMBOL = ":";
	public static final String UNDERSCORE = "_";
	public static final String CHAIN_SPLIT_SYMBOL = ".";

	private static final String FILE_SEPARATOR = System.getProperty("file.separator");

	String path;

	String cachePath;
	
	// make sure IDs are loaded uniquely
	Collection<String> currentlyLoading = Collections.synchronizedCollection(new TreeSet<String>());

	private  ScopDatabase scopInstallation ;

	protected PDPProvider pdpprovider;

	boolean autoFetch;
	boolean isSplit;
	boolean strictSCOP;
	
	protected FileParsingParameters params;

	private boolean fetchFileEvenIfObsolete;

	private boolean fetchCurrent;

	public static final String PDP_DOMAIN_IDENTIFIER = "PDP:";
	public static final String BIOL_ASSEMBLY_IDENTIFIER = "BIO:";


	/**
	 * Default AtomCache constructor.
	 * 
	 * Usually stores files in a temp directory, but this can be overriden
	 * by setting the PDB_DIR variable at runtime.
	 * 
	 * @see UserConfiguration#UserConfiguration()
	 */
	public AtomCache() {
		this(new UserConfiguration());
	}

	/** Creates an instance of an AtomCache that is pointed to the a particular
	 * path in the file system.
	 * 
	 * @param pdbFilePath a directory in the file system to use as a location to cache files.
	 * @param isSplit a flag to indicate if the directory organisation is "split" as on the PDB ftp servers, or if all files are contained in one directory.
	 */
	public AtomCache(String pdbFilePath, boolean isSplit){

		if ( ! pdbFilePath.endsWith(FILE_SEPARATOR))
			pdbFilePath += FILE_SEPARATOR;

		// we are caching the binary files that contain the PDBs gzipped
		// that is the most memory efficient way of caching...
		// set the input stream provider to caching mode
		System.setProperty(InputStreamProvider.CACHE_PROPERTY, "true");

		path = pdbFilePath;		
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);
		
		
		String tmpCache = System.getProperty(AbstractUserArgumentProcessor.CACHE_DIR);
		if ( tmpCache == null || tmpCache.equals(""))
			tmpCache = pdbFilePath;
		
		cachePath = tmpCache;
		System.setProperty(AbstractUserArgumentProcessor.CACHE_DIR, cachePath);
		
		//this.cache = cache;
		this.isSplit = isSplit;

		autoFetch = true;
		fetchFileEvenIfObsolete = false;
		fetchCurrent = false;

		currentlyLoading.clear();
		params = new FileParsingParameters();

		// we don't need this here
		params.setAlignSeqRes(false);
		// no secstruc either 
		params.setParseSecStruc(false);
		// 


		this.strictSCOP = true;

		scopInstallation = null;
	}


	/** Creates a new AtomCache object based on the provided UserConfiguration.
	 * 
	 * @param config the UserConfiguration to use for this cache.
	 */
	public AtomCache(UserConfiguration config){
		this(config.getPdbFilePath(),config.isSplit());
		autoFetch = config.getAutoFetch();
	}


	/** Get the path that is used to cache PDB files.
	 * 
	 * @return path to a directory
	 */
	public String getPath() {
		return path;
	}

	/** Set the path that is used to cache PDB files.
	 * 
	 * @param path to a directory
	 */
	public void setPath(String path) {
		System.setProperty(AbstractUserArgumentProcessor.PDB_DIR,path);
		this.path = path;
	}
	
	
	/** Returns the path that contains the caching file for utility data, such as domain definitons.
	 *  
	 * @return
	 */
	public String getCachePath(){
		if ( cachePath == null || cachePath.equals(""))
			return getPath();
		return cachePath;
	}

	/** set the location at which utility data should be cached.
	 * 
	 * @param cachePath
	 */
	public void setCachePath(String cachePath){
		this.cachePath = cachePath;
		System.setProperty(AbstractUserArgumentProcessor.CACHE_DIR, cachePath);
		
	}
	
	/** Is the organization of files within the directory split, as on the PDB FTP servers,
	 * or are all files contained in one directory.
	 * @return flag 
	 */
	public boolean isSplit() {
		return isSplit;
	}

	/** Is the organization of files within the directory split, as on the PDB FTP servers,
	 * or are all files contained in one directory.
	 * @param isSplit flag 
	 */
	public void setSplit(boolean isSplit) {
		this.isSplit = isSplit;
	}

	/** Does the cache automatically download files that are missing from the local installation from the PDB FTP site?
	 * 
	 * @return flag
	 */
	public boolean isAutoFetch() {
		return autoFetch;
	}

	/** Does the cache automatically download files that are missing from the local installation from the PDB FTP site?
	 * 
	 * @param autoFetch flag
	 */
	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}



	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchFileEvenIfObsolete the fetchFileEvenIfObsolete to set
	 */
	public void setFetchFileEvenIfObsolete(boolean fetchFileEvenIfObsolete) {
		this.fetchFileEvenIfObsolete = fetchFileEvenIfObsolete;
	}


	/**forces the cache to fetch the file if its status is OBSOLETE.
	 * This feature has a higher priority than {@link #setFetchCurrent(boolean)}.<br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchFileEvenIfObsolete
	 * @author Amr AL-Hossary
	 * @see #fetchCurrent
	 * @since 3.0.2
	 */
	public boolean isFetchFileEvenIfObsolete() {
		return fetchFileEvenIfObsolete;
	}


	/**if enabled, the reader searches for the newest possible PDB ID, if not present in he local installation.
	 * The {@link #setFetchFileEvenIfObsolete(boolean)} function has a higher priority than this function.<br>
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @param fetchCurrent the fetchCurrent to set
	 * @author Amr AL-Hossary
	 * @see #setFetchFileEvenIfObsolete(boolean)
	 * @since 3.0.2
	 */
	public void setFetchCurrent(boolean fetchNewestCurrent) {
		this.fetchCurrent = fetchNewestCurrent;
	}

	/**
	 * <b>N.B.</b> This feature won't work unless the structure wasn't found & autoFetch is set to <code>true</code>.
	 * @return the fetchCurrent
	 */
	public boolean isFetchCurrent() {
		return fetchCurrent;
	}


	/**
	 * Reports whether strict scop naming will be enforced, or whether this AtomCache
	 * should try to guess some simple variants on scop domains.
	 * @return true if scop names should be used strictly with no guessing
	 */
	public boolean isStrictSCOP() {
		return strictSCOP;
	}

	/**
	 * When strictSCOP is enabled, SCOP domain identifiers (eg 'd1gbga_') are
	 * matched literally to the SCOP database.
	 * 
	 * When disabled, some simple mistakes are corrected automatically.
	 * For instance, the invalid identifier 'd1gbg__' would be corrected to 'd1gbga_' automatically.
	 * @param strictSCOP Indicates whether strict scop names should be used.
	 */
	public void setStrictSCOP(boolean strictSCOP) {
		this.strictSCOP = strictSCOP;
	}

	/** Returns the representation of a ScopDomain as a BioJava Structure object
	 * 
	 * @param domain a scop domain
	 * @return a Structure object.
	 * @throws IOException
	 * @throws StructureException
	 */

	public Structure getStructureForDomain(ScopDomain domain) throws IOException, StructureException{
		if ( scopInstallation == null) {
			scopInstallation = ScopFactory.getSCOP();
		}

		Structure s = null;

		String pdbId = domain.getPdbId();

		try {
			s = getStructure(pdbId);

		} catch (StructureException ex){
			System.err.println("error getting Structure for " + pdbId);

			throw new StructureException(ex);
		}


		StringWriter range = new StringWriter();
		range.append("(");
		int rangePos = 0;
		for ( String r : domain.getRanges()) {
			rangePos++;
			range.append(r);
			if ( ( domain.getRanges().size()> 1) && (rangePos < domain.getRanges().size())){
				range.append(",");
			}

		}
		range.append(")");
		//System.out.println("getting range for "+ pdbId + " " + range);

		Structure n = StructureTools.getSubRanges(s, range.toString());


		// get ligands of chain
		// add the ligands of the chain...
		StructureName structureName = new StructureName(domain.getScopId());
		
		// this does not work for multi-chain domains as of yet...
		// ignore multi-chain domains
		if (! structureName.getChainId().equals(".")){
			Chain newChain  = n.getChainByPDB(structureName.getChainId());
			Chain origChain = s.getChainByPDB(structureName.getChainId());
			List<Group> ligands = origChain.getAtomLigands();

			for(Group g: ligands){
				if ( ! newChain.getAtomGroups().contains(g)) {
					newChain.addGroup(g);
				}
			}
		}

		n.setName(domain.getScopId());
		n.setPDBCode(domain.getScopId());


		StringWriter d = new StringWriter();
		d.append(domain.getClassificationId());



		if ( scopInstallation != null) {
			int sf = domain.getSuperfamilyId();

			ScopDescription scopDesc = scopInstallation.getScopDescriptionBySunid(sf);

			if( scopDesc != null){
				d.append( " | " );
				d.append(scopDesc.getDescription());

			}
		}
		n.getPDBHeader().setDescription(d.toString());
		return n;
	}


	/** Returns the CA atoms for the provided name. See {@link #getStructure(String)} for supported naming conventions.
	 * 
	 * @param name
	 * @return an array of Atoms. 
	 * @throws IOException
	 * @throws StructureException
	 */
	public  Atom[] getAtoms(String name) throws IOException,StructureException{

		Atom[] atoms = null;

		//System.out.println("loading " + name);
		Structure s = null;
		try {

			s = getStructure(name);

		} catch (StructureException ex){
			System.err.println("error getting Structure for " + name);
			throw new StructureException(ex.getMessage(),ex);
		}

		atoms =  StructureTools.getAtomCAArray(s);

		/*synchronized (cache){
			cache.put(name, atoms);
		}*/


		return atoms;
	}


	/** Returns the CA atoms for the provided name. See {@link #getStructure(String)} for supported naming conventions.
	 * 
	 * @param name
	 * @param clone flag to make  sure that the atoms are getting coned
	 * @return an array of Atoms. 
	 * @throws IOException
	 * @throws StructureException
	 * @deprecated does the same as {@link #getAtoms(String)} ;
	 */
	public  Atom[] getAtoms(String name, boolean clone)throws IOException,StructureException{
		Atom[] atoms =  getAtoms(name);

		if ( clone)
			return StructureTools.cloneCAArray(atoms);
		return atoms; 

	}




	/** Request a Structure based on a <i>name</i>.
	 * 
	 * <pre>
		Formal specification for how to specify the <i>name</i>:

		name     := pdbID
		               | pdbID '.' chainID
		               | pdbID '.' range
		               | scopID
		range         := '('? range (',' range)? ')'?
		               | chainID
		               | chainID '_' resNum '-' resNum
		pdbID         := [0-9][a-zA-Z0-9]{3}
		chainID       := [a-zA-Z0-9]
		scopID        := 'd' pdbID [a-z_][0-9_]
		resNum        := [-+]?[0-9]+[A-Za-z]?


		Example structures:
		1TIM     #whole structure
		4HHB.C     #single chain
		4GCR.A_1-83     #one domain, by residue number
		3AA0.A,B     #two chains treated as one structure
		d2bq6a1     #scop domain
		</pre>
	 * 
	 * With the additional set of rules:
	 *  
	 *  <ul>
	 *  <li>If only a PDB code is provided, the whole structure will be return including ligands, but the first model only (for NMR).
	 *	<li>Chain IDs are case sensitive, PDB ids are not. To specify a particular chain write as: 4hhb.A or 4HHB.A </li>
	 *  <li>To specify a SCOP domain write a scopId e.g. d2bq6a1. Some flexibility can be allowed in SCOP domain names, see {@link #setStrictSCOP(boolean)}</li>
	 *  <li>URLs are accepted as well</li>
	 *  </ul>
	 *  
	 * @param name
	 * @return a Structure object, or null if name appears improperly formated (eg too short, etc)
	 * @throws IOException The PDB file cannot be cached due to IO errors
	 * @throws StructureException The name appeared valid but did not correspond to a structure.
	 * 	Also thrown by some submethods upon errors, eg for poorly formatted subranges.
	 */
	public Structure getStructure(String name) throws IOException, StructureException{

		if ( name.length() < 4)
			throw new IllegalArgumentException("Can't interpret IDs that are shorter than 4 residues!");

		Structure n = null;

		boolean useChainNr = false;
		boolean useDomainInfo = false;
		String range = null;
		int chainNr = -1;

		try {

			StructureName structureName = new StructureName(name);

			String pdbId   = null;
			String chainId = null;

			if ( name.length() == 4){

				pdbId = name; 
				Structure s = loadStructureFromByPdbId(pdbId);
				return s;
			} else if ( structureName.isScopName()){

				// return based on SCOP domain ID
				return getStructureFromSCOPDomain(name);
			} else if ( structureName.isCathID()) {
				return getStructureFromCATHDomain(structureName);
			} else if (name.length() == 6){
				// name is PDB.CHAINID style (e.g. 4hhb.A)

				pdbId = name.substring(0,4);
				if ( name.substring(4,5).equals(CHAIN_SPLIT_SYMBOL)) {
					chainId = name.substring(5,6);
				} else if ( name.substring(4,5).equals(CHAIN_NR_SYMBOL)) {

					useChainNr = true;	
					chainNr = Integer.parseInt(name.substring(5,6));
				}



			}  else if ( name.startsWith("file:/") || name.startsWith("http:/") ) {

				// this is a URL
				try {

					URL url = new URL(name);

					return getStructureFromURL(url);

				} catch (Exception e){
					e.printStackTrace();
					return null;
				}


			} else if ( structureName.isPDPDomain()){

				// this is a PDP domain definition

				try {
					return getPDPStructure(name);
				} catch (Exception e){
					e.printStackTrace();
					return null;
				}
			} else if ( name.startsWith(BIOL_ASSEMBLY_IDENTIFIER)){

				try {

					return getBioAssembly(name);

				} catch (Exception e){
					e.printStackTrace();
					return null;
				}
			} else if ( (name.length() > 6) 
					&&  ( ! name.startsWith(PDP_DOMAIN_IDENTIFIER) ) 
					&& (name.contains(CHAIN_NR_SYMBOL) || name.contains(UNDERSCORE)) 
					&& ( ! (name.startsWith("file:/") || name.startsWith("http:/")))

					) {

				// this is a name + range 

				pdbId = name.substring(0,4);
				// this ID has domain split information...
				useDomainInfo = true;
				range = name.substring(5);

			}

			//System.out.println("got: >" + name + "< " + pdbId + " " + chainId + " useChainNr:" + useChainNr + " " +chainNr + " useDomainInfo:" + useDomainInfo + " " + range);

			if (pdbId == null) {

				return null;
			}

			while ( checkLoading(pdbId) ){
				// waiting for loading to be finished...

				try {
					Thread.sleep(100);
				} catch (InterruptedException e){
					System.err.println(e.getMessage());
				}

			}


			//long start  = System.currentTimeMillis();

			Structure s = loadStructureFromByPdbId(pdbId);

			//long end  = System.currentTimeMillis();
			//System.out.println("time to load " + pdbId + " " + (end-start) + "\t  size :" + StructureTools.getNrAtoms(s) + "\t cached: " + cache.size());
			
			if ( chainId == null && chainNr < 0 && range == null) {								
				// we only want the 1st model in this case
				n = StructureTools.getReducedStructure(s,-1);
			}
			else {

				if ( useChainNr) {
					//System.out.println("using ChainNr");
					n = StructureTools.getReducedStructure(s, chainNr);
				} else if ( useDomainInfo) {
					//System.out.println("calling getSubRanges");
					n = StructureTools.getSubRanges(s, range);
				} else  {
					//System.out.println("reducing Chain Id " + chainId);
					n = StructureTools.getReducedStructure(s, chainId);
				}
			}


		} catch (Exception e){
			System.err.println("problem loading:" + name);
			e.printStackTrace();

			throw new StructureException(e.getMessage() + " while parsing " + name,e);

		}

		n.setName(name);
		return n;


	}

	protected Structure loadStructureFromByPdbId(String pdbId)
			throws StructureException {
		
		
		Structure s;
		flagLoading(pdbId);
		try {
			PDBFileReader reader = new PDBFileReader();
			reader.setPath(path);
			reader.setPdbDirectorySplit(isSplit);
			reader.setAutoFetch(autoFetch);
			reader.setFetchFileEvenIfObsolete(fetchFileEvenIfObsolete);
			reader.setFetchCurrent(fetchCurrent);

			reader.setFileParsingParameters(params);

			s = reader.getStructureById(pdbId.toLowerCase());

		} catch (Exception e){
			flagLoadingFinished(pdbId);
			throw new StructureException(e.getMessage() + " while parsing " + pdbId,e);
		}
		flagLoadingFinished(pdbId);
		
		return s;
	}


	private Structure getBioAssembly(String name) throws IOException, StructureException {

		// can be specified as:
		// BIO:1fah   - first one
		// BIO:1fah:0 - asym unit
		// BIO:1fah:1 - first one
		// BIO:1fah:2 - second one 

		String pdbId = name.substring(4,8);
		int biolNr = 1;
		if ( name.length() > 8)
			biolNr = Integer.parseInt(name.substring(9,name.length()));

		return StructureIO.getBiologicalAssembly(pdbId,biolNr);

	}

	private Structure getPDPStructure(String pdpDomainName) {

		//System.out.println("loading PDP domain from server " + pdpDomainName);
		if (pdpprovider == null)
			pdpprovider = new RemotePDPProvider(true);

		return pdpprovider.getDomain(pdpDomainName,this);

	}

	private Structure getStructureFromSCOPDomain(String name)
			throws IOException, StructureException {
		// looks like a SCOP domain!
		ScopDomain domain;
		if( this.strictSCOP) {
			domain = getScopDomain(name);
		} else {
			domain = guessScopDomain(name);
		}

		System.out.println(domain);
		if ( domain != null){
			Structure s = getStructureForDomain(domain);
			return s;
		}

		// Guessing didn't work, so just use the PDBID and Chain from name
		if( !this.strictSCOP) {
			Matcher scopMatch = scopIDregex.matcher(name);
			if( scopMatch.matches() ) {
				String pdbID = scopMatch.group(1);
				String chainID = scopMatch.group(2);

				// None of the actual SCOP domains match. Guess that '_' means 'whole chain'
				if( !chainID.equals("_") ) {
					//Add chain identifier
					pdbID += "."+scopMatch.group(2);
				}
				// Fetch the structure by pdb id
				Structure struct = getStructure(pdbID);
				if(struct != null) {
					System.err.println("Trying chain "+pdbID);
				}

				return struct;
			}
		}

		throw new StructureException("Unable to get structure for SCOP domain: "+name);
	}


	private Structure getStructureFromCATHDomain(StructureName structureName)
			throws IOException, StructureException {



		CathInstallation cathInstall = new CathInstallation(path);

		CathDomain cathDomain = cathInstall.getDomainByCathId(structureName.getName());

		List<CathSegment> segments = cathDomain.getSegments();

		StringWriter range = new StringWriter();

		int rangePos = 0;
		String chainId = structureName.getChainId();
		for ( CathSegment segment : segments) {
			rangePos++;

			range.append(chainId);
			range.append("_");

			range.append(segment.getStart());
			range.append("-");
			range.append(segment.getStop());
			if ( ( segments.size()> 1) && (rangePos < segments.size())){
				range.append(",");
			}
		}

		String pdbId = structureName.getPdbId();

		Structure s = null;

		try {
			s = getStructure(pdbId);

		} catch (StructureException ex){
			System.err.println("error getting Structure for " + pdbId);

			throw new StructureException(ex);
		}

		String rangeS = range.toString();
		System.out.println(rangeS);
		Structure n = StructureTools.getSubRanges(s, rangeS);

		// add the ligands of the chain...

		Chain newChain = n.getChainByPDB(structureName.getChainId());
		Chain origChain = s.getChainByPDB(structureName.getChainId());
		List<Group> ligands = origChain.getAtomLigands();

		for(Group g: ligands){
			if ( ! newChain.getAtomGroups().contains(g)) {
				newChain.addGroup(g);
			}
		}

		// set new Header..
		n.setName(structureName.getName());
		n.setPDBCode(structureName.getPdbId());

		n.getPDBHeader().setDescription(cathDomain.getDomainName());

		return n;
	}


	private Structure getStructureFromURL(URL url) throws IOException, StructureException {
		// looks like a URL for a file was provided:
		System.out.println("fetching structure from URL:" + url);

		String queryS = url.getQuery();

		String chainId = null;
		if ( queryS != null && (queryS.startsWith("chainId="))) {
			chainId = queryS.substring(8);

			String fullu = url.toString();

			if (fullu.startsWith("file:") && fullu.endsWith("?"+queryS)) {
				// for windowze, drop the query part from the URL again
				// otherwise there will be a "file not found error" ...

				String newu = fullu.substring(0,(fullu.length()-(("?"+queryS).length())));
				//System.out.println(newu);
				url = new URL(newu);
			}
		}


		PDBFileReader reader = new PDBFileReader();
		reader.setPath(path);
		reader.setPdbDirectorySplit(isSplit);
		reader.setAutoFetch(autoFetch);
		reader.setFetchFileEvenIfObsolete(fetchFileEvenIfObsolete);
		reader.setFetchCurrent(fetchCurrent);

		reader.setFileParsingParameters(params);

		Structure s = reader.getStructure(url);
		if ( chainId == null)
			return StructureTools.getReducedStructure(s,-1);
		else 
			return StructureTools.getReducedStructure(s,chainId);
	}


	public static final Pattern scopIDregex = Pattern.compile("d(....)(.)(.)" );
	/**
	 * <p>Guess a scop domain. If an exact match is found, return that.
	 * 
	 * <p>Otherwise, return the first scop domain found for the specified protein
	 * such that<ul>
	 *   <li>The chains match, or one of the chains is '_' or '.'.
	 *   <li>The domains match, or one of the domains is '_'.
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
		if ( domain != null){
			return domain;
		}

		// Didn't work. Guess it!
		System.err.println("Warning, could not find SCOP domain: " + name);

		Matcher scopMatch = scopIDregex.matcher(name);
		if( scopMatch.matches() ) {
			String pdbID = scopMatch.group(1);
			String chainID = scopMatch.group(2);
			String domainID = scopMatch.group(3);

			if ( scopInstallation == null) {
				scopInstallation = ScopFactory.getSCOP();
			}

			for( ScopDomain potentialSCOP : scopInstallation.getDomainsForPDB(pdbID) ) {
				Matcher potMatch = scopIDregex.matcher(potentialSCOP.getScopId());
				if(potMatch.matches()) {
					if( chainID.equals(potMatch.group(2)) ||
							chainID.equals("_") || chainID.equals(".") ||
							potMatch.group(2).equals("_") || potMatch.group(2).equals(".") ) {
						if( domainID.equals(potMatch.group(3)) || domainID.equals("_") || potMatch.group(3).equals("_") ) {
							// Match, or near match
							matches.add(potentialSCOP);
						}
					}
				}
			}
		}

		Iterator<ScopDomain> match = matches.iterator();
		if( match.hasNext() ) {
			ScopDomain bestMatch = match.next();
			System.err.print("Trying domain "+bestMatch.getScopId()+".");
			if( match.hasNext() ) {
				System.err.print(" Other possibilities: ");
				while(match.hasNext()) {
					System.err.print(match.next().getScopId() + " ");
				}
			}
			System.err.println();
			return bestMatch;
		} else {
			return null;
		}
	}

	private  boolean checkLoading(String name) {
		return  currentlyLoading.contains(name);

	}

	protected  void flagLoading(String name){
		if ( ! currentlyLoading.contains(name))	
			currentlyLoading.add(name);
	}

	protected  void flagLoadingFinished(String name){
		currentlyLoading.remove(name);   
	}

	private ScopDomain getScopDomain(String scopId)
	{

		if ( scopInstallation == null) {
			scopInstallation = ScopFactory.getSCOP();
		}

		return scopInstallation.getDomainByScopID(scopId);
	}
	public ScopDatabase getScopInstallation() {
		if ( scopInstallation == null) {
			scopInstallation = ScopFactory.getSCOP();
		}

		return scopInstallation;
	}

	public FileParsingParameters getFileParsingParams()
	{
		return params;
	}

	public void setFileParsingParams(FileParsingParameters params)
	{
		this.params = params;
	}


	/** 
	 * Loads the default biological unit (*.pdb1.gz) file. If it is not available, the original
	 * PDB file will be loaded, i.e., for NMR structures, where the original files is also the 
	 * biological assembly.
	 * 
	 * @param pdbId the PDB ID
	 * @return a structure object
	 * @throws IOException 
	 * @throws StructureException 
	 * @since 3.2
	 */
	public Structure getBiologicalUnit(String pdbId) throws StructureException, IOException{
		int bioAssemblyId = 1;
		boolean bioAssemblyFallback = true;
		return getBiologicalAssembly(pdbId, bioAssemblyId, bioAssemblyFallback);
	}

	/** 
	 * Loads the biological assembly for a given PDB ID and bioAssemblyId. 
	 * If a bioAssemblyId > 0 is specified, the corresponding biological assembly file will be loaded. Note, the
	 * number of available biological unit files varies. Many entries don't have a biological assembly specified (i.e. NMR structures),
	 * many entries have only one biological assembly (bioAssemblyId=1), and a few structures have multiple biological assemblies.
	 * Set bioAssemblyFallback to true, to download the original PDB file in cases that a biological assembly file is not available.
	 * 
	 * @param pdbId the PDB ID
	 * @param bioAssemblyId the ID of the biological assembly
	 * @param bioAssemblyFallback if true, try reading original PDB file in case the biological assembly file is not available
	 * @return a structure object
	 * @throws IOException 
	 * @throws StructureException 
	 * @author Peter Rose
	 * @since 3.2
	 */
	public Structure getBiologicalAssembly(String pdbId, int bioAssemblyId, boolean bioAssemblyFallback) throws StructureException, IOException {
		Structure s;
		if (bioAssemblyId < 1) {
			throw new StructureException("bioAssemblyID must be greater than zero: " + pdbId + 
					" bioAssemblyId " + bioAssemblyId);
		}
		PDBFileReader reader = new PDBFileReader();
		reader.setPath(path);
		reader.setPdbDirectorySplit(isSplit);
		reader.setAutoFetch(autoFetch);
		reader.setFetchFileEvenIfObsolete(fetchFileEvenIfObsolete);
		reader.setFetchCurrent(fetchCurrent);
		reader.setFileParsingParameters(params);
		reader.setBioAssemblyId(bioAssemblyId);
		reader.setBioAssemblyFallback(bioAssemblyFallback);
		s = reader.getStructureById(pdbId.toLowerCase());
		s.setPDBCode(pdbId);
		return s;
	}



	/** Send a signal to the cache that the system is shutting down.
	 * Notifies underlying SerializableCache instances to flush themselves...
	 */
	public void notifyShutdown(){
		//System.out.println(" AtomCache got notify shutdown..");
		if ( pdpprovider != null ){
			if ( pdpprovider instanceof RemotePDPProvider){
				RemotePDPProvider remotePDP = (RemotePDPProvider) pdpprovider;
				remotePDP.flushCache();
			}
		}


		// todo: use a SCOP implementation that is backed by SerializableCache 
		if ( scopInstallation != null ){
			if ( scopInstallation instanceof CachedRemoteScopInstallation ){
				CachedRemoteScopInstallation cacheScop = (CachedRemoteScopInstallation) scopInstallation;			
				cacheScop.flushCache();
			} 
		} 

	}

	public PDPProvider getPdpprovider() {
		return pdpprovider;
	}

	public void setPdpprovider(PDPProvider pdpprovider) {
		this.pdpprovider = pdpprovider;
	}


}
