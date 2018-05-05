/**
 *                  BioJava development code
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
 * Author: Andreas Prlic
 *
 *
 */

package org.biojava.nbio.structure.scop;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.core.util.FileDownloadUtils;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;


/**
 * This class provides access to the SCOP protein structure classification.
 *
 * For more information about SCOP see here:
 *  <ul>
 *   <li>SCOP: <a href="http://scop.mrc-lmb.cam.ac.uk/scop/">http://scop.mrc-lmb.cam.ac.uk/scop/</a></li>

 *  <li> Introduction: <a href="http://scop.mrc-lmb.cam.ac.uk/scop/intro.html">http://scop.mrc-lmb.cam.ac.uk/scop/intro.html</a> </li>

 *   <li> SCOP parsable files: <a href="http://scop.mrc-lmb.cam.ac.uk/scop/parse/">http://scop.mrc-lmb.cam.ac.uk/scop/parse/</a> </li>
 * </ul>

 *
 * This class can automatically download missing files from the SCOP classification.
 *
 * @author Andreas Prlic
 *
 */
public class ScopInstallation implements LocalScopDatabase {

	public static final String DEFAULT_VERSION = "1.75";

	private static final Logger logger = LoggerFactory.getLogger(ScopInstallation.class);

	private String scopVersion;

	// Stores URLs for cla, des, hie, and com files
	private final List<ScopMirror> mirrors;

	// Cache filenames (with version appended)
	public static final String claFileName = "dir.cla.scop.txt_";
	public static final String desFileName = "dir.des.scop.txt_";
	public static final String hieFileName = "dir.hie.scop.txt_";
	public static final String comFileName = "dir.com.scop.txt_";

	// Download locations
	public static final String SCOP_DOWNLOAD = "http://scop.berkeley.edu/downloads/parse/";
	public static final String SCOP_DOWNLOAD_ALTERNATE = "http://scop.berkeley.edu/downloads/parse/";

	//public static final String NEWLINE = System.getProperty("line.separator");
	public static final String FILESPLIT = System.getProperty("file.separator");

	private String cacheLocation ;

	private AtomicBoolean installedCla;
	private AtomicBoolean installedDes;
	private AtomicBoolean installedHie;
	private AtomicBoolean installedCom;

	private Map<Integer, List<String>> commentsMap;
	private Map<String, List<ScopDomain>> domainMap;
	private Map<Integer, ScopDescription> sunidMap;
	private Map<Integer, ScopNode> scopTree;


	/**
	 * Create a new SCOP installation.
	 *
	 * @param cacheLocation where the SCOP files are stored. If they can't be found at that location they will get automatically downloaded and installed there.
	 */
	public ScopInstallation(String cacheLocation){

		setCacheLocation(cacheLocation);

		installedCla = new AtomicBoolean();
		installedCla.set(false);
		installedDes = new AtomicBoolean();
		installedDes.set(false);
		installedHie = new AtomicBoolean();
		installedHie.set(false);
		installedCom = new AtomicBoolean();
		installedCom.set(false);

		scopVersion = DEFAULT_VERSION;
		mirrors = new ArrayList<ScopMirror>(1);

		domainMap = new HashMap<String, List<ScopDomain>>();

		sunidMap  = new HashMap<Integer, ScopDescription>();
		scopTree  = new TreeMap<Integer, ScopNode>();

	}

	/**
	 * Removes all of the comments (dir.com file) in order to free memory. The file will need to be reloaded if {@link #getComments(int)} is called subsequently.
	 */
	public void nullifyComments() {
		commentsMap = null;
		installedCom.set(false);
	}

	/**
	 * Create a new SCOP installation, downloading the file to "the right place".
	 * This will first check for system properties or environmental variables
	 * called {@link UserConfiguration#PDB_CACHE_DIR}, or else will use a temporary directory
	 */
	public ScopInstallation() {
		this((new UserConfiguration()).getCacheFilePath());
	}

	public void ensureClaInstalled() throws IOException {
		if (installedCla.get()) return;
		if (!claFileAvailable()) downloadClaFile();
		parseClassification();
		installedCla.set(true);
	}

	public void ensureDesInstalled() throws IOException {
		if (installedDes.get()) return;
		if (!desFileAvailable()) downloadDesFile();
		parseDescriptions();
		installedDes.set(true);
	}

	public void ensureComInstalled() throws IOException {
		if (installedCom.get()) return;
		if (!comFileAvailable()) downloadComFile();
		parseComments();
		installedCom.set(true);
	}

	public void ensureHieInstalled() throws IOException {
		if ( installedHie.get()) return;
		if ( ! hieFileAvailable()) downloadHieFile();
		parseHierarchy();
		installedHie.set(true);
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getByCategory(org.biojava.nbio.structure.scop.ScopCategory)
	 */
	@Override
	public List<ScopDescription> getByCategory(ScopCategory category){

		try {
			ensureDesInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);
			if ( sc.getCategory().equals(category))

				try {
					matches.add((ScopDescription)sc.clone());
				} catch (CloneNotSupportedException e) {
					throw new RuntimeException("Could not clone " + ScopDescription.class + " subclass", e);
				}

		}
		return matches;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#filterByClassificationId(java.lang.String)
	 */
	@Override
	public List<ScopDescription> filterByClassificationId(String query){

		try {
			ensureDesInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);


			if( sc.getClassificationId().startsWith(query)){
				matches.add(sc);
			}
		}

		return matches;
	}


	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getTree(org.biojava.nbio.structure.scop.ScopDomain)
	 */
	@Override
	public List<ScopNode> getTree(ScopDomain domain){
		ScopNode node = getScopNode(domain.getSunid());


		List<ScopNode> tree = new ArrayList<ScopNode>();
		while (node != null){

			//System.out.println("This node: sunid:" + node.getSunid() );
			//System.out.println(getScopDescriptionBySunid(node.getSunid()));
			node = getScopNode(node.getParentSunid());
			if ( node != null)
				tree.add(node);
		}
		Collections.reverse(tree);
		return tree;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#filterByDomainName(java.lang.String)
	 */
	@Override
	public List<ScopDomain> filterByDomainName(String query) {

		List<ScopDomain > domains = new ArrayList<ScopDomain>();
		if (query.length() <5){
			return domains;
		}

		String pdbId = query.substring(1,5);

		List<ScopDomain> doms = getDomainsForPDB(pdbId);


		if ( doms == null)
			return domains;

		query = query.toLowerCase();
		for ( ScopDomain d: doms){
			if ( d.getScopId().toLowerCase().contains(query)){
				domains.add(d);
			}
		}

		return domains;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#filterByDescription(java.lang.String)
	 */
	@Override
	public List<ScopDescription> filterByDescription(String query) throws ScopIOException {
		try {
			ensureDesInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		query = query.toLowerCase();
		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);

			if( sc.getDescription().toLowerCase().startsWith(query)){
				matches.add(sc);
			}
		}

		return matches;
	}


	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getScopDescriptionBySunid(int)
	 */
	@Override
	public ScopDescription getScopDescriptionBySunid(int sunid) {
		try {
			ensureDesInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}
		return sunidMap.get(sunid);
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getDomainsForPDB(java.lang.String)
	 */
	@Override
	public  List<ScopDomain> getDomainsForPDB(String pdbId) {

		try {
			ensureClaInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		List<ScopDomain> doms = domainMap.get(pdbId.toLowerCase());

		List<ScopDomain> retdoms = new ArrayList<ScopDomain>();

		if ( doms == null)
			return retdoms;

		for ( ScopDomain d : doms){
			try {
				ScopDomain n = (ScopDomain) d.clone();
				retdoms.add(n);
			}  catch (CloneNotSupportedException e){
				throw new RuntimeException(ScopDomain.class + " subclass does not support clone()", e);
			}


		}
		return retdoms;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getDomainByScopID(java.lang.String)
	 */
	@Override
	public ScopDomain getDomainByScopID(String scopId) {

		try {
			ensureClaInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		if ( scopId.length() < 6) {
			throw new ScopIOException("Does not look like a scop ID! " + scopId);
		}
		String pdbId = scopId.substring(1,5);
		List<ScopDomain> doms = getDomainsForPDB(pdbId);
		if ( doms == null)
			return null;
		for ( ScopDomain d : doms){
			if ( d.getScopId().equalsIgnoreCase(scopId))
				return d;
		}

		return null;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getScopNode(int)
	 */
	@Override
	public ScopNode getScopNode(int sunid){

		try {
			ensureHieInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		return scopTree.get(sunid);
	}


	private void parseClassification() throws IOException{

		File file = new File(getClaFilename());

		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseClassification(buffer);

	}

	private void parseHierarchy() throws IOException{

		File file = new File(getHieFilename());


		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseHierarchy(buffer);

	}

	private void parseHierarchy(BufferedReader buffer) throws IOException {
		String line;

		int counter =0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 3 ) {
				throw new IOException("parseHierarchy: Can't parse line " + line +" (length: " + spl.length+")");
			}
			counter++;
			int sunid       = Integer.parseInt(spl[0]);
			int parentSunid = -1;

			if ( sunid != 0)
				parentSunid = Integer.parseInt(spl[1]);

			String children = spl[2];
			String[] childIds = children.split(",");

			List<Integer> chis = new ArrayList<Integer>();

			for ( String id : childIds){
				if ( id.equals("-"))
					continue;
				chis.add(Integer.parseInt(id));
			}

			ScopNode node = new ScopNode();

			node.setSunid(sunid);
			node.setParentSunid(parentSunid);
			node.setChildren(chis);

			scopTree.put(sunid, node);
		}
		logger.info("Parsed {} SCOP sunid nodes.", counter);
	}


	private void parseDescriptions() throws IOException{

		File file = new File(getDesFilename());

		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseDescriptions(buffer);

	}

	private void parseComments() throws IOException{

		File file = new File(getComFilename());

		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseComments(buffer);

	}

	private void parseComments(BufferedReader buffer) throws IOException {

		commentsMap = new HashMap<Integer,List<String>>();

		int counter = 0;
		String line;
		while ((line = buffer.readLine ()) != null) {
			if (line.startsWith("#")) continue;
			String[] parts = line.split("!");
			int sunId = Integer.parseInt(parts[0].trim());
			if (parts.length == 1) {
				commentsMap.put(sunId, new ArrayList<String>(1));
				continue;
			}
			List<String> comments = new ArrayList<String>(parts.length - 1);
			for (int i = 1; i < parts.length; i++) {
				String trimmed = parts[i].trim();
				if( !trimmed.isEmpty() ) {
					comments.add(trimmed);
				}
			}
			commentsMap.put(sunId, comments);
			counter++;
		}
		logger.info("Parsed {} SCOP comments.", counter);

	}

	private void parseDescriptions(BufferedReader buffer) throws IOException {
		String line = null;

		int counter = 0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 5 ) {
				throw new IOException("parseDescriptions: Can't parse line " + line +" (length: " + spl.length+")");
			}
			counter++;

			//46464  dm  a.1.1.2 -   Hemoglobin I
			int sunID = Integer.parseInt(spl[0]);
			ScopCategory category =  ScopCategory.fromString(spl[1]);
			String classificationId = spl[2];
			String name = spl[3];
			String desc = spl[4];

			ScopDescription c = new ScopDescription();
			c.setSunID(sunID);
			c.setCategory(category);
			c.setClassificationId(classificationId);
			c.setName(name);
			c.setDescription(desc);

			sunidMap.put(sunID, c);

		}
		logger.info("Parsed {} SCOP sunid descriptions.", counter);
	}



	private void parseClassification(BufferedReader buffer) throws IOException {
		String line = null;

		int counter = 0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 6){
				throw new IOException("Can't parse line " + line);
			}
			counter++;

			String scopId = spl[0];
			String pdbId = spl[1];
			String range = spl[2];
			String classificationId = spl[3];
			Integer sunid = Integer.parseInt(spl[4]);
			String tree = spl[5];



			ScopDomain d = new ScopDomain();
			d.setScopId(scopId);
			d.setPdbId(pdbId);

			d.setRanges(extractRanges(range));

			d.setClassificationId(classificationId);
			d.setSunid(sunid);

			String[] treeSplit = tree.split(",");

			if (  treeSplit.length != 7 ) {
				throw new IOException("Can't process: " + line );
			}

			int classId =Integer.parseInt(treeSplit[0].substring(3));
			int foldId = Integer.parseInt(treeSplit[1].substring(3));
			int superfamilyId = Integer.parseInt(treeSplit[2].substring(3));
			int familyId = Integer.parseInt(treeSplit[3].substring(3));
			int domainId = Integer.parseInt(treeSplit[4].substring(3));
			int speciesId = Integer.parseInt(treeSplit[5].substring(3));
			int px = Integer.parseInt(treeSplit[6].substring(3));

			d.setClassId(classId);
			d.setFoldId(foldId);
			d.setSuperfamilyId(superfamilyId);
			d.setFamilyId(familyId);
			d.setDomainId(domainId);
			d.setSpeciesId(speciesId);
			d.setPx(px);

			List<ScopDomain> domainList;
			if ( domainMap.containsKey(pdbId)){
				domainList = domainMap.get(pdbId);
			} else {
				domainList = new ArrayList<ScopDomain>();
				domainMap.put(pdbId,domainList);
			}

			domainList.add(d);
		}
		logger.info("Parsed {} SCOP sunid domains.", counter);

	}

	/**
	 * Converts the SCOP range field into a list of subranges suitable for
	 * storage in a ScopDomain object. Each range should be of a format
	 * compatible with {@link StructureTools#getSubRanges(Structure,String)}.
	 * @param range
	 * @return
	 */
	private List<String> extractRanges(String range) {
		List<String> ranges;
		String[] rangeSpl = range.split(",");

		// Recent versions of scop always specify a chain, so no processing is needed
		if(scopVersion.compareTo("1.73") < 0 ) {
			for(int i=0; i<rangeSpl.length;i++) {
				String subRange = rangeSpl[i];

				// Allow single-chains, as well as the '-' special case
				if(subRange.length()<2) {
					continue;
				}

				// Allow explicit chain syntax
				if(subRange.charAt(1) != ':') {
					// Early versions sometimes skip the chain identifier for single-chain domains
					// Indicate this with a chain "_"
					rangeSpl[i] = "_:"+subRange;
				}
			}
		}
		ranges = Arrays.asList(rangeSpl);
		return ranges;
	}

	protected void downloadClaFile() throws FileNotFoundException, IOException{
		if(mirrors.size()<1) {
			initScopURLs();
		}
		IOException exception = null;
		for(ScopMirror mirror:mirrors) {
			try {
				URL url = new URL(mirror.getClaURL(scopVersion));

				String localFileName = getClaFilename();
				File localFile = new File(localFileName);

				downloadFileFromRemote(url, localFile);
				return;
			} catch(IOException e ) {
				exception = e;
			}
		}
		throw new IOException("Unable to download SCOP .cla file",exception);
	}

	protected void downloadDesFile() throws FileNotFoundException, IOException{
		if(mirrors.size()<1) {
			initScopURLs();
		}
		IOException exception = null;
		for(ScopMirror mirror:mirrors) {
			try {
				URL url = new URL(mirror.getDesURL( scopVersion));

				String localFileName = getDesFilename();
				File localFile = new File(localFileName);

				downloadFileFromRemote(url, localFile);
				return;
			} catch(IOException e ) {
				exception = e;
			}
		}
		throw new IOException("Unable to download SCOP .des file",exception);
	}

	protected void downloadHieFile() throws IOException{
		if(mirrors.size()<1) {
			initScopURLs();
		}
		IOException exception = null;
		for(ScopMirror mirror:mirrors) {
			try {
				URL url = new URL(mirror.getHieURL( scopVersion));

				String localFileName = getHieFilename();
				File localFile = new File(localFileName);

				downloadFileFromRemote(url, localFile);
				return;
			} catch(IOException e ) {
				exception = e;
			}
		}
		throw new IOException("Unable to download SCOP .hie file",exception);

	}

	protected void downloadComFile() throws FileNotFoundException, IOException{
		if(mirrors.size()<1) {
			initScopURLs();
		}
		IOException exception = null;
		for(ScopMirror mirror:mirrors) {
			try {
				URL url = new URL(mirror.getComURL(scopVersion));

				String localFileName = getComFilename();
				File localFile = new File(localFileName);

				downloadFileFromRemote(url, localFile);
				return;
			} catch (IOException e ) {
				exception = e;
			}
		}
		throw new IOException("Unable to download SCOP .com file",exception);
	}

	protected void downloadFileFromRemote(URL remoteURL, File localFile) throws IOException{
		logger.info("Downloading " + remoteURL + " to: " + localFile);
		FileDownloadUtils.downloadFile(remoteURL, localFile);
	}

	private boolean claFileAvailable(){
		String fileName = getClaFilename();

		File f = new File(fileName);

		return f.exists() && f.length()>0;
	}

	private boolean desFileAvailable(){
		String fileName = getDesFilename();

		File f = new File(fileName);
		return f.exists() && f.length()>0;
	}

	private boolean hieFileAvailable(){
		String fileName = getHieFilename();

		File f = new File(fileName);

		return f.exists() && f.length()>0;
	}

	private boolean comFileAvailable(){
		String fileName = getComFilename();

		File f = new File(fileName);

		return f.exists() && f.length()>0;
	}

	protected String getClaFilename(){
		return cacheLocation + claFileName + scopVersion;
	}

	protected String getDesFilename(){
		return cacheLocation + desFileName + scopVersion;

	}

	protected String getHieFilename(){
		return cacheLocation + hieFileName + scopVersion;

	}

	protected String getComFilename(){
		return cacheLocation + comFileName + scopVersion;
	}

	public String getCacheLocation() {
		return cacheLocation;
	}

	public void setCacheLocation(String cacheLocation) {

		if (! cacheLocation.endsWith(FILESPLIT))
			cacheLocation += FILESPLIT;
		this.cacheLocation = cacheLocation;


	}
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getScopVersion()
	 */
	@Override
	public String getScopVersion() {
		return scopVersion;
	}
	@Override
	public void setScopVersion(String scopVersion) {
		if(scopVersion == null)
			throw new NullPointerException("Null scop version");
		if(this.scopVersion.equals(scopVersion))
			return;
		this.scopVersion = scopVersion;
		// reset installation flags
		installedCla.set(false);
		installedDes.set(false);
		installedHie.set(false);
		installedCom.set(false);

	}

	public void addMirror(String scopDownloadURL) {
		mirrors.add(new ScopMirror(scopDownloadURL));
	}
	void addMirror(ScopMirror scopURLs) {
		mirrors.add(scopURLs);
	}
	public List<ScopMirror> getMirrors() {
		if(mirrors.isEmpty()) {
			this.initScopURLs();
		}
		return mirrors;
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.scop.ScopDatabase#getScopDomainsBySunid(java.lang.Integer)
	 */
	@Override
	public List<ScopDomain> getScopDomainsBySunid(Integer sunid)
	{

		try {
			ensureClaInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}

		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		for (String pdbId: domainMap.keySet()){
			for (ScopDomain d : domainMap.get(pdbId)){
				try {
					if ( d.getPx() == sunid) {
						domains.add((ScopDomain)d.clone());
					} else if ( d.getSpeciesId() == sunid ){
						domains.add((ScopDomain)d.clone());
					}else if ( d.getDomainId() == sunid ){
						domains.add((ScopDomain)d.clone());
					}else if ( d.getFamilyId() == sunid ){
						domains.add((ScopDomain)d.clone());
					}else if ( d.getSuperfamilyId() == sunid ){
						domains.add((ScopDomain)d.clone());
					}else if ( d.getFoldId() == sunid ){
						domains.add((ScopDomain)d.clone());
					}else if ( d.getClassId() == sunid ){
						domains.add((ScopDomain)d.clone());
					} else {
						throw new RuntimeException("Type " + d + " not recognized"); // only possible if SCOP changes
					}
				} catch (CloneNotSupportedException e){
					throw new RuntimeException(ScopDomain.class + " subclass does not support clone()", e);
				}
			}
		}
		return domains;

	}

	@Override
	public List<String> getComments(int sunid) {
		try {
			ensureComInstalled();
		} catch (IOException e) {
			throw new ScopIOException(e);
		}
		if (!commentsMap.containsKey(sunid)) return new ArrayList<String>(1);
		return commentsMap.get(sunid);
	}


	private void initScopURLs() {
		if(!this.mirrors.isEmpty()) {
			return;
		}

		// first, try default scop
		ScopMirror primary = new ScopMirror();
		// If unreachable, try alternate Berkeley location
		ScopMirror alt = new ScopMirror(
				SCOP_DOWNLOAD_ALTERNATE,
				"dir.cla.scop.%s.txt","dir.des.scop.%s.txt",
				"dir.hie.scop.%s.txt","dir.com.scop.%s.txt");
		mirrors.add(primary);
		mirrors.add(alt);
	}
}
