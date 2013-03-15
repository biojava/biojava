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

package org.biojava.bio.structure.scop;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava3.core.util.InputStreamProvider;


/** This class provides access to the SCOP protein structure classification.
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
public class ScopInstallation implements ScopDatabase {

	public static final String DEFAULT_VERSION = "1.75";

	protected String scopVersion;
	
	public static final String claFileName = "dir.cla.scop.txt_";
	public static final String desFileName = "dir.des.scop.txt_";
	public static final String hieFileName = "dir.hie.scop.txt_";

	public static final String SCOP_DOWNLOAD = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/";

	protected String scopDownloadURL ;
	
	public static final String NEWLINE;
	public static final String FILESPLIT ;

	static {

		NEWLINE     = System.getProperty("line.separator");
		FILESPLIT   = System.getProperty("file.separator");
	}

	String cacheLocation ;

	AtomicBoolean installedCla;
	AtomicBoolean installedDes;
	AtomicBoolean installedHie;

	Map<String, List<ScopDomain>> domainMap;
	Map<Integer, ScopDescription> sunidMap;
	Map<Integer, ScopNode> scopTree;

	/** Create a new SCOP installation. 
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

		scopVersion = DEFAULT_VERSION;
		scopDownloadURL = SCOP_DOWNLOAD;
		
		
		domainMap = new HashMap<String, List<ScopDomain>>();

		sunidMap  = new HashMap<Integer, ScopDescription>();
		scopTree  = new TreeMap<Integer, ScopNode>();

	}

	/**
	 * Create a new SCOP installation, downloading the file to "the right place".
	 * This will first check for system properties or environmental variables
	 * called 'CACHE_PDB_DIR', or else will use a temporary directory
	 */
	public ScopInstallation() {
		this((new UserConfiguration()).getCacheFilePath());
	}
	public void ensureClaInstalled(){
		if ( installedCla.get())
			return;

		if ( ! claFileAvailable()){
			try {
				downloadClaFile();
			} catch (Exception e){
				e.printStackTrace();
				installedCla.set(false);
				return;
			}			
		} 

		try {
			parseClassification();

		} catch (Exception e){
			e.printStackTrace();
			installedCla.set(false);
			return;
		}
		installedCla.set(true);

	}

	public void ensureDesInstalled(){
		if ( installedDes.get()) {
			return;
		}

		if ( ! desFileAvailable()){
			try {
				downloadDesFile();
			} catch (Exception e){
				e.printStackTrace();
				installedDes.set(false);
				return;
			}   
		}

		try {

			parseDescriptions();
		} catch (Exception e){
			e.printStackTrace();
			installedDes.set(false);
			return;
		}
		installedDes.set(true);


	}

	public void ensureHieInstalled(){
		if ( installedHie.get()) {
			return;
		}

		if ( ! hieFileAvailable()){
			try {
				downloadHieFile();
			} catch (Exception e){
				e.printStackTrace();
				installedHie.set(false);
				return;
			}   
		}

		try {

			parseHierarchy();
		} catch (Exception e){
			e.printStackTrace();
			installedHie.set(false);
			return;
		}
		installedHie.set(true);


	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getByCategory(org.biojava.bio.structure.scop.ScopCategory)
	 */
	public List<ScopDescription> getByCategory(ScopCategory category){

		ensureDesInstalled();

		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);
			if ( sc.getCategory().equals(category))
				try {
					matches.add((ScopDescription)sc.clone());
				} catch (CloneNotSupportedException e){
					e.printStackTrace();
				}
		}
		return matches;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#filterByClassificationId(java.lang.String)
	 */
	public List<ScopDescription> filterByClassificationId(String query){
		ensureDesInstalled();

		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);


			if( sc.getClassificationId().startsWith(query)){
				matches.add(sc);
				continue;
			}
		}

		return matches;
	}


	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getTree(org.biojava.bio.structure.scop.ScopDomain)
	 */
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
	 * @see org.biojava.bio.structure.scop.ScopDatabase#filterByDomainName(java.lang.String)
	 */
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
	 * @see org.biojava.bio.structure.scop.ScopDatabase#filterByDescription(java.lang.String)
	 */
	public List<ScopDescription> filterByDescription(String query){
		ensureDesInstalled();

		query = query.toLowerCase();
		List<ScopDescription> matches = new ArrayList<ScopDescription>();
		for (Integer i : sunidMap.keySet()){
			ScopDescription sc = sunidMap.get(i);


			if( sc.getDescription().toLowerCase().startsWith(query)){
				matches.add(sc);
				continue;
			}
		}

		return matches;
	}


	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getScopDescriptionBySunid(int)
	 */
	public ScopDescription getScopDescriptionBySunid(int sunid){
		ensureDesInstalled();
		return sunidMap.get(sunid);
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getDomainsForPDB(java.lang.String)
	 */
	public  List<ScopDomain> getDomainsForPDB(String pdbId){
		ensureClaInstalled();


		List<ScopDomain> doms = domainMap.get(pdbId.toLowerCase());
		
		List<ScopDomain> retdoms = new ArrayList<ScopDomain>();
		
		if ( doms == null)
			return retdoms;

		for ( ScopDomain d : doms){
			try {
				ScopDomain n = (ScopDomain) d.clone();
				retdoms.add(n);
			}  catch (CloneNotSupportedException e){
				e.printStackTrace();
			}


		}
		return retdoms;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getDomainByScopID(java.lang.String)
	 */
	public ScopDomain getDomainByScopID(String scopId) {
		ensureClaInstalled();

		if ( scopId.length() < 6) {
			throw new IllegalArgumentException("Does not look like a scop ID! " + scopId);
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
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getScopNode(int)
	 */
	public ScopNode getScopNode(int sunid){
		ensureHieInstalled();
		ScopNode node = scopTree.get(sunid);

		return node;
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
		String line = null;

		int counter =0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 3 ) {
				System.err.println("parseHierarchy: Can't parse line " + line +" (length: " + spl.length+")");
				continue;
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
		System.out.println("parsed " + counter + " scop sunid nodes.");
	}


	private void parseDescriptions() throws IOException{

		File file = new File(getDesFilename());


		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseDescriptions(buffer);

	}
	private void parseDescriptions(BufferedReader buffer) throws IOException {
		String line = null;

		int counter = 0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 5 ) {
				System.err.println("parseDescriptions: Can't parse line " + line +" (length: " + spl.length+")");
				continue;
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

			sunidMap.put(new Integer(sunID), c);

		}
		System.out.println("parsed " + counter + " scop sunid descriptions.");
	}



	private void parseClassification(BufferedReader buffer) throws IOException {
		String line = null;

		int counter = 0;
		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");

			if ( spl.length != 6){
				System.err.println("Can't parse line " + line);
				continue;

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
				System.err.println("can't process: " + line );
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
			if ( sunid == 47763)
				System.out.println("FOUND DOMAIN!!!! " + sunid + " " + d);
		}
		System.out.println("parsed "+ counter + " scop sunid domains.");

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
				if(subRange.charAt(1) == ':') {
					continue;
				}
				else {
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
		String remoteFilename = claFileName + scopVersion;
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getClaFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	protected void downloadDesFile() throws FileNotFoundException, IOException{
		String remoteFilename = desFileName + scopVersion;
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getDesFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	protected void downloadHieFile() throws FileNotFoundException, IOException{
		String remoteFilename = hieFileName + scopVersion;
		URL url = new URL(scopDownloadURL + remoteFilename);

		String localFileName = getHieFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	protected void downloadFileFromRemote(URL remoteURL, File localFile) throws FileNotFoundException, IOException{
		System.out.println("downloading " + remoteURL + " to: " + localFile);
		FileOutputStream out = new FileOutputStream(localFile);

		InputStream in = remoteURL.openStream();
		byte[] buf = new byte[4 * 1024]; // 4K buffer
		int bytesRead;
		while ((bytesRead = in.read(buf)) != -1) {
			out.write(buf, 0, bytesRead);
		}
		in.close();
		out.close();


	}

	private boolean claFileAvailable(){
		String fileName = getClaFilename();

		File f = new File(fileName);

		return f.exists();
	}

	private boolean desFileAvailable(){
		String fileName = getDesFilename();

		File f = new File(fileName);

		return f.exists();
	}

	private boolean hieFileAvailable(){
		String fileName = getHieFilename();

		File f = new File(fileName);

		return f.exists();
	}

	protected String getClaFilename(){
		String f = cacheLocation + claFileName + scopVersion;
		return f;
	}

	protected String getDesFilename(){
		String f = cacheLocation + desFileName + scopVersion;
		return f;

	}

	protected String getHieFilename(){
		String f = cacheLocation + hieFileName + scopVersion;
		return f;

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
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getScopVersion()
	 */
	public String getScopVersion() {
		return scopVersion;
	}
	public void setScopVersion(String scopVersion) {
		this.scopVersion = scopVersion;
	}
	

	public String getScopDownloadURL() {
		return scopDownloadURL;
	}

	public void setScopDownloadURL(String scopDownloadURL) {
		this.scopDownloadURL = scopDownloadURL;
	}

	/* (non-Javadoc)
	 * @see org.biojava.bio.structure.scop.ScopDatabase#getScopDomainsBySunid(java.lang.Integer)
	 */
	public List<ScopDomain> getScopDomainsBySunid(Integer sunid)
	{

		ensureClaInstalled();

		List<ScopDomain> domains = new ArrayList<ScopDomain>();

		for (String pdbId: domainMap.keySet()){
			for (ScopDomain d : domainMap.get(pdbId)){
				try {
					if ( d.getPx() == sunid) {
						domains.add((ScopDomain)d.clone());
						continue;
					} else if ( d.getSpeciesId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}else if ( d.getDomainId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}else if ( d.getFamilyId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}else if ( d.getSuperfamilyId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}else if ( d.getFoldId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}else if ( d.getClassId() == sunid ){
						domains.add((ScopDomain)d.clone());
						continue;
					}
				} catch (CloneNotSupportedException e){
					e.printStackTrace();
				}
			}
		}
		return domains;

	}
}
