package org.biojava.bio.structure.scop;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.utils.io.InputStreamProvider;

public class ScopInstallation {

	public static final String DEFAULT_VERSION = "1.75";

	public static final String fileName = "dir.cla.scop.txt_";

	public static final String SCOP_DOWNLOAD = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/";



	public static final String NEWLINE;
	public static final String FILESPLIT ;

	static {

		NEWLINE     = System.getProperty("line.separator");
		FILESPLIT   = System.getProperty("file.separator");
	}

	String cacheLocation ;

	String scopVersion ;

	AtomicBoolean installed;

	
	Map<String, List<ScopDomain>> domainMap;
	/** Create a new SCOP installation. 
	 * 
	 * @param cacheLocation where the SCOP files are stored. If they can't be found at that location they will get automatically downloaded and installed there. 
	 */
	public ScopInstallation(String cacheLocation){

		setCacheLocation(cacheLocation);

		installed = new AtomicBoolean();
		installed.set(false);

		scopVersion = DEFAULT_VERSION;

		domainMap = new HashMap<String, List<ScopDomain>>();
	}

	public void ensureInstalled(){
		if ( installed.get())
			return;

		if ( ! filesAvailable()){
			try {
				downloadSCOPFiles();
			} catch (Exception e){
				e.printStackTrace();
				installed.set(false);
				return;
			}			
		} 
		try {
			parseClassification();
		} catch (Exception e){
			e.printStackTrace();
			installed.set(false);
			return;
		}
		installed.set(true);

	}


	public  List<ScopDomain> getDomainsForPDB(String pdbId){
		ensureInstalled();

		return domainMap.get(pdbId.toLowerCase());
	}


	private void parseClassification() throws IOException{

		File file = new File(getFilename());


		InputStreamProvider ips = new InputStreamProvider();
		BufferedReader buffer = new BufferedReader (new InputStreamReader(ips.getInputStream(file)));

		parseClassification(buffer);

	}




	private void parseClassification(BufferedReader buffer) throws IOException {
		String line = null;


		while ((line = buffer.readLine ()) != null) {
			if ( line.startsWith("#"))
				continue;

			String[] spl  = line.split("\t");
			
			if ( spl.length != 6){
				System.err.println("Can't parse line " + line);
				continue;
				
			}

			String scopId = spl[0];
			String pdbId = spl[1];
			String range = spl[2];
			String classificationId = spl[3];
			String sunid = spl[4];
			String tree = spl[5];
			
			String[] rangeSpl = range.split(",");
			
			
			ScopDomain d = new ScopDomain();
			d.setScopId(scopId);
			d.setPdbId(pdbId);
			
			d.setRanges(Arrays.asList(rangeSpl));
			
			d.setClassificationId(classificationId);
			d.setSunid(sunid);
			
			String[] treeSplit = tree.split(",");
			
			if (  treeSplit.length != 7 ) {
				System.err.println("can't process: " + tree );
			}
			
			int classId =Integer.parseInt(treeSplit[0].substring(3));
			int foldId = Integer.parseInt(treeSplit[1].substring(3));
			int familyId = Integer.parseInt(treeSplit[2].substring(3));
			int superfamilyId = Integer.parseInt(treeSplit[3].substring(3));			
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

	}

	private void downloadSCOPFiles() throws FileNotFoundException, IOException{
		String remoteFilename = fileName + scopVersion;
		URL url = new URL(SCOP_DOWNLOAD + remoteFilename);

		String localFileName = getFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}

	private void downloadFileFromRemote(URL remoteURL, File localFile) throws FileNotFoundException, IOException{
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

	private boolean filesAvailable(){
		String fileName = getFilename();

		File f = new File(fileName);

		return f.exists();
	}

	private String getFilename(){
		String f = cacheLocation + fileName + scopVersion;
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
	public String getScopVersion() {
		return scopVersion;
	}
	public void setScopVersion(String scopVersion) {
		this.scopVersion = scopVersion;
	}

	/** get a ScopDomain by its SCOP ID (warning, they are not stable between releasese!)
	 * 
	 *
	 * @param string e.g. d2bq6a1
	 * @return a ScopDomain or null if no domain with the particular ID could be found
	 */
	public ScopDomain getDomainByScopID(String scopId) {
		
		if ( scopId.length() < 6) {
			throw new IllegalArgumentException("Does not look like a scop ID! " + scopId);
		}
		String pdbId = scopId.substring(1,5);
		List<ScopDomain> doms = getDomainsForPDB(pdbId);
		for ( ScopDomain d : doms){
			if ( d.getScopId().equalsIgnoreCase(scopId)) 
				return d;
		}
		
		return null;
	}


}
