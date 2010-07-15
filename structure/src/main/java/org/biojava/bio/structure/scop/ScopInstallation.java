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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicBoolean;

import org.biojava.utils.io.InputStreamProvider;


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
public class ScopInstallation {

	public static final String DEFAULT_VERSION = "1.75";

	public static final String claFileName = "dir.cla.scop.txt_";
	public static final String desFileName = "dir.des.scop.txt_";
	public static final String hieFileName = "dir.hie.scop.txt_";
	
	public static final String SCOP_DOWNLOAD = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/";

	public static final String NEWLINE;
	public static final String FILESPLIT ;

	static {

		NEWLINE     = System.getProperty("line.separator");
		FILESPLIT   = System.getProperty("file.separator");
	}

	String cacheLocation ;

	String scopVersion ;

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

		domainMap = new HashMap<String, List<ScopDomain>>();
		
		sunidMap  = new HashMap<Integer, ScopDescription>();
		scopTree  = new TreeMap<Integer, ScopNode>();
		
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
	
	/** Get all records of a particular classification.
	 * 
	 * @param category e.g. "superfamily"
	 * @return all records of this type
	 */
	public List<ScopDescription> getByCategory(ScopCategory category){
	   
	   ensureDesInstalled();
	   
	   List<ScopDescription> matches = new ArrayList<ScopDescription>();
	   for (Integer i : sunidMap.keySet()){
	      ScopDescription sc = sunidMap.get(i);
	      if ( sc.getCategory().equals(category))
	         matches.add(sc);
	   }
	   return matches;
	}
	
	/** Return the SCOP description for a node in the hierarchy by its "sunid" id.
	 * 
	 * @param sunid
	 * @return a ScopDescription object
	 */
	public ScopDescription getScopDescriptionBySunid(int sunid){
	   ensureDesInstalled();
	   return sunidMap.get(sunid);
	}
	
	/** Get a list of ScopDomains that have been assigned to a PDB ID
	 * 
	 * @param pdbId the PDB entry
	 * @return a list of ScopDomains
	 */
	public  List<ScopDomain> getDomainsForPDB(String pdbId){
		ensureClaInstalled();

		return domainMap.get(pdbId.toLowerCase());
	}
	
	/** get a ScopDomain by its SCOP ID (warning, they are not stable between releasese!)
     * 
     *
     * @param string e.g. d2bq6a1
     * @return a ScopDomain or null if no domain with the particular ID could be found
     */
    public ScopDomain getDomainByScopID(String scopId) {
       ensureClaInstalled();
       
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
    
    /** Access a particular ScopNode. The scopNode then allows to traverse through the scop hierarchy...
     * 
     * @param sunid
     * @return
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
			if ( sunid == 47763)
			   System.out.println("FOUND DOMAIN!!!! " + sunid + " " + d);
		}
		System.out.println("parsed "+ counter + " scop sunid domains.");

	}

	private void downloadClaFile() throws FileNotFoundException, IOException{
		String remoteFilename = claFileName + scopVersion;
		URL url = new URL(SCOP_DOWNLOAD + remoteFilename);

		String localFileName = getClaFilename();
		File localFile = new File(localFileName);

		downloadFileFromRemote(url, localFile);

	}
	
	private void downloadDesFile() throws FileNotFoundException, IOException{
       String remoteFilename = desFileName + scopVersion;
       URL url = new URL(SCOP_DOWNLOAD + remoteFilename);

       String localFileName = getDesFilename();
       File localFile = new File(localFileName);

       downloadFileFromRemote(url, localFile);

   }
	
	private void downloadHieFile() throws FileNotFoundException, IOException{
       String remoteFilename = hieFileName + scopVersion;
       URL url = new URL(SCOP_DOWNLOAD + remoteFilename);

       String localFileName = getHieFilename();
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
	
	private String getClaFilename(){
		String f = cacheLocation + claFileName + scopVersion;
		return f;
	}
	
	private String getDesFilename(){
	   String f = cacheLocation + desFileName + scopVersion;
	   return f;
	   
	}
	
	private String getHieFilename(){
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
	public String getScopVersion() {
		return scopVersion;
	}
	public void setScopVersion(String scopVersion) {
		this.scopVersion = scopVersion;
	}

	/** Get a SCOP domain by its sunid
	 * 
	 * @param sun
	 * @return
	 */
   public List<ScopDomain> getScopDomainsBySunid(Integer sunid)
   {
      
      ensureClaInstalled();
     
      List<ScopDomain> domains = new ArrayList<ScopDomain>();

      for (String pdbId: domainMap.keySet()){
         for (ScopDomain d : domainMap.get(pdbId)){
            if ( d.getPx() == sunid) {
               domains.add(d);
               continue;
            } else if ( d.getSpeciesId() == sunid ){
               domains.add(d);
               continue;
            }else if ( d.getDomainId() == sunid ){
               domains.add(d);
               continue;
            }else if ( d.getFamilyId() == sunid ){
               domains.add(d);
               continue;
            }else if ( d.getSuperfamilyId() == sunid ){
               domains.add(d);
               continue;
            }else if ( d.getFoldId() == sunid ){
               domains.add(d);
               continue;
            }else if ( d.getClassId() == sunid ){
               domains.add(d);
               continue;
            }
         }
      }
      return domains;
     
   }

	


}
