package org.biojava.bio.structure.scop;

import java.util.HashMap;
import java.util.Map;



/**
 * Controls the global ScopDatabase being used.
 * 
 * <p>Defaults to a {@link RemoteScopInstallation}, which is fast for small numbers
 * of queries. For many queries, a {@link ScopInstallation} instance may be faster,
 * since it makes only one network request.
 * 
 * <p>Example: Fetch the structure corresponding to an old version of scop
 * 
 * <pre>
 * ScopInstallation scop = new ScopInstallation();
 * scop.setScopVersion("1.69");
 * ScopFactory.setScopDatabase(scop);
 * AtomCache cache = new AtomCache();
 * cache.setFetchFileEvenIfObsolete(true); //fetch older PDBs
 * cache.setStrictSCOP(false); // correct simple errors in domain names
 * Structure s = cache.getStructure("d3hbia_");
 * @author sbliven
 *
 */
public class ScopFactory {

	//private static ScopDatabase scop = new ScopInstallation();

	// by default we now request SCOP via web services -> much less memory consumption
	private static ScopDatabase scop = new RemoteScopInstallation();

	public static ScopDatabase getSCOP(){
		return scop;
	}

	public static void setScopDatabase(ScopDatabase s){
		System.out.println("Setting ScopDatabase to type: " + s.getClass().getName());
		scop = s;
	}


	// berkeley
	public static final String VERSION_1_75A = "1.75A";	
	public static final String VERSION_1_75B = "1.75B";	


	// original SCOP
	// latest SCOP release from SCOP website = 1.75;
	public static final String VERSION_1_75 = "1.75";
	
	static Map<String,ScopDatabase> versionedScopDBs = new HashMap<String, ScopDatabase>(); 
	
	
	public static ScopDatabase getSCOP(String version, boolean useLocalData){
		if ( useLocalData) {
			if ( version.equalsIgnoreCase(VERSION_1_75A)) {
				return getBerkeley_1_75A();

			}
			else if ( version.equalsIgnoreCase(VERSION_1_75B)) {
				return getBerkeley_1_75B();
			} else if  ( version.equalsIgnoreCase(VERSION_1_75)){
				getScop_1_75();
			} else {
				getScop_1_75();
			}
		} else {
			// should to get proxied via Domain service servlet
			return scop;
		}
		return scop;
	}
	
	
	private static ScopDatabase getScop_1_75() {
		ScopInstallation scop = (ScopInstallation)versionedScopDBs.get(VERSION_1_75);
		if ( scop == null) {
			scop = new ScopInstallation();
			scop.setScopVersion(VERSION_1_75);
			versionedScopDBs.put(VERSION_1_75, scop);
		}
		
		return scop;
		
	}

	/** requests a particular version of SCOP
	 * 
	 * @param version
	 * @return
	 */
	public static ScopDatabase getSCOP(String version){
		if ( version.equalsIgnoreCase(VERSION_1_75A)) {
			return getBerkeley_1_75A();

		}
		else if ( version.equalsIgnoreCase(VERSION_1_75B)) {
			return getBerkeley_1_75B();
		} else if ( version.equalsIgnoreCase(VERSION_1_75)) {
			return getScop_1_75();
		} else {
			System.err.println("Unknown SCOP version " + version + " . Returning default");
			
			return scop;
		}
	}

	private static ScopDatabase getBerkeley_1_75A() {
		BerkeleyScopInstallation berkeley_1_75A = (BerkeleyScopInstallation) versionedScopDBs.get(VERSION_1_75A);
		
		if ( berkeley_1_75A == null) {
			berkeley_1_75A  = new BerkeleyScopInstallation();
			berkeley_1_75A.setScopVersion(VERSION_1_75A);				
			versionedScopDBs.put(VERSION_1_75A, berkeley_1_75A);
		}
		return berkeley_1_75A;
	}

	private static ScopDatabase getBerkeley_1_75B() {
		BerkeleyScopInstallation berkeley_1_75B = (BerkeleyScopInstallation) versionedScopDBs.get(VERSION_1_75B);
		
		if ( berkeley_1_75B == null) {
			berkeley_1_75B  = new BerkeleyScopInstallation();
			berkeley_1_75B.setScopVersion(VERSION_1_75B);				
			versionedScopDBs.put(VERSION_1_75B, berkeley_1_75B);
		}
		return berkeley_1_75B;
	}
}
