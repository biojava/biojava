package org.biojava.bio.structure.scop;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;

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
		scop = s;
	}
		
		
		
}
