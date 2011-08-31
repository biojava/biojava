package org.biojava.bio.structure.align.util;

import org.biojava3.core.util.SoftHashMap;


/** provides a SoftHashMap singleton.
 * 
 * 
 * @Deprecated find better ways for caching or use a SoftHashMap directly
 */

public class CacheFactory  {

	@SuppressWarnings("rawtypes")
	private static SoftHashMap  cache  = new SoftHashMap ();
	
	// no public constructor;
	private CacheFactory(){
		
	}
	
	@SuppressWarnings("rawtypes")
	public static SoftHashMap getCache(){
		return cache;
	}
	
}
