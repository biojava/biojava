package org.biojava.bio.structure.align.util;

import org.biojava3.core.util.SoftHashMap;


public class CacheFactory {

	private static SoftHashMap cache  = new SoftHashMap();
	
	// no public constructor;
	private CacheFactory(){
		
	}
	
	public static SoftHashMap getCache(){
		return cache;
	}
	
}
