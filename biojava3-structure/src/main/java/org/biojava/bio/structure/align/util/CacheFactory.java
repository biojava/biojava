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
