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
 * Created on Oct 1, 2009
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.core.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;

/** 
 * Provides a cache for storing multiple small files in memory. Can be used to e.g cache gzip compressed PDB files 
 * for avoiding disk IO bottlenecks.
 * Note this is just a wrapper for the singleton cache.
 * 
 * @author Andreas Prlic.
 *
 */
public class FlatFileCache {

	private final static Logger logger = LoggerFactory.getLogger(FlatFileCache.class);

	/**
	 * The cache singleton.
	 */
	private static SoftHashMap<String, byte[]> cache = new SoftHashMap<String, byte[]>(0);


	// no public constructor;
	private FlatFileCache(){

	}


	public  static void addToCache(String key, File fileToCache){
		//logger.debug("storing " + key + " on file cache (cache size: " + cache.size() + ")");
		try {
			InputStream is = new FileInputStream(fileToCache);
			// Get the size of the file
			long length = fileToCache.length();

			// You cannot create an array using a long type.
			// It needs to be an int type.
			// Before converting to an int type, check
			// to ensure that file is not larger than Integer.MAX_VALUE.
			if (length > Integer.MAX_VALUE) {
				// File is too large
			}

			// Create the byte array to hold the data
			byte[] bytes = new byte[(int)length];

			// Read in the bytes
			int offset = 0;
			int numRead = 0;
			while (offset < bytes.length
					&& (numRead=is.read(bytes, offset, bytes.length-offset)) >= 0) {
				offset += numRead;
			}

			// Ensure all the bytes have been read in
			if (offset < bytes.length) {
				is.close();
				throw new IOException("Could not completely read file "+fileToCache.getName());
			}

			// Close the input stream and return bytes
			is.close();

			cache.put(key,bytes);

		} catch (Exception e){
			logger.error("Error adding to cache! " + e.getMessage(), e);
		}
	}

	public  static InputStream getInputStream(String key){
		//logger.debug("returning " + key + " from file cache (cache size: " + cache.size() + ")");
		byte[] bytes = cache.get(key);
		if ( bytes == null)
			return null;

		return new ByteArrayInputStream(bytes);

	}

	public static int size() {
		if ( cache != null)
			return cache.size();
		else
			return -1;
	}

	public static void clear(){
	   cache.clear();
	}

	

}
