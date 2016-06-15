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
package org.biojava.nbio.structure.domain;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.HashMap;
import java.util.Map;


/** A class that provides all that is necessary to create a Serializable Cache
 *
 * @author Andreas Prlic
 *
 * @param <K> The key type of the cache
 * @param <V> the value type to be stored on the cache
 *
 * @since 3.0.3
 */
public class SerializableCache <K,V>{

	private static final Logger logger = LoggerFactory.getLogger(SerializableCache.class);

	protected String cacheFileName;
	protected Map<K,V> serializedCache ;


	/** set cacheFileName to null to disable caching
	 *
	 * @param cacheFileName
	 */
	public SerializableCache(String cacheFileName ) {
		this.cacheFileName = cacheFileName;

		if ( cacheFileName != null) {
			reloadFromFile();
		}

	}

	public boolean isCacheEnabled(){
		return ( serializedCache != null );
	}

	/** This will not cache null values.
	 *  Null means not cached yet.
	 *  If you want to cache "no data exists" use e.g. empty collections to represent this.
	 *
	 * @param name
	 * @param data
	 */
	public void cache(K name, V data) {

		if ( data == null){
			return;
		}
		if ( serializedCache != null){


			logger.debug("Caching {}  {}", name, data);

			serializedCache.put(name,data);


			// every 1000 objects we are writing to disk
			if ( serializedCache.keySet().size() % 1000 == 0 ) {

				flushCache();

			}

		}


	}

	public V get(K name) {
		if ( serializedCache == null)
			return null;
		return (serializedCache.get(name));
	}

	public void disableCache(){
		//flushCache();
		serializedCache = null;
	}

	public void enableCache(){
		reloadFromFile();
	}



	@SuppressWarnings("unchecked")
	public Map<K,V> reloadFromFile() {

		File f = getCacheFile();

		serializedCache = new HashMap<K,V>();

		// has never been cached here before
		if( ! f.exists()) {
			logger.info("Creating new cache " + f.getAbsolutePath());
			return serializedCache;
		}

		try{

			logger.debug("Reloading from cache " + f.getAbsolutePath());

			FileInputStream fis = new FileInputStream(f);
			ObjectInputStream ois = new ObjectInputStream(fis);
			serializedCache = (HashMap<K,V>) ois.readObject();
			ois.close();
		} catch (IOException e){
			// TODO shouldn't this be thrown forward?
			logger.error("Exception caught while reading serialized file",e);
			return null;
		} catch (ClassNotFoundException e) {
			logger.error("Exception caught while reading serialized file",e);
			return null;
		}

		//if ( debug )
		logger.info("Reloaded from cache: " + f.getName()+ " size: " + serializedCache.keySet().size() + " cached records.");
		return serializedCache;
	}

	private File getCacheFile() {
		AtomCache cache =new AtomCache();
		String path = cache.getCachePath();
		File f = new File(path + System.getProperty("file.separator") + cacheFileName);

		logger.debug(f.getAbsolutePath());
		return f;
	}

	public void flushCache(){
		if ( serializedCache == null)
			return;
		synchronized(serializedCache){

			File f = getCacheFile();

			try {

				FileOutputStream fos = new FileOutputStream(f);
				ObjectOutputStream oos = new ObjectOutputStream(fos);
				oos.writeObject(serializedCache);
				oos.close();
			} catch (IOException e){
				logger.error("Exception caught", e);
			}
		}
	}

}
